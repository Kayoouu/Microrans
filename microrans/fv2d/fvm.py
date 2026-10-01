"""Opérateurs volumes finis colocalisés, 2D (polygones) et 3D (polyèdres : mesh3d), écrits
face par face : la dimension n'intervient que par le nombre de composantes des vecteurs.

Conventions (proches d'OpenFOAM) :
- équation d'une cellule P :  Σ_f F_f φ_f − Σ_f Γ_f (∇φ)_f·S_f = S_exp V − D_imp V φ_P (+ temps)
- matrice creuse = diagonale + coefficient « upper » (ligne owner, colonne neighbour)
  + coefficient « lower » (ligne neighbour, colonne owner) pour chaque face interne ;
- condition aux limites d'une face frontière : φ_b = α φ_P + β et
  ∂φ/∂n|_b = γ φ_P + δ (valueInternalCoeffs / valueBoundaryCoeffs / gradient… d'OpenFOAM).

Axisymétrique (x = axe, y = rayon r ≥ 0) : volumes et surfaces sont ceux d'un secteur d'un
radian, V = ∫ r dA = A·r_centre et S_f = S_f,plan·r_face (exact pour des arêtes droites) —
c'est la formulation « wedge » d'OpenFOAM avec un angle infinitésimal. Les faces sur l'axe
ont une surface nulle. Le gradient de Green-Gauss inclut la contribution des faces latérales
du secteur (−φ_P Σ S_f), sans laquelle le gradient d'un champ constant serait non nul.
"""
from __future__ import annotations

import numpy as np

from ..linalg import LinearSolver, array_module
from ..mesh2d.mesh import Mesh2D


class FVM:
    def __init__(self, mesh: Mesh2D, backend=None, axisymmetric: bool = False, numba=False,
                 threads=1):
        from ..backend import get_backend
        m = mesh
        self.mesh = m
        self.dim = int(getattr(m, "dim", m.Sf.shape[1]))
        self.axisymmetric = bool(axisymmetric)
        if self.axisymmetric and self.dim != 2:
            raise ValueError("Axisymétrique : maillage 2D seulement (le maillage est en 3D).")
        self.backend = be = get_backend(backend)
        self.xp = xp = be.xp
        self.nc, self.ni, self.nf = m.n_cells, m.n_internal, m.n_faces
        self.nb = self.nf - self.ni
        ni = self.ni
        P, N = m.owner[:ni], m.neighbour
        C = m.cell_centers
        rows = np.concatenate([np.arange(self.nc), P, N])
        cols = np.concatenate([np.arange(self.nc), N, P])
        # structure CSR figée : chaque matrice = une somme pondérée vectorisée (bincount)
        key = rows.astype(np.int64) * self.nc + cols
        uk, csr_pos = np.unique(key, return_inverse=True)
        self._csr_nnz = len(uk)
        indptr = np.concatenate([[0], np.cumsum(np.bincount(uk // self.nc, minlength=self.nc))])
        nonorth = np.linalg.norm(m.nonorth_vec, axis=1) / np.maximum(m.magSf[:ni], 1e-300)
        self.orthogonal = bool(np.all(nonorth < 1e-8))
        Sf, magSf, V = m.Sf, m.magSf, m.cell_volumes
        g, kvec = m.orth_coeff, m.nonorth_vec
        if self.axisymmetric:
            rc = m.cell_centers[:, 1]
            scale = max(float(np.ptp(m.points[:, 1])), 1e-300)
            if np.any(rc <= 0.0) or np.any(m.points[:, 1] < -1e-9 * scale):
                raise ValueError("Axisymétrique : le domaine doit être dans le demi-plan y ≥ 0 "
                                 "(x = axe de révolution, y = rayon). Couper le maillage à "
                                 "l'axe ([mesh] cut_axis = true) ou le décaler.")
            rf = np.maximum(m.face_centers[:, 1], 0.0)
            Sf, magSf, V = Sf * rf[:, None], magSf * rf, V * rc
            g, kvec = g * rf[:ni], kvec * rf[:ni, None]
        self.total_area = float(np.sum(magSf))
        # solveurs linéaires (hiérarchie AMG bâtie sur le graphe CPU, poids |S|²/(d·S))
        # LU creuse : remplissage bien plus fort en 3D (mesuré, laplacien, rtol 1e-6 : LU
        # 0.12 s contre AMG 0.004 s à 4 096 cellules, 7 s contre 0.04 s à 32 768)
        self.lin = LinearSolver(self.nc, P, N, g, direct_max=20000 if self.dim == 2 else 2000)
        # tableaux transférés sur le matériel de calcul (CPU : aucune copie)
        A = be.asarray
        self.P, self.N, self.Pb = A(P), A(N), A(m.owner[ni:])
        self.w, self.g, self.kvec = A(m.weights), A(g), A(kvec)
        self.Si, self.Sb = A(Sf[:ni]), A(Sf[ni:])
        self.magSb, self.nb_hat = A(magSf[ni:]), A(m.nf[ni:])
        self.V, self.dperp = A(V), A(m.d_perp_b)
        if self.axisymmetric:
            self.radius = A(m.cell_centers[:, 1])
            # Σ_f S_f / V (= ê_r / r_P) : faces latérales du secteur, pour le gradient
            sx = np.bincount(P, Sf[:ni, 0], self.nc) - np.bincount(N, Sf[:ni, 0], self.nc)
            sy = np.bincount(P, Sf[:ni, 1], self.nc) - np.bincount(N, Sf[:ni, 1], self.nc)
            if self.nb:
                sx = sx + np.bincount(m.owner[ni:], Sf[ni:, 0], self.nc)
                sy = sy + np.bincount(m.owner[ni:], Sf[ni:, 1], self.nc)
            self._side = A(np.column_stack([sx, sy]) / V[:, None])
        self.rP = A(m.face_centers[:ni] - C[P])
        self.rN = A(m.face_centers[:ni] - (C[N] - m.shift))
        self.rows, self.cols = A(rows), A(cols)
        self._csr_pos = A(csr_pos.ravel())
        self._csr_indices = A((uk % self.nc).astype(np.int32))
        self._csr_indptr = A(indptr.astype(np.int32))
        # noyaux multi-cœurs (Numba, facultatif, CPU seulement) : structures « cellule →
        # faces » calculées une fois (voir kernels.py)
        from . import kernels
        self.fast = be.name == "cpu" and kernels.wanted(numba)
        if self.fast:
            self.threads = kernels.set_threads(threads)
            self._k = kernels
            nb = self.nb
            cells = np.concatenate([P, N, m.owner[ni:]])
            fac = np.concatenate([np.arange(ni), np.arange(ni), ni + np.arange(nb)])
            sgn = np.concatenate([np.ones(ni), -np.ones(ni), np.ones(nb)])
            ptr, order = kernels.segments(cells, self.nc)
            self._k_all = (ptr, fac[order].astype(np.int64), sgn[order])
            self._k_P = kernels.segments(P, self.nc)
            self._k_N = kernels.segments(N, self.nc)
            self._k_B = kernels.segments(m.owner[ni:], self.nc)
            self._k_M = kernels.segments(csr_pos.ravel(), self._csr_nnz)
            self.Si, self.Sb = np.ascontiguousarray(self.Si), np.ascontiguousarray(self.Sb)
            self.w = np.ascontiguousarray(self.w, dtype=float)

    def zeros(self, *shape):
        return self.xp.zeros(shape if len(shape) > 1 else shape[0])

    def ones(self, n):
        return self.xp.ones(n)

    # ------------------------------------------------------------- opérateurs explicites
    def interp(self, phi):
        """Interpolation linéaire aux faces internes (scalaire ou tableau (nc, k))."""
        if self.fast and phi.dtype == np.float64:
            if phi.ndim == 1:
                out = np.empty(self.ni)
                self._k.interp1(self.w, np.ascontiguousarray(phi), self.P, self.N, out)
                return out
            p2 = np.ascontiguousarray(phi.reshape(self.nc, -1))
            out = np.empty((self.ni, p2.shape[1]))
            self._k.interp2(self.w, p2, self.P, self.N, out)
            return out.reshape((self.ni,) + phi.shape[1:])
        w = self.w.reshape((-1,) + (1,) * (phi.ndim - 1))
        return w * phi[self.P] + (1.0 - w) * phi[self.N]

    def sum_faces(self, fi, fb):
        """Σ des flux sortants par cellule (face interne : + owner, − neighbour)."""
        if self.fast:
            out = np.empty(self.nc)
            ptr, fac, sgn = self._k_all
            self._k.facesum(ptr, fac, sgn, np.ascontiguousarray(fi, dtype=float),
                            np.ascontiguousarray(fb, dtype=float), self.ni, out)
            return out
        bc = self.xp.bincount
        out = bc(self.P, weights=fi, minlength=self.nc) - bc(self.N, weights=fi,
                                                             minlength=self.nc)
        if self.nb:
            out = out + bc(self.Pb, weights=fb, minlength=self.nc)
        return out

    def grad(self, phi, phi_b):
        """Gradient de Green-Gauss (interpolation linéaire, valeurs frontières imposées)."""
        if self.fast and phi.dtype == np.float64:
            g = np.empty((self.nc, self.dim))
            ptr, fac, sgn = self._k_all
            self._k.green_gauss(ptr, fac, sgn, self.w, self.P, self.N,
                                np.ascontiguousarray(phi), np.ascontiguousarray(phi_b, dtype=float),
                                self.Si, self.Sb, self.V, self.ni, g)
            if self.axisymmetric:
                g = g - phi[:, None] * self._side
            return g
        pf = self.interp(phi)
        g = self.xp.stack([self.sum_faces(pf * self.Si[:, k], phi_b * self.Sb[:, k])
                           for k in range(self.dim)], axis=1) / self.V[:, None]
        if self.axisymmetric:
            g = g - phi[:, None] * self._side
        return g

    def limit_grad(self, phi, grad, phi_b):
        """Limiteur de Barth-Jespersen (« cellLimited Gauss linear 1 » d'OpenFOAM) : le
        gradient est réduit pour que les valeurs reconstruites aux faces restent entre le
        minimum et le maximum de la cellule et de ses voisines (valeurs frontières
        comprises). Uniquement des lectures indexées : même code sur tous les backends."""
        xp = self.xp
        if getattr(self, "_lim", None) is None:
            self._lim = self._limiter_stencil()
        idx, R = self._lim
        ext = xp.concatenate([phi, phi_b])[idx]                  # (nc, k), cellule incluse
        dmax = xp.max(ext, axis=1) - phi
        dmin = xp.min(ext, axis=1) - phi
        d = xp.sum(grad[:, None, :] * R, axis=2)                 # variation cellule → face
        pos, neg = d > 0.0, d < 0.0
        r = xp.where(pos, dmax[:, None] / xp.where(pos, d, 1.0),
                     xp.where(neg, dmin[:, None] / xp.where(neg, d, -1.0), 1.0))
        return grad * xp.minimum(xp.min(r, axis=1), 1.0)[:, None]

    def _limiter_stencil(self):
        m, ni, nc = self.mesh, self.ni, self.nc
        P, N, Pb = m.owner[:ni], m.neighbour, m.owner[ni:]
        C = m.cell_centers
        cells = np.concatenate([P, N, Pb])
        other = np.concatenate([N, P, nc + np.arange(self.nb)])
        R = np.concatenate([m.face_centers[:ni] - C[P],
                            m.face_centers[:ni] - (C[N] - m.shift),
                            m.face_centers[ni:] - C[Pb]])
        order = np.argsort(cells, kind="stable")
        cs = cells[order]
        counts = np.bincount(cells, minlength=nc)
        start = np.concatenate([[0], np.cumsum(counts)[:-1]])
        pos = np.arange(len(cs)) - start[cs]
        k = int(counts.max()) + 1                                # + la cellule elle-même
        idx = np.tile(np.arange(nc)[:, None], (1, k))
        Rk = np.zeros((nc, k, self.dim))
        idx[cs, pos] = other[order]
        Rk[cs, pos] = R[order]
        A = self.backend.asarray
        return A(idx), A(Rk)

    def div(self, flux_i, flux_b):
        return self.sum_faces(flux_i, flux_b)

    # ------------------------------------------------------------- assemblage
    def _sum(self, idx, w):
        if self.fast:
            seg = (self._k_P if idx is self.P else self._k_N if idx is self.N
                   else self._k_B if idx is self.Pb else None)
            if seg is not None:
                out = np.empty(self.nc)
                self._k.segsum(seg[0], seg[1], np.ascontiguousarray(w, dtype=float), out)
                return out
        return self.xp.bincount(idx, weights=w, minlength=self.nc)

    def _sum_b(self, w):
        return self._sum(self.Pb, w) if self.nb else self.xp.zeros(self.nc)

    def nonorth_flux(self, gam_i, grad_phi, phi, limit: float | None = 0.5):
        """Partie non orthogonale du flux diffusif Γ (∇φ)_f·k, limitée comme `limited ψ` d'OpenFOAM."""
        xp = self.xp
        corr = gam_i * xp.sum(self.interp(grad_phi) * self.kvec, axis=1)
        if limit is not None and limit < 1.0:
            orth = gam_i * self.g * xp.abs(phi[self.N] - phi[self.P])
            cap = limit / (1.0 - limit) * orth
            corr = xp.clip(corr, -cap, cap)
        return corr

    def assemble(self, F_i, F_b, gam_i, gam_b, bc, grad_phi=None, scheme="upwind",
                 bounded=False, phi=None, nonorth_limit: float | None = 0.5, grad_conv=None):
        """Convection (flux F) + diffusion (Γ) d'un scalaire. Retourne (diag, upper, lower, rhs).

        grad_conv : gradient utilisé par la correction linearUpwind (par défaut grad_phi ;
        gradient limité pour linearUpwindLimited), grad_phi servant à la correction non
        orthogonale de la diffusion."""
        xp = self.xp
        P, N = self.P, self.N
        alpha, beta, gamma, delta = bc
        Fp = xp.maximum(F_i, 0.0)
        Fm = xp.minimum(F_i, 0.0)
        dcoef = gam_i * self.g
        diag = (self._sum(P, Fp + dcoef) + self._sum(N, -Fm + dcoef)
                + self._sum_b(F_b * alpha - gam_b * self.magSb * gamma))
        upper = Fm - dcoef
        lower = -Fp - dcoef
        rhs = self._sum_b(-F_b * beta + gam_b * self.magSb * delta)
        if bounded:
            diag = diag - self.sum_faces(F_i, F_b)
        if grad_phi is not None:
            if scheme in ("linearUpwind", "linearUpwindLimited"):
                gc = grad_phi if grad_conv is None else grad_conv
                corr = xp.where(F_i >= 0.0, xp.sum(gc[P] * self.rP, axis=1),
                                xp.sum(gc[N] * self.rN, axis=1))
                fc = F_i * corr
                rhs = rhs - (self._sum(P, fc) - self._sum(N, fc))
            if not self.orthogonal and bool(xp.any(gam_i)):
                lim = nonorth_limit if phi is not None else None
                nf = self.nonorth_flux(gam_i, grad_phi, phi, lim)
                rhs = rhs + self._sum(P, nf) - self._sum(N, nf)
        return diag, upper, lower, rhs

    def matrix(self, diag, upper, lower):
        if self.fast:
            data = np.empty(self._csr_nnz)
            vals = np.concatenate([diag, upper, lower]).astype(float, copy=False)
            self._k.segsum(self._k_M[0], self._k_M[1], vals, data)
            return self.backend.sparse.csr_matrix((data, self._csr_indices, self._csr_indptr),
                                                  shape=(self.nc, self.nc))
        data = self.xp.bincount(self._csr_pos, weights=self.xp.concatenate([diag, upper, lower]),
                                minlength=self._csr_nnz)
        return self.backend.sparse.csr_matrix((data, self._csr_indices, self._csr_indptr),
                                              shape=(self.nc, self.nc))


# ----------------------------------------------------------------------- solveurs linéaires
def normalized_residual(A, x, b, scale: float = 0.0) -> float:
    """Résidu normalisé à la manière d'OpenFOAM (indépendant de l'échelle de x).

    `scale` (ordre de grandeur attendu du champ) ajoute un plancher au dénominateur pour
    qu'un champ identiquement nul (ex. v = 0) ne donne pas « bruit / bruit ≈ 1 ».
    """
    xp = array_module(x)
    Ax = A @ x
    Axbar = A @ (xp.ones_like(x) * x.mean())
    floor = 1e-8 * scale * float(xp.sum(xp.abs(A.diagonal())))
    norm = float(xp.sum(xp.abs(Ax - Axbar) + xp.abs(b - Axbar))) + floor + 1e-300
    return float(xp.sum(xp.abs(b - Ax))) / norm


def solve_linear(A, b, x0=None, method: str = "auto", rtol: float = 1e-6,
                 maxiter: int = 1000, symmetric: bool = False):
    """Résolution ponctuelle (sans hiérarchie AMG réutilisée) ; voir `linalg.LinearSolver`."""
    A = A.tocsr()
    off = A.tocoo()
    sel = off.row < off.col
    return LinearSolver(A.shape[0], off.row[sel], off.col[sel], np.abs(off.data[sel])).solve(
        A, b, x0, method, rtol, maxiter, symmetric)
