"""Opérateurs volumes finis 2D colocalisés (maillages non structurés polygonaux).

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
    def __init__(self, mesh: Mesh2D, backend=None, axisymmetric: bool = False):
        from ..backend import get_backend
        m = mesh
        self.mesh = m
        self.axisymmetric = bool(axisymmetric)
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
        self.lin = LinearSolver(self.nc, P, N, g)
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

    def zeros(self, *shape):
        return self.xp.zeros(shape if len(shape) > 1 else shape[0])

    def ones(self, n):
        return self.xp.ones(n)

    # ------------------------------------------------------------- opérateurs explicites
    def interp(self, phi):
        """Interpolation linéaire aux faces internes (scalaire ou tableau (nc, k))."""
        w = self.w.reshape((-1,) + (1,) * (phi.ndim - 1))
        return w * phi[self.P] + (1.0 - w) * phi[self.N]

    def sum_faces(self, fi, fb):
        """Σ des flux sortants par cellule (face interne : + owner, − neighbour)."""
        bc = self.xp.bincount
        out = bc(self.P, weights=fi, minlength=self.nc) - bc(self.N, weights=fi,
                                                             minlength=self.nc)
        if self.nb:
            out = out + bc(self.Pb, weights=fb, minlength=self.nc)
        return out

    def grad(self, phi, phi_b):
        """Gradient de Green-Gauss (interpolation linéaire, valeurs frontières imposées)."""
        pf = self.interp(phi)
        gx = self.sum_faces(pf * self.Si[:, 0], phi_b * self.Sb[:, 0])
        gy = self.sum_faces(pf * self.Si[:, 1], phi_b * self.Sb[:, 1])
        g = self.xp.stack([gx, gy], axis=1) / self.V[:, None]
        if self.axisymmetric:
            g = g - phi[:, None] * self._side
        return g

    def div(self, flux_i, flux_b):
        return self.sum_faces(flux_i, flux_b)

    # ------------------------------------------------------------- assemblage
    def _sum(self, idx, w):
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
                 bounded=False, phi=None, nonorth_limit: float | None = 0.5):
        """Convection (flux F) + diffusion (Γ) d'un scalaire. Retourne (diag, upper, lower, rhs)."""
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
            if scheme == "linearUpwind":
                corr = xp.where(F_i >= 0.0, xp.sum(grad_phi[P] * self.rP, axis=1),
                                xp.sum(grad_phi[N] * self.rN, axis=1))
                fc = F_i * corr
                rhs = rhs - (self._sum(P, fc) - self._sum(N, fc))
            if not self.orthogonal and bool(xp.any(gam_i)):
                lim = nonorth_limit if phi is not None else None
                nf = self.nonorth_flux(gam_i, grad_phi, phi, lim)
                rhs = rhs + self._sum(P, nf) - self._sum(N, nf)
        return diag, upper, lower, rhs

    def matrix(self, diag, upper, lower):
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
