"""Opérateurs volumes finis 2D colocalisés (maillages non structurés polygonaux).

Conventions (proches d'OpenFOAM) :
- équation d'une cellule P :  Σ_f F_f φ_f − Σ_f Γ_f (∇φ)_f·S_f = S_exp V − D_imp V φ_P (+ temps)
- matrice creuse = diagonale + coefficient « upper » (ligne owner, colonne neighbour)
  + coefficient « lower » (ligne neighbour, colonne owner) pour chaque face interne ;
- condition aux limites d'une face frontière : φ_b = α φ_P + β et
  ∂φ/∂n|_b = γ φ_P + δ (valueInternalCoeffs / valueBoundaryCoeffs / gradient… d'OpenFOAM).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from ..linalg import LinearSolver
from ..mesh2d.mesh import Mesh2D


class FVM:
    def __init__(self, mesh: Mesh2D):
        m = mesh
        self.mesh = m
        self.nc, self.ni, self.nf = m.n_cells, m.n_internal, m.n_faces
        self.nb = self.nf - self.ni
        ni = self.ni
        self.P = m.owner[:ni]
        self.N = m.neighbour
        self.Pb = m.owner[ni:]
        self.w = m.weights
        self.g = m.orth_coeff
        self.kvec = m.nonorth_vec
        self.Si = m.Sf[:ni]
        self.Sb = m.Sf[ni:]
        self.magSb = m.magSf[ni:]
        self.nb_hat = m.nf[ni:]
        self.V = m.cell_volumes
        self.dperp = m.d_perp_b
        C = m.cell_centers
        self.rP = m.face_centers[:ni] - C[self.P]
        self.rN = m.face_centers[:ni] - (C[self.N] - m.shift)
        self.rows = np.concatenate([np.arange(self.nc), self.P, self.N])
        self.cols = np.concatenate([np.arange(self.nc), self.N, self.P])
        # structure CSR figée : chaque matrice = une somme pondérée vectorisée (bincount)
        key = self.rows.astype(np.int64) * self.nc + self.cols
        uk, self._csr_pos = np.unique(key, return_inverse=True)
        self._csr_nnz = len(uk)
        self._csr_indices = (uk % self.nc).astype(np.int32)
        self._csr_indptr = np.concatenate(
            [[0], np.cumsum(np.bincount(uk // self.nc, minlength=self.nc))]).astype(np.int32)
        # solveurs linéaires (hiérarchie AMG bâtie sur le graphe, poids |S|²/(d·S))
        self.lin = LinearSolver(self.nc, self.P, self.N, self.g)
        nonorth = np.linalg.norm(self.kvec, axis=1) / np.maximum(m.magSf[:ni], 1e-300)
        self.orthogonal = bool(np.all(nonorth < 1e-8))

    # ------------------------------------------------------------- opérateurs explicites
    def interp(self, phi):
        """Interpolation linéaire aux faces internes (scalaire ou tableau (nc, k))."""
        w = self.w.reshape((-1,) + (1,) * (np.ndim(phi) - 1))
        return w * phi[self.P] + (1.0 - w) * phi[self.N]

    def sum_faces(self, fi, fb):
        """Σ des flux sortants par cellule (face interne : + owner, − neighbour)."""
        return (np.bincount(self.P, fi, self.nc) - np.bincount(self.N, fi, self.nc)
                + np.bincount(self.Pb, fb, self.nc))

    def grad(self, phi, phi_b):
        """Gradient de Green-Gauss (interpolation linéaire, valeurs frontières imposées)."""
        pf = self.interp(phi)
        gx = self.sum_faces(pf * self.Si[:, 0], phi_b * self.Sb[:, 0])
        gy = self.sum_faces(pf * self.Si[:, 1], phi_b * self.Sb[:, 1])
        return np.column_stack([gx, gy]) / self.V[:, None]

    def div(self, flux_i, flux_b):
        return self.sum_faces(flux_i, flux_b)

    # ------------------------------------------------------------- assemblage
    def nonorth_flux(self, gam_i, grad_phi, phi, limit: float | None = 0.5):
        """Partie non orthogonale du flux diffusif Γ (∇φ)_f·k, limitée comme `limited ψ` d'OpenFOAM."""
        corr = gam_i * np.sum(self.interp(grad_phi) * self.kvec, axis=1)
        if limit is not None and limit < 1.0:
            orth = gam_i * self.g * np.abs(phi[self.N] - phi[self.P])
            cap = limit / (1.0 - limit) * orth
            corr = np.clip(corr, -cap, cap)
        return corr

    def assemble(self, F_i, F_b, gam_i, gam_b, bc, grad_phi=None, scheme="upwind",
                 bounded=False, phi=None, nonorth_limit: float | None = 0.5):
        """Convection (flux F) + diffusion (Γ) d'un scalaire. Retourne (diag, upper, lower, rhs)."""
        P, N, Pb, nc = self.P, self.N, self.Pb, self.nc
        alpha, beta, gamma, delta = bc
        Fp = np.maximum(F_i, 0.0)
        Fm = np.minimum(F_i, 0.0)
        dcoef = gam_i * self.g
        diag = (np.bincount(P, Fp + dcoef, nc) + np.bincount(N, -Fm + dcoef, nc)
                + np.bincount(Pb, F_b * alpha - gam_b * self.magSb * gamma, nc))
        upper = Fm - dcoef
        lower = -Fp - dcoef
        rhs = np.bincount(Pb, -F_b * beta + gam_b * self.magSb * delta, nc)
        if bounded:
            diag -= self.sum_faces(F_i, F_b)
        if grad_phi is not None:
            if scheme == "linearUpwind":
                corr = np.where(F_i >= 0.0, np.sum(grad_phi[P] * self.rP, axis=1),
                                np.sum(grad_phi[N] * self.rN, axis=1))
                fc = F_i * corr
                rhs -= np.bincount(P, fc, nc) - np.bincount(N, fc, nc)
            if not self.orthogonal and np.any(gam_i):
                lim = nonorth_limit if phi is not None else None
                nf = self.nonorth_flux(gam_i, grad_phi, phi, lim)
                rhs += np.bincount(P, nf, nc) - np.bincount(N, nf, nc)
        return diag, upper, lower, rhs

    def matrix(self, diag, upper, lower):
        data = np.bincount(self._csr_pos, np.concatenate([diag, upper, lower]), self._csr_nnz)
        return sp.csr_matrix((data, self._csr_indices, self._csr_indptr),
                             shape=(self.nc, self.nc))


# ----------------------------------------------------------------------- solveurs linéaires
def normalized_residual(A, x, b, scale: float = 0.0) -> float:
    """Résidu normalisé à la manière d'OpenFOAM (indépendant de l'échelle de x).

    `scale` (ordre de grandeur attendu du champ) ajoute un plancher au dénominateur pour
    qu'un champ identiquement nul (ex. v = 0) ne donne pas « bruit / bruit ≈ 1 ».
    """
    Ax = A @ x
    xbar = np.full_like(x, x.mean())
    Axbar = A @ xbar
    floor = 1e-8 * scale * np.sum(np.abs(A.diagonal()))
    norm = np.sum(np.abs(Ax - Axbar) + np.abs(b - Axbar)) + floor + 1e-300
    return float(np.sum(np.abs(b - Ax)) / norm)


def solve_linear(A, b, x0=None, method: str = "auto", rtol: float = 1e-6,
                 maxiter: int = 1000, symmetric: bool = False):
    """Résolution ponctuelle (sans hiérarchie AMG réutilisée) ; voir `linalg.LinearSolver`."""
    A = A.tocsr()
    off = A.tocoo()
    sel = off.row < off.col
    return LinearSolver(A.shape[0], off.row[sel], off.col[sel], np.abs(off.data[sel])).solve(
        A, b, x0, method, rtol, maxiter, symmetric)
