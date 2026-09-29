"""Opérateurs discrets 1D sur maillage non uniforme et solveur tridiagonal implicite."""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from .grid import Grid


def ddy(grid: Grid, f: np.ndarray) -> np.ndarray:
    """Dérivée première aux noeuds, ordre 2 (centrée à l'intérieur, décentrée aux parois)."""
    y = grid.y
    out = np.empty_like(f)
    h1 = y[1:-1] - y[:-2]
    h2 = y[2:] - y[1:-1]
    out[1:-1] = (-h2 / (h1 * (h1 + h2)) * f[:-2]
                 + (h2 - h1) / (h1 * h2) * f[1:-1]
                 + h1 / (h2 * (h1 + h2)) * f[2:])
    a, b = y[1] - y[0], y[2] - y[1]
    out[0] = (-(2 * a + b) / (a * (a + b)) * f[0] + (a + b) / (a * b) * f[1]
              - a / (b * (a + b)) * f[2])
    a, b = y[-1] - y[-2], y[-2] - y[-3]
    out[-1] = ((2 * a + b) / (a * (a + b)) * f[-1] - (a + b) / (a * b) * f[-2]
               + a / (b * (a + b)) * f[-3])
    return out


def d2dy2(grid: Grid, f: np.ndarray) -> np.ndarray:
    """Dérivée seconde aux noeuds intérieurs (ordre 2) ; valeurs pariétales recopiées."""
    y = grid.y
    out = np.empty_like(f)
    h1 = y[1:-1] - y[:-2]
    h2 = y[2:] - y[1:-1]
    out[1:-1] = 2.0 * (h2 * f[:-2] - (h1 + h2) * f[1:-1] + h1 * f[2:]) / (h1 * h2 * (h1 + h2))
    out[0], out[-1] = out[1], out[-2]
    return out


def integrate(grid: Grid, f: np.ndarray) -> float:
    """Intégrale sur [0, 2h] par la méthode des trapèzes."""
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * grid.dy))


def solve_transport(grid: Grid, gamma: np.ndarray, diag: np.ndarray | float,
                    rhs: np.ndarray, wall_values=(0.0, 0.0)) -> np.ndarray:
    """Résout  diag·φ − d/dy(Γ dφ/dy) = rhs  aux noeuds intérieurs, Dirichlet aux parois.

    Discrétisation volumes finis centrés noeuds : flux aux faces
    Γ_{i+1/2} (φ_{i+1} − φ_i) / Δy_{i+1/2}, avec Γ_{i+1/2} moyenne arithmétique.
    """
    n = grid.n
    gamma_f = 0.5 * (gamma[:-1] + gamma[1:])
    cond = gamma_f / grid.dy                   # conductance de chaque face
    vol = grid.volumes[1:-1]
    a_w = cond[:-1] / vol
    a_e = cond[1:] / vol

    ab = np.zeros((3, n))
    ab[1, 1:-1] = np.broadcast_to(diag, (n,))[1:-1] + a_w + a_e
    ab[1, 0] = ab[1, -1] = 1.0
    ab[0, 2:] = -a_e        # sur-diagonale : A[i, i+1]
    ab[2, :-2] = -a_w       # sous-diagonale : A[i, i-1]

    b = np.array(rhs, dtype=float, copy=True)
    b[0], b[-1] = wall_values
    phi = solve_banded((1, 1), ab, b)
    # Le pivotage LAPACK peut laisser ~1e-14 aux parois : on impose la valeur exacte
    # (sinon √k pariétal ≈ 1e-7 pollue le terme D = 2ν(∂√k/∂y)² du k-ε).
    phi[0], phi[-1] = wall_values
    return phi


def diffusion(grid: Grid, gamma: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """d/dy(Γ dφ/dy) explicite, même discrétisation que `solve_transport` (0 aux parois)."""
    gamma_f = 0.5 * (gamma[:-1] + gamma[1:])
    flux = gamma_f * np.diff(phi) / grid.dy
    out = np.zeros_like(phi)
    out[1:-1] = (flux[1:] - flux[:-1]) / grid.volumes[1:-1]
    return out
