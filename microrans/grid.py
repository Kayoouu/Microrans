"""Maillage 1D d'un canal plan : y ∈ [0, 2h], parois en y = 0 et y = 2h.

Les inconnues sont stockées aux noeuds (volumes finis centrés sur les noeuds).
Les noeuds 0 et N sont sur les parois (conditions de Dirichlet).
"""
from __future__ import annotations

import numpy as np


class Grid:
    """Maillage 1D non uniforme avec les grandeurs géométriques dérivées."""

    is_1d = True

    def __init__(self, y):
        y = np.asarray(y, dtype=float)
        if y.ndim != 1 or y.size < 5:
            raise ValueError("Il faut un tableau 1D d'au moins 5 noeuds.")
        if np.any(np.diff(y) <= 0.0):
            raise ValueError("Les noeuds doivent être strictement croissants.")
        self.y = y
        self.n = y.size
        self.h = 0.5 * (y[-1] - y[0])
        # Pas entre noeuds (taille n-1) : la face i+1/2 est entre les noeuds i et i+1
        self.dy = np.diff(y)
        # Volumes de contrôle duaux (demi-volumes aux parois, non utilisés)
        vol = np.empty_like(y)
        vol[1:-1] = 0.5 * (y[2:] - y[:-2])
        vol[0] = 0.5 * self.dy[0]
        vol[-1] = 0.5 * self.dy[-1]
        self.volumes = vol
        self.wall_distance = np.minimum(y - y[0], y[-1] - y)

    @property
    def first_cell_height(self) -> float:
        return float(min(self.dy[0], self.dy[-1]))

    def __repr__(self) -> str:
        return (f"Grid(n={self.n}, h={self.h:g}, dy_min={self.dy.min():.3e}, "
                f"dy_max={self.dy.max():.3e})")


def tanh_grid(n_cells: int, gamma: float, h: float = 1.0) -> Grid:
    """Maillage symétrique resserré aux deux parois par une loi en tanh.

    y_j = h (1 + tanh(γ ξ_j) / tanh(γ)),  ξ_j = -1 + 2 j / N.
    γ → 0 donne un maillage uniforme.
    """
    if n_cells < 4:
        raise ValueError("n_cells doit être >= 4.")
    xi = np.linspace(-1.0, 1.0, n_cells + 1)
    if gamma < 1e-6:
        y = h * (1.0 + xi)
    else:
        y = h * (1.0 + np.tanh(gamma * xi) / np.tanh(gamma))
    y[0], y[-1] = 0.0, 2.0 * h
    return Grid(y)


def channel_grid(n_cells: int, re_tau: float, y1_plus: float = 0.5,
                 h: float = 1.0) -> Grid:
    """Maillage en tanh dont la première maille vaut `y1_plus` en unités de paroi.

    L'unité de paroi est ν/u_τ = h/Re_τ (u_τ nominale). Si un maillage uniforme
    donne déjà y1+ <= y1_plus, le maillage uniforme est renvoyé.
    """
    target = y1_plus * h / re_tau

    xi1 = -1.0 + 2.0 / n_cells

    def first_height(gamma: float) -> float:
        if gamma < 1e-6:
            return 2.0 * h / n_cells
        return h * (1.0 + np.tanh(gamma * xi1) / np.tanh(gamma))

    if first_height(0.0) <= target:
        return tanh_grid(n_cells, 0.0, h)
    lo, hi = 1e-6, 10.0
    if first_height(hi) > target:
        raise ValueError(
            f"Impossible d'atteindre y1+={y1_plus} avec {n_cells} mailles : "
            "augmenter n_cells.")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if first_height(mid) > target:
            lo = mid
        else:
            hi = mid
    return tanh_grid(n_cells, hi, h)
