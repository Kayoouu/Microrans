"""Champ moyen U(y), ses dérivées, diagnostics pariétaux et champ initial."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid import Grid
from .numerics import d2dy2, ddy, integrate


@dataclass
class FlowField:
    """Vitesse moyenne et dérivées aux noeuds (canal : seul dU/dy est non nul)."""
    U: np.ndarray
    dudy: np.ndarray
    d2udy2: np.ndarray

    @classmethod
    def from_velocity(cls, grid: Grid, U: np.ndarray) -> "FlowField":
        return cls(U=U, dudy=ddy(grid, U), d2udy2=d2dy2(grid, U))

    @property
    def strain(self) -> np.ndarray:
        """S = sqrt(2 S_ij S_ij) = |dU/dy| (égal au module du rotationnel en canal)."""
        return np.abs(self.dudy)


def wall_shear(grid: Grid, nu: float, U: np.ndarray) -> tuple[float, float]:
    """Contraintes pariétales (paroi basse, paroi haute), positives pour un écoulement vers +x."""
    dudy = ddy(grid, U)
    return float(nu * dudy[0]), float(-nu * dudy[-1])


def bulk_velocity(grid: Grid, U: np.ndarray) -> float:
    return integrate(grid, U) / (grid.y[-1] - grid.y[0])


def centerline_velocity(grid: Grid, U: np.ndarray) -> float:
    return float(np.interp(grid.h, grid.y, U))


def cess_eddy_viscosity(grid: Grid, nu: float, u_tau: float = 1.0,
                        kappa: float = 0.426, a_plus: float = 25.4) -> np.ndarray:
    """Viscosité turbulente algébrique de Cess (forme de Reynolds & Tiederman 1967).

    Sert uniquement de champ initial : elle donne un profil turbulent réaliste.
    """
    xi = grid.wall_distance / grid.h
    re_tau = u_tau * grid.h / nu
    damping = 1.0 - np.exp(-xi * re_tau / a_plus)
    arg = (1.0 + (kappa * re_tau) ** 2 / 9.0 * (2 * xi - xi ** 2) ** 2
           * (3 - 4 * xi + 2 * xi ** 2) ** 2 * damping ** 2)
    return nu * (0.5 * np.sqrt(arg) - 0.5)


def initial_guess(grid: Grid, nu: float, u_tau: float = 1.0):
    """Profil initial (U, ν_t) cohérent avec l'équilibre (ν + ν_t) dU/dy = u_τ² (1 − d/h)."""
    nut = cess_eddy_viscosity(grid, nu, u_tau)
    d = grid.wall_distance
    g = u_tau ** 2 * np.clip(1.0 - d / grid.h, 0.0, None) / (nu + nut)
    g = np.where(grid.y <= grid.h, g, -g)
    U = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * grid.dy)))
    U[-1] = 0.0
    return U, nut
