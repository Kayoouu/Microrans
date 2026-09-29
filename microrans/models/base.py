"""Interface commune des modèles de turbulence (1D et 2D).

Un modèle ne décrit que la physique locale. La discrétisation est fournie par :
  - `self.ops`  : opérateurs de gradient (|∇f|², ∇f·∇g) du maillage 1D ou 2D ;
  - `flow`      : grandeurs du champ moyen (S = √(2SᵢⱼSᵢⱼ), Ω = |rot U|, dérivées secondes) ;
  - `step`      : résolution implicite d'une équation de transport (schéma en temps,
                  convection, conditions aux limites, planchers).
Le même code de modèle sert donc au solveur 1D (canal) et au solveur volumes finis 2D.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ..numerics import ddy


class Ops1D:
    """Opérateurs de gradient sur le maillage 1D (différences finies d'ordre 2)."""

    def __init__(self, grid):
        self.grid = grid

    def grad_sq(self, f, name=None):
        return ddy(self.grid, f) ** 2

    def grad_dot(self, f, g, fname=None, gname=None):
        return ddy(self.grid, f) * ddy(self.grid, g)


class TurbulenceModel(ABC):
    """Modèle de turbulence à 0, 1 ou 2 équations de transport.

    Chaque équation s'écrit  ∂φ/∂t + convection = S_exp − D_imp·φ + ∇·(Γ ∇φ)
    avec S_exp ≥ 0 (explicite) et D_imp ≥ 0 (puits implicite) : voir `linearize_source`.
    """

    name: str = "base"
    label: str = "base"
    variables: tuple[str, ...] = ()
    floors: dict[str, float] = {}

    def __init__(self, grid, nu: float):
        self.grid = grid
        self.nu = float(nu)
        d = np.array(grid.wall_distance, dtype=float, copy=True)
        if getattr(grid, "is_1d", False):
            # 1D : noeuds 0 et N sur les parois (Dirichlet), distance factice non nulle
            d[0] = d[-1] = grid.first_cell_height
            self.wall_nodes = np.array([0, len(d) - 1])
            self.ops = Ops1D(grid)
        else:
            self.wall_nodes = np.array([], dtype=int)
            self.ops = None          # fourni par le solveur 2D
        self.d = np.maximum(d, 1e-300)

    # -- à implémenter -----------------------------------------------------------
    @abstractmethod
    def initial_state(self, flow, nut0: np.ndarray) -> dict[str, np.ndarray]:
        """État initial à partir d'une estimation de la viscosité turbulente."""

    @abstractmethod
    def eddy_viscosity(self, state: dict, flow) -> np.ndarray:
        """Viscosité turbulente ν_t (nulle aux noeuds pariétaux en 1D)."""

    @abstractmethod
    def update(self, state: dict, flow, step) -> dict[str, np.ndarray]:
        """Avance les variables de transport d'une (sous-)itération implicite."""

    def freestream_values(self, velocity: float, intensity: float = 0.001,
                          viscosity_ratio: float = 0.1, length: float = 1.0) -> dict[str, float]:
        """Valeurs amont (entrée / champ lointain) à partir de I = u'/U et ν_t/ν."""
        return {}

    # -- utilitaires -------------------------------------------------------------
    def wall_value(self, name: str, d1) -> np.ndarray:
        """Valeur imposée à la paroi pour la variable `name` (d1 : distance du 1er point)."""
        return float_array(d1) * 0.0

    def wall_values(self, name: str) -> tuple[float, float]:
        """(1D) valeurs aux deux parois du canal."""
        dy = self.grid.dy
        return (float(self.wall_value(name, dy[0])), float(self.wall_value(name, dy[-1])))

    def extra_fields(self, state: dict, flow) -> dict[str, np.ndarray]:
        """Champs dérivés à écrire en sortie (ex. ε pour un modèle k-ω)."""
        return {}

    def _solve(self, step, name, gamma, source, sink):
        return step.solve(name, gamma, source, sink, model=self)

    def _zero_at_walls(self, a: np.ndarray) -> np.ndarray:
        a[self.wall_nodes] = 0.0
        return a

    def __repr__(self) -> str:
        return f"{type(self).__name__}(nu={self.nu:g})"


class Laminar(TurbulenceModel):
    """Pas de modèle : ν_t = 0 (sert à la vérification contre les solutions exactes)."""

    name = "laminar"
    label = "Laminaire (ν_t = 0)"

    def initial_state(self, flow, nut0):
        return {}

    def eddy_viscosity(self, state, flow):
        return self.d * 0.0

    def update(self, state, flow, step):
        return {}


def float_array(x):
    """Tableau de réels sans changer de module (NumPy ou CuPy) ; accepte les scalaires."""
    if hasattr(x, "astype"):
        return x.astype(float, copy=False)
    return np.asarray(x, dtype=float)


def linearize_source(phi: np.ndarray, q: np.ndarray, dq: np.ndarray):
    """Linéarisation de Newton d'un terme source local Q(φ) ≈ Q(φ⁰) + Q'(φ⁰)(φ − φ⁰).

    Renvoie (S_exp, D_imp) avec S_exp ≥ 0, D_imp ≥ 0 et S_exp − D_imp φ⁰ = Q(φ⁰) :
    la solution convergée est inchangée, mais la pente négative de Q est implicitée
    (indispensable pour les puits quadratiques βω², C2 ε²/k, c_w1 f_w (ν̃/d)²,
    qu'une linéarisation « Picard » rend oscillants aux grands pas de temps).
    """
    sink = np.maximum(-dq, 0.0)
    source = q + sink * phi
    negative = source < 0.0
    if np.any(negative):
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            extra = np.where(phi > 0.0, -source / phi, 0.0)
        sink = np.where(negative, sink + np.minimum(extra, 1e30), sink)
        source = np.where(negative, 0.0, source)
    return source, sink


def k_omega_guess(model: TurbulenceModel, nut0: np.ndarray, beta: float = 0.075,
                  beta_star: float = 0.09, kappa: float = 0.41, u_tau: float = 1.0):
    """Estimation initiale (k, ω) cohérente avec ν_t = k/ω.

    ω mélange la solution proche paroi 6ν/(β d²) et la zone logarithmique
    u_τ/(√β* κ d) ; k = ν_t ω donne k ≈ u_τ²/√β* dans la zone log et k ~ d² à la paroi.
    """
    d = model.d
    omega = np.sqrt((6.0 * model.nu / (beta * d ** 2)) ** 2
                    + (u_tau / (np.sqrt(beta_star) * kappa * d)) ** 2)
    k = nut0 * omega
    k[model.wall_nodes] = 0.0
    return k, omega


def k_omega_freestream(nu, velocity, intensity, viscosity_ratio):
    k = 1.5 * (intensity * velocity) ** 2
    k = max(k, 1e-12 * max(velocity, 1e-12) ** 2)
    omega = k / (viscosity_ratio * nu)
    return k, omega
