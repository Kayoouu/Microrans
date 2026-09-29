"""Interface commune des modèles de turbulence."""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ..flow import FlowField
from ..grid import Grid


class TurbulenceModel(ABC):
    """Modèle de turbulence à 0, 1 ou 2 équations de transport.

    Chaque équation s'écrit sous la forme
        ∂φ/∂t = S_exp − D_imp·φ + ∂/∂y(Γ ∂φ/∂y)
    avec S_exp ≥ 0 (termes explicites) et D_imp ≥ 0 (puits implicités, linéarisation
    de type Patankar) pour préserver la positivité. Le schéma en temps est fourni par
    l'objet `step` (voir solver.ImplicitStep) : le modèle ne décrit que la physique.
    """

    name: str = "base"
    label: str = "base"
    variables: tuple[str, ...] = ()
    floors: dict[str, float] = {}

    def __init__(self, grid: Grid, nu: float):
        self.grid = grid
        self.nu = float(nu)
        d = grid.wall_distance.copy()
        # Aux noeuds pariétaux (d = 0) on met une valeur factice : ces noeuds sont en Dirichlet.
        d[0] = d[-1] = grid.first_cell_height
        self.d = d

    # -- à implémenter -----------------------------------------------------------
    @abstractmethod
    def initial_state(self, flow: FlowField, nut0: np.ndarray) -> dict[str, np.ndarray]:
        """État initial à partir d'une estimation de la viscosité turbulente."""

    @abstractmethod
    def eddy_viscosity(self, state: dict, flow: FlowField) -> np.ndarray:
        """Viscosité turbulente ν_t aux noeuds (nulle aux parois)."""

    @abstractmethod
    def update(self, state: dict, flow: FlowField, step) -> dict[str, np.ndarray]:
        """Avance les variables de transport d'une (sous-)itération implicite."""

    # -- utilitaires -------------------------------------------------------------
    def wall_values(self, name: str) -> tuple[float, float]:
        return (0.0, 0.0)

    def extra_fields(self, state: dict, flow: FlowField) -> dict[str, np.ndarray]:
        """Champs dérivés à écrire en sortie (ex. ε pour un modèle k-ω)."""
        return {}

    def _solve(self, step, name, gamma, source, sink):
        phi = step.solve(name, gamma, source, sink, self.wall_values(name))
        floor = self.floors.get(name)
        if floor is not None:
            phi[1:-1] = np.maximum(phi[1:-1], floor)
        return phi

    @staticmethod
    def _zero_at_walls(a: np.ndarray) -> np.ndarray:
        a[0] = a[-1] = 0.0
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
        return np.zeros_like(flow.U)

    def update(self, state, flow, step):
        return {}


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
    k[0] = k[-1] = 0.0
    return k, omega
