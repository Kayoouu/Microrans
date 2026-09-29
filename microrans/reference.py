"""Solutions et corrélations de référence (exactes ou empiriques).

Aucune donnée DNS n'est embarquée : pour comparer à une DNS (ex. Lee & Moser 2015,
https://turbulence.oden.utexas.edu), fournir un fichier via l'option --reference du CLI.
"""
from __future__ import annotations

import numpy as np

KAPPA = 0.41
B_LOG = 5.2


def log_law(y_plus, kappa: float = KAPPA, b: float = B_LOG):
    """Loi logarithmique U+ = ln(y+)/κ + B (valable ~ 30 < y+ < 0.2 Re_τ)."""
    return np.log(y_plus) / kappa + b


def reichardt(y_plus, kappa: float = KAPPA):
    """Loi de paroi composite de Reichardt (1951), de la sous-couche visqueuse à la zone log."""
    y_plus = np.asarray(y_plus, dtype=float)
    return (np.log1p(kappa * y_plus) / kappa
            + 7.8 * (1.0 - np.exp(-y_plus / 11.0) - y_plus / 11.0 * np.exp(-y_plus / 3.0)))


def dean_bulk_velocity_plus(re_tau: float) -> float:
    """U_b/u_τ déduit de la corrélation de Dean (1978) : C_f = 0.073 Re_b^(-1/4).

    Avec C_f = 2/(U_b+)² et Re_b = 2 Re_τ U_b+ (Re_b basé sur 2h et U_b).
    Corrélation empirique (précision de quelques %), domaine ~6000 < Re_b < 6e5.
    """
    return float(((2.0 / 0.073) * (2.0 * re_tau) ** 0.25) ** (1.0 / 1.75))


def dean_cf(re_b: float) -> float:
    return 0.073 * re_b ** -0.25


def laminar_poiseuille(y, forcing: float, nu: float, h: float = 1.0):
    """Écoulement de Poiseuille plan : U = f y (2h − y) / (2ν)."""
    y = np.asarray(y, dtype=float)
    return forcing * y * (2.0 * h - y) / (2.0 * nu)


def stokes_oscillation(y, omega: float, nu: float, h: float = 1.0) -> np.ndarray:
    """Réponse complexe Û(y) à un forçage e^{iωt} d'amplitude unité, en laminaire.

    Solution de iω Û = 1 + ν Û'' avec Û(0) = Û(2h) = 0 :
        Û = (1/iω) [1 − cosh(λ(y−h)) / cosh(λh)],  λ = (1+i) √(ω/2ν).
    Forme numériquement stable (pas de débordement de cosh pour λh grand).
    """
    y = np.asarray(y, dtype=float)
    lam = (1.0 + 1.0j) * np.sqrt(omega / (2.0 * nu))
    ratio = (np.exp(lam * (y - 2.0 * h)) + np.exp(-lam * y)) / (1.0 + np.exp(-2.0 * lam * h))
    return (1.0 - ratio) / (1.0j * omega)


def womersley_channel(y, t: float, forcing_mean: float, amplitude: float, omega: float,
                      nu: float, h: float = 1.0):
    """Solution laminaire exacte (régime périodique établi) pour f(t) = f0 + A sin(ωt)."""
    # A sin(ωt) = Re[−iA e^{iωt}]
    osc = np.real(-1.0j * amplitude * stokes_oscillation(y, omega, nu, h) * np.exp(1.0j * omega * t))
    return laminar_poiseuille(y, forcing_mean, nu, h) + osc


def load_reference_profile(path: str):
    """Lit un fichier texte à 2 colonnes (y+, U+) ; lignes commençant par # ou % ignorées."""
    data = np.loadtxt(path, comments=("#", "%"))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"{path} : il faut au moins 2 colonnes (y+, U+).")
    return data[:, 0], data[:, 1]
