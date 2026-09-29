"""Spalart-Allmaras (forme « standard » du NASA Turbulence Modeling Resource).

Réf. : Spalart & Allmaras, AIAA Paper 92-0439 (1992) ;
       https://turbmodels.larc.nasa.gov/spalart.html
Par défaut la variante « SA-noft2 » (f_t2 = 0) est utilisée ; `ft2=True` active le terme
de transition f_t2. La limitation de S̃ (note c du TMR, c_v2 = 0.7, c_v3 = 0.9) est incluse.
"""
from __future__ import annotations

import numpy as np

from ..numerics import ddy
from .base import TurbulenceModel, linearize_source


class SpalartAllmaras(TurbulenceModel):
    name = "sa"
    label = "Spalart-Allmaras"
    variables = ("nu_tilde",)
    floors = {"nu_tilde": 0.0}

    cb1, cb2, sigma, kappa = 0.1355, 0.622, 2.0 / 3.0, 0.41
    cw2, cw3, cv1 = 0.3, 2.0, 7.1
    ct3, ct4 = 1.2, 0.5
    cv2, cv3 = 0.7, 0.9
    cw1 = cb1 / kappa ** 2 + (1.0 + cb2) / sigma

    def __init__(self, grid, nu, ft2: bool = False):
        super().__init__(grid, nu)
        self.use_ft2 = ft2
        if ft2:
            self.label = "Spalart-Allmaras (avec f_t2)"
        else:
            self.label = "Spalart-Allmaras (SA-noft2)"

    def fv1(self, chi):
        chi3 = chi ** 3
        return chi3 / (chi3 + self.cv1 ** 3)

    def initial_state(self, flow, nut0):
        # Inversion de ν_t = ν̃ f_v1(ν̃/ν) par dichotomie (fonction monotone croissante).
        target = nut0 / self.nu
        lo = np.zeros_like(target)
        hi = np.maximum(2.0 * target, 20.0)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            too_big = mid * self.fv1(mid) > target
            hi = np.where(too_big, mid, hi)
            lo = np.where(too_big, lo, mid)
        nt = 0.5 * (lo + hi) * self.nu
        return {"nu_tilde": self._zero_at_walls(nt)}

    def eddy_viscosity(self, state, flow):
        nt = state["nu_tilde"]
        return self._zero_at_walls(nt * self.fv1(nt / self.nu))

    def s_tilde(self, nt, omega):
        """Vorticité modifiée S̃ avec la limitation de S̄ (note c du TMR)."""
        chi = nt / self.nu
        fv1 = self.fv1(chi)
        fv2 = 1.0 - chi / (1.0 + chi * fv1)
        s_bar = nt * fv2 / (self.kappa ** 2 * self.d ** 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            s_lim = omega + omega * (self.cv2 ** 2 * omega + self.cv3 * s_bar) / (
                (self.cv3 - 2.0 * self.cv2) * omega - s_bar)
        return np.where(s_bar >= -self.cv2 * omega, omega + s_bar, s_lim)

    def local_source(self, nt, omega):
        """Terme source local net Q(ν̃) = production − destruction."""
        k2, d = self.kappa ** 2, self.d
        s_tilde = self.s_tilde(nt, omega)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(s_tilde > 1e-300, nt / (s_tilde * k2 * d ** 2), 10.0)
        r = np.minimum(r, 10.0)
        g = r + self.cw2 * (r ** 6 - r)
        fw = g * ((1.0 + self.cw3 ** 6) / (g ** 6 + self.cw3 ** 6)) ** (1.0 / 6.0)
        ft2 = self.ct3 * np.exp(-self.ct4 * (nt / self.nu) ** 2) if self.use_ft2 else 0.0
        production = self.cb1 * (1.0 - ft2) * s_tilde * nt
        destruction = (self.cw1 * fw - self.cb1 / k2 * ft2) * (nt / d) ** 2
        return production - destruction

    def update(self, state, flow, step):
        nt = state["nu_tilde"]
        omega = flow.strain
        # Jacobienne locale dQ/dν̃ par différence finie (dépendances via S̃, r, f_w, f_t2).
        q = self.local_source(nt, omega)
        dnt = 1e-7 * nt + 1e-30
        dq = (self.local_source(nt + dnt, omega) - q) / dnt
        source, sink = linearize_source(nt, q, dq)
        # Terme non conservatif c_b2/σ (∂ν̃/∂y)² ≥ 0 : explicite.
        source = source + self.cb2 / self.sigma * ddy(self.grid, nt) ** 2
        gamma = (self.nu + nt) / self.sigma
        return {"nu_tilde": self._solve(step, "nu_tilde", gamma, source, sink)}
