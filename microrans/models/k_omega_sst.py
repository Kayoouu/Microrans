"""k-ω SST de Menter, version 2003.

Réf. : Menter, Kuntz & Langtry, « Ten Years of Industrial Experience with the SST
       Turbulence Model », Turbulence, Heat and Mass Transfer 4 (2003) ;
       https://turbmodels.larc.nasa.gov/sst.html
ν_t = a1 k / max(a1 ω, S F2), limiteur de production 10 β* k ω, CD_kω borné à 1e-10.
Condition pariétale : ω_w = 10 · 6ν / (β1 Δy₁²).
"""
from __future__ import annotations

import numpy as np

from ..numerics import ddy
from .base import TurbulenceModel, k_omega_guess, linearize_source


class MenterSST(TurbulenceModel):
    name = "sst"
    label = "k-ω SST (Menter 2003)"
    variables = ("k", "omega")
    floors = {"k": 1e-14, "omega": 1e-10}

    a1, beta_star = 0.31, 0.09
    sigma_k1, sigma_w1, beta1, gamma1 = 0.85, 0.5, 0.075, 5.0 / 9.0
    sigma_k2, sigma_w2, beta2, gamma2 = 1.0, 0.856, 0.0828, 0.44

    def omega_wall(self, dy1: float) -> float:
        return 60.0 * self.nu / (self.beta1 * dy1 ** 2)

    def wall_values(self, name):
        if name == "omega":
            return (self.omega_wall(self.grid.dy[0]), self.omega_wall(self.grid.dy[-1]))
        return (0.0, 0.0)

    def initial_state(self, flow, nut0):
        k, omega = k_omega_guess(self, nut0, beta=self.beta1)
        omega[0], omega[-1] = self.wall_values("omega")
        return {"k": k, "omega": omega}

    def blending(self, k, w, dk, dw):
        """Fonctions de raccordement F1 et F2."""
        nu, d, bs = self.nu, self.d, self.beta_star
        sqk = np.sqrt(np.maximum(k, 0.0))
        cd_kw = np.maximum(2.0 * self.sigma_w2 * dk * dw / w, 1e-10)
        arg1 = np.minimum(np.maximum(sqk / (bs * w * d), 500.0 * nu / (d ** 2 * w)),
                          4.0 * self.sigma_w2 * k / (cd_kw * d ** 2))
        arg2 = np.maximum(2.0 * sqk / (bs * w * d), 500.0 * nu / (d ** 2 * w))
        return np.tanh(arg1 ** 4), np.tanh(arg2 ** 2)

    def _nut(self, k, w, f2, flow):
        return self._zero_at_walls(self.a1 * k / np.maximum(self.a1 * w, flow.strain * f2))

    def eddy_viscosity(self, state, flow):
        k, w = state["k"], state["omega"]
        _, f2 = self.blending(k, w, ddy(self.grid, k), ddy(self.grid, w))
        return self._nut(k, w, f2, flow)

    def update(self, state, flow, step):
        nu, bs = self.nu, self.beta_star
        k, w = state["k"], state["omega"]
        dk, dw = ddy(self.grid, k), ddy(self.grid, w)
        f1, f2 = self.blending(k, w, dk, dw)
        nut = self._nut(k, w, f2, flow)
        s2 = flow.dudy ** 2

        def blend(c1, c2):
            return f1 * c1 + (1.0 - f1) * c2

        sigma_k = blend(self.sigma_k1, self.sigma_k2)
        sigma_w = blend(self.sigma_w1, self.sigma_w2)
        beta = blend(self.beta1, self.beta2)
        gamma = blend(self.gamma1, self.gamma2)

        prod_k = np.minimum(nut * s2, 10.0 * bs * k * w)
        k_new = self._solve(step, "k", nu + sigma_k * nut, prod_k, bs * w)

        # Diffusion croisée (∝ 1/ω) et −βω² linéarisés par Newton.
        cd = 2.0 * (1.0 - f1) * self.sigma_w2 * ddy(self.grid, k_new) * dw / w
        q = gamma * s2 + cd - beta * w ** 2
        dq = -2.0 * beta * w - cd / w
        source, sink = linearize_source(w, q, dq)
        w_new = self._solve(step, "omega", nu + sigma_w * nut, source, sink)
        return {"k": k_new, "omega": w_new}

    def extra_fields(self, state, flow):
        k, w = state["k"], state["omega"]
        f1, _ = self.blending(k, w, ddy(self.grid, k), ddy(self.grid, w))
        return {"eps_total": self.beta_star * k * w, "F1": f1}
