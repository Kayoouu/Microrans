"""k-ω de Wilcox (2006).

Réf. : Wilcox, AIAA J. 46(11), 2823-2838 (2008) ; « Turbulence Modeling for CFD », 3e éd. (2006) ;
       https://turbmodels.larc.nasa.gov/wilcox.html
Inclut le limiteur de contrainte (C_lim = 7/8) et la diffusion croisée (σ_do = 1/8).
Condition pariétale : ω_w = 10 · 6ν / (β Δy₁²) (Menter 1994), Δy₁ = distance du 1er point.
"""
from __future__ import annotations

import numpy as np

from .base import (TurbulenceModel, float_array, k_omega_freestream, k_omega_guess,
                   linearize_source)


class WilcoxKOmega2006(TurbulenceModel):
    name = "kw"
    label = "k-ω Wilcox 2006"
    variables = ("k", "omega")
    floors = {"k": 1e-14, "omega": 1e-10}

    alpha, beta, beta_star = 13.0 / 25.0, 0.0708, 0.09
    sigma, sigma_star, sigma_do, c_lim = 0.5, 0.6, 1.0 / 8.0, 7.0 / 8.0

    def wall_value(self, name, d1):
        d1 = float_array(d1)
        if name == "omega":
            return 60.0 * self.nu / (self.beta * d1 ** 2)
        return d1 * 0.0

    def initial_state(self, flow, nut0):
        k, omega = k_omega_guess(self, nut0, beta=self.beta)
        if len(self.wall_nodes):
            omega[0], omega[-1] = self.wall_values("omega")
        return {"k": k, "omega": omega}

    def freestream_values(self, velocity, intensity=0.001, viscosity_ratio=0.1, length=1.0):
        k, omega = k_omega_freestream(self.nu, velocity, intensity, viscosity_ratio)
        return {"k": k, "omega": omega}

    def _omega_tilde(self, omega, flow):
        return np.maximum(omega, self.c_lim * flow.strain / np.sqrt(self.beta_star))

    def eddy_viscosity(self, state, flow):
        return self._zero_at_walls(state["k"] / self._omega_tilde(state["omega"], flow))

    def update(self, state, flow, step):
        nu = self.nu
        k, w = state["k"], state["omega"]
        s2 = flow.strain ** 2
        w_t = self._omega_tilde(w, flow)
        nut = self._zero_at_walls(k / w_t)

        k_new = self._solve(step, "k", nu + self.sigma_star * k / w, nut * s2,
                            self.beta_star * w)

        cross = self.ops.grad_dot(k_new, w, "k", "omega")
        cross_diff = np.where(cross > 0.0, self.sigma_do * cross / w, 0.0)
        # α (ω/k) P_k = α (ω/ω̃) S² ; Newton sur −βω² et sur la diffusion croisée (∝ 1/ω).
        q = self.alpha * (w / w_t) * s2 + cross_diff - self.beta * w ** 2
        dq = -2.0 * self.beta * w - cross_diff / w
        source, sink = linearize_source(w, q, dq)
        w_new = self._solve(step, "omega", nu + self.sigma * k_new / w, source, sink)
        return {"k": k_new, "omega": w_new}

    def extra_fields(self, state, flow):
        return {"eps_total": self.beta_star * state["k"] * state["omega"]}
