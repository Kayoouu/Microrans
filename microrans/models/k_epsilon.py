"""k-ε bas-Reynolds de Launder & Sharma (1974), intégré jusqu'à la paroi.

Réf. : Launder & Sharma, Letters in Heat and Mass Transfer 1, 131-138 (1974) ;
       https://turbmodels.larc.nasa.gov/ke-ls.html

Un k-ε « standard » (haut-Reynolds) ne peut pas être intégré jusqu'à la paroi : il faut
soit des fonctions d'amortissement (ce qui est fait ici), soit des lois de paroi.
La variable transportée est la dissipation « isotrope » ε̃ = ε − D, nulle à la paroi.
"""
from __future__ import annotations

from ._xp import anp

from .base import TurbulenceModel, k_omega_freestream, k_omega_guess, linearize_source


class LaunderSharmaKE(TurbulenceModel):
    name = "ke"
    label = "k-ε Launder-Sharma"
    variables = ("k", "eps")
    floors = {"k": 1e-14, "eps": 1e-14}

    c_mu, c1, c2, sigma_k, sigma_e = 0.09, 1.44, 1.92, 1.0, 1.3

    def _rt(self, k, eps):
        return k ** 2 / (self.nu * anp.maximum(eps, self.floors["eps"]))

    def f_mu(self, rt):
        return anp.exp(-3.4 / (1.0 + rt / 50.0) ** 2)

    def initial_state(self, flow, nut0):
        k, omega = k_omega_guess(self, nut0)
        eps = self._zero_at_walls(0.09 * k * omega)
        return {"k": k, "eps": eps}

    def freestream_values(self, velocity, intensity=0.001, viscosity_ratio=0.1, length=1.0):
        k, omega = k_omega_freestream(self.nu, velocity, intensity, viscosity_ratio)
        return {"k": k, "eps": self.c_mu * k * omega}

    def eddy_viscosity(self, state, flow):
        k, eps = state["k"], state["eps"]
        nut = self.c_mu * self.f_mu(self._rt(k, eps)) * k ** 2 / anp.maximum(
            eps, self.floors["eps"])
        return self._zero_at_walls(nut)

    def _d_term(self, k):
        """D = 2ν |∇√k|² (dissipation pariétale)."""
        return 2.0 * self.nu * self.ops.grad_sq(anp.sqrt(anp.maximum(k, 0.0)), "sqrt_k")

    def update(self, state, flow, step):
        nu = self.nu
        k, eps = state["k"], state["eps"]
        s2 = flow.strain ** 2
        nut = self.eddy_viscosity(state, flow)

        # Équation de k : P_k − ε̃ − D
        d_term = self._d_term(k)
        k_safe = anp.maximum(k, self.floors["k"])
        k_new = self._solve(step, "k", nu + nut / self.sigma_k, nut * s2,
                            (eps + d_term) / k_safe)

        # Équation de ε̃ : C1 (ε̃/k) P_k − C2 f2 ε̃²/k + E, E = 2 ν ν_t (∂²Uᵢ/∂xⱼ∂xₖ)².
        # C1 (ε̃/k) P_k est écrit C1 C_μ f_μ k S² (identique, sans division par k).
        rt = self._rt(k_new, eps)
        f2 = 1.0 - 0.3 * anp.exp(-rt ** 2)
        e_term = 2.0 * nu * nut * flow.second_derivative_sq
        k_safe = anp.maximum(k_new, self.floors["k"])
        # Newton sur le puits quadratique −C2 f2 ε̃²/k.
        q = self.c1 * self.c_mu * self.f_mu(rt) * k_new * s2 + e_term \
            - self.c2 * f2 * eps ** 2 / k_safe
        dq = -2.0 * self.c2 * f2 * eps / k_safe
        source, sink = linearize_source(eps, q, dq)
        eps_new = self._solve(step, "eps", nu + nut / self.sigma_e, source, sink)
        return {"k": k_new, "eps": eps_new}

    def extra_fields(self, state, flow):
        return {"eps_total": state["eps"] + self._d_term(state["k"])}
