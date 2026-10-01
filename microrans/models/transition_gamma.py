"""k-ω SST (2003) couplé au modèle de transition à une équation γ de Menter et al. (2015).

Réf. : F.R. Menter, P.E. Smirnov, T. Liu, R. Avancha, « A One-Equation Local
       Correlation-Based Transition Model », Flow Turbulence Combust. 95 (2015) 583-619.
Forme implémentée (corrélation fondée sur Tu, sans rugosité ni écoulement transverse) :

  Dγ/Dt = P_γ − E_γ + ∇·[(ν + ν_t/σ_f) ∇γ]
  P_γ = F_length S γ (1 − γ) F_onset          E_γ = c_a2 Ω F_turb γ (c_e2 γ − 1)
  F_onset1 = Re_V / (2.2 Re_θc)   F_onset2 = min(F_onset1, 2)
  F_onset3 = max(1 − (R_T/3.5)³, 0)   F_onset = max(F_onset2 − F_onset3, 0)
  F_turb = exp(−(R_T/2)⁴),  Re_V = d² S/ν,  R_T = k/(ν ω)
  Re_θc = C_TU1 + C_TU2 exp(−C_TU3 Tu_L F_PG(λ_θL))
  Tu_L = min(100 √(2k/3) / (ω d), 100)                       (Tu local, Galiléen)
  λ_θL = −7.57e-3 (dV/dy) d²/ν + 0.0128, borné à [−1, 1] ; dV/dy = n·∇(n·U), n = ∇d
  F_PG = min(1 + C_PG1 λ, C_PG1_lim) si λ ≥ 0,
         min(1 + C_PG2 λ + C_PG3 min(λ + 0.0681, 0), C_PG2_lim) sinon ; F_PG ≥ 0
  F_length = 100, c_a2 = 0.06, c_e2 = 50, σ_f = 1 ; C_TU1 = 100, C_TU2 = 1000, C_TU3 = 1 ;
  C_PG1 = 14.68, C_PG1_lim = 1.5, C_PG2 = −7.34, C_PG2_lim = 3, C_PG3 = 0.

Couplage à l'équation de k du SST :
  P̃_k = γ P_k + P_k^lim,   D̃_k = max(γ, 0.1) β* k ω,
  P_k^lim = 5 C_k max(γ − 0.2, 0)(1 − γ) F_on^lim max(3 C_SEP ν − ν_t, 0) S Ω,
  F_on^lim = min(max(Re_V/(2.2 Re_θc^lim) − 1, 0), 3), Re_θc^lim = 1100, C_k = C_SEP = 1 ;
  F1 = max(F1_SST, F3), F3 = exp(−(R_y/120)⁸), R_y = d √k/ν (comme γ-Re_θ) ; équation de ω
  inchangée. P_k = min(ν_t S², 10 β* k ω) (SST 2003), ou ν_t S Ω (Kato-Launder,
  option `kato_launder = true`, souvent associée au modèle γ, p. ex. dans SU2).
Conditions aux limites : γ = 1 en entrée / champ lointain, gradient nul aux parois.

L'article n'ayant pas pu être consulté pendant l'écriture, les équations ont été écrites de
mémoire puis vérifiées (formes et constantes) contre l'implémentation libre de SU2
(« SLM / MENTER_SLM », PR su2code/SU2#1901). Non implémentés : rugosité, écoulement
transverse, variante Spalart-Allmaras. Solveur 1D : γ à la paroi = valeur du premier noeud
intérieur (gradient nul décalé d'une itération), dV/dy = 0.

Usage : y⁺ ≈ 1 (lois de paroi refusées), convection de la turbulence au 2e ordre
(`[solver] convection_turb = "linearUpwindLimited"` : l'upwind avance de 36 % la transition
du cas T3A-), turbulence amont réglée sur la décroissance mesurée (`freestream_decay`).
Validation et limites : docs/notes_transition.md.
"""
from __future__ import annotations

from ._xp import anp

from ..linalg import array_module
from .base import float_array, linearize_source
from .k_omega_sst import MenterSST


class MenterSSTGamma(MenterSST):
    name = "sst_gamma"
    label = "k-ω SST + transition γ (Menter 2015)"
    variables = ("k", "omega", "gamma")
    floors = {"k": 1e-14, "omega": 1e-10, "gamma": 0.0}
    wall_zero_gradient = ("gamma",)

    f_length, c_a2, c_e2, sigma_f = 100.0, 0.06, 50.0, 1.0
    c_tu1, c_tu2, c_tu3 = 100.0, 1000.0, 1.0
    c_pg1, c_pg1_lim, c_pg2, c_pg2_lim, c_pg3 = 14.68, 1.5, -7.34, 3.0, 0.0
    re_thc_lim, c_k, c_sep = 1100.0, 1.0, 1.0

    def __init__(self, grid, nu, kato_launder: bool = False):
        super().__init__(grid, nu)
        self.kato_launder = bool(kato_launder)
        if self.kato_launder:
            self.label = self.label + ", production de Kato-Launder"
        # normale à la paroi (2D) : sert à dV/dy du paramètre de gradient de pression λ_θL
        self._normal_host = None if getattr(grid, "is_1d", False) else grid.wall_normal()
        self._normal = (None, None)
        self._gamma_wall = (1.0, 1.0)

    # -- conditions limites, état initial ---------------------------------------------
    def wall_value(self, name, d1):
        if name == "gamma":
            return float_array(d1) * 0.0 + 1.0      # non utilisé en 2D (gradient nul)
        return super().wall_value(name, d1)

    def wall_values(self, name):
        if name == "gamma":
            return self._gamma_wall                 # 1D : gradient nul (décalé)
        return super().wall_values(name)

    def initial_state(self, flow, nut0):
        st = super().initial_state(flow, nut0)
        st["gamma"] = st["k"] * 0.0 + 1.0
        return st

    def freestream_values(self, velocity, intensity=0.001, viscosity_ratio=0.1, length=1.0):
        out = super().freestream_values(velocity, intensity, viscosity_ratio, length)
        out["gamma"] = 1.0
        return out

    # -- fonctions du modèle --------------------------------------------------------------
    def blending(self, k, w, cross):
        f1, f2 = super().blending(k, w, cross)
        ry = self.d * anp.sqrt(anp.maximum(k, 0.0)) / self.nu
        f3 = anp.exp(-(ry / 120.0) ** 8)
        return anp.maximum(f1, f3), f2

    def _wall_normal(self):
        """Normale pariétale sur le même matériel que les champs (transfert au 1er appel)."""
        xp = array_module(self.d)
        if self._normal[0] is not xp:
            self._normal = (xp, xp.asarray(self._normal_host))
        return self._normal[1]

    def _dvdy(self, flow):
        """Dérivée normale de la vitesse normale à la paroi, n·∇(n·U) = nᵢ nⱼ ∂uᵢ/∂xⱼ
        (n est constant le long de la normale : n·∇n = 0)."""
        g = getattr(flow, "gradU", None)
        if g is None or self._normal_host is None:
            return self.d * 0.0                     # 1D (canal) : V = 0
        n = self._wall_normal()
        if n.shape[1] == 3:
            xp = array_module(self.d)
            return xp.einsum("ci,cij,cj->c", n, g, n)
        nx, ny = n[:, 0], n[:, 1]
        return nx * nx * g[:, 0, 0] + nx * ny * (g[:, 0, 1] + g[:, 1, 0]) + ny * ny * g[:, 1, 1]

    def re_theta_c(self, k, w, flow):
        """Reynolds critique Re_θc(Tu_L, λ_θL) ; renvoie aussi Tu_L et λ_θL."""
        d, nu = self.d, self.nu
        tu_l = anp.minimum(100.0 * anp.sqrt(2.0 * anp.maximum(k, 0.0) / 3.0) / (w * d), 100.0)
        lam = anp.clip(-7.57e-3 * self._dvdy(flow) * d ** 2 / nu + 0.0128, -1.0, 1.0)
        f_pos = anp.minimum(1.0 + self.c_pg1 * lam, self.c_pg1_lim)
        f_neg = anp.minimum(1.0 + self.c_pg2 * lam + self.c_pg3 * anp.minimum(lam + 0.0681, 0.0),
                           self.c_pg2_lim)
        f_pg = anp.maximum(anp.where(lam >= 0.0, f_pos, f_neg), 0.0)
        re_thc = self.c_tu1 + self.c_tu2 * anp.exp(-self.c_tu3 * tu_l * f_pg)
        return re_thc, tu_l, lam

    def onset_functions(self, k, w, flow, strain=None):
        """(F_onset, F_turb, Re_V, Re_θc, Tu_L)."""
        nu = self.nu
        re_v = self.d ** 2 * (flow.strain if strain is None else strain) / nu
        r_t = k / (nu * w)
        re_thc, tu_l, _ = self.re_theta_c(k, w, flow)
        f_onset2 = anp.minimum(re_v / (2.2 * re_thc), 2.0)
        f_onset3 = anp.maximum(1.0 - (r_t / 3.5) ** 3, 0.0)
        f_onset = anp.maximum(f_onset2 - f_onset3, 0.0)
        f_turb = anp.exp(-(r_t / 2.0) ** 4)
        return f_onset, f_turb, re_v, re_thc, tu_l

    # -- une itération ------------------------------------------------------------------
    def update(self, state, flow, step):
        nu, bs = self.nu, self.beta_star
        k, w, g = state["k"], state["omega"], state["gamma"]
        f1, f2 = self.blending(k, w, self.ops.grad_dot(k, w, "k", "omega"))
        nut = self._nut(k, w, f2, flow)
        S, Om = flow.strain, flow.vorticity
        s2 = S ** 2

        # 1) intermittence γ : sources linéarisées par Newton
        f_onset, f_turb, re_v, _, _ = self.onset_functions(k, w, flow, S)
        prod = self.f_length * S * f_onset
        dest = self.c_a2 * Om * f_turb
        q = prod * g * (1.0 - g) - dest * g * (self.c_e2 * g - 1.0)
        dq = prod * (1.0 - 2.0 * g) - dest * (2.0 * self.c_e2 * g - 1.0)
        source, sink = linearize_source(g, q, dq)
        if len(self.wall_nodes):
            self._gamma_wall = (float(g[1]), float(g[-2]))
        g_new = self._solve(step, "gamma", nu + nut / self.sigma_f, source, sink)
        g_new = anp.minimum(g_new, 1.0)

        # 2) k : production × γ (+ P_k^lim), destruction × max(γ, 0.1)
        def blend(c1, c2):
            return f1 * c1 + (1.0 - f1) * c2

        sigma_k = blend(self.sigma_k1, self.sigma_k2)
        sigma_w = blend(self.sigma_w1, self.sigma_w2)
        beta = blend(self.beta1, self.beta2)
        gamma_w = blend(self.gamma1, self.gamma2)

        p_raw = nut * S * Om if self.kato_launder else nut * s2
        prod_k = anp.minimum(p_raw, 10.0 * bs * k * w)
        f_on_lim = anp.minimum(anp.maximum(re_v / (2.2 * self.re_thc_lim) - 1.0, 0.0), 3.0)
        p_lim = (5.0 * self.c_k * anp.maximum(g_new - 0.2, 0.0) * (1.0 - g_new) * f_on_lim
                 * anp.maximum(3.0 * self.c_sep * nu - nut, 0.0) * S * Om)
        k_new = self._solve(step, "k", nu + sigma_k * nut, g_new * prod_k + p_lim,
                            anp.maximum(g_new, 0.1) * bs * w)

        # 3) ω : équation du SST inchangée
        cd = 2.0 * (1.0 - f1) * self.sigma_w2 * self.ops.grad_dot(k_new, w, "k", "omega") / w
        q = gamma_w * s2 + cd - beta * w ** 2
        dq = -2.0 * beta * w - cd / w
        source, sink = linearize_source(w, q, dq)
        w_new = self._solve(step, "omega", nu + sigma_w * nut, source, sink)
        return {"k": k_new, "omega": w_new, "gamma": g_new}

    def extra_fields(self, state, flow):
        out = super().extra_fields(state, flow)
        k, w = state["k"], state["omega"]
        f_onset, _, _, re_thc, tu_l = self.onset_functions(k, w, flow)
        out.update({"F_onset": f_onset, "Re_theta_c": re_thc, "Tu_L": tu_l})
        return out


def freestream_decay(distance, velocity, nu, intensity, viscosity_ratio, beta=0.0828,
                     beta_star=0.09):
    """Décroissance de la turbulence amont prédite par le SST (F1 = 0 hors couche limite,
    donc β = β₂) : écoulement uniforme U, sans production.

      ω(x) = ω₀ / (1 + β ω₀ x/U),   k(x) = k₀ (1 + β ω₀ x/U)^(−β*/β)
      Tu(x) = Tu₀ (1 + β ω₀ x/U)^(−β*/(2β)),   k₀ = 1.5 (Tu₀ U)²,  ω₀ = k₀ / (ν · ν_t/ν)

    `distance` : distance parcourue depuis l'entrée. Renvoie (Tu, ν_t/ν) à cette distance.
    Sert à choisir Tu et ν_t/ν d'entrée pour retrouver le Tu mesuré au bord d'attaque et
    sa décroissance le long de la plaque (voir docs/notes_transition.md)."""
    x = anp.asarray(distance, dtype=float)
    k0 = 1.5 * (intensity * velocity) ** 2
    w0 = k0 / (nu * viscosity_ratio)
    f = 1.0 + beta * w0 * x / velocity
    tu = intensity * f ** (-beta_star / (2.0 * beta))
    return tu, viscosity_ratio * f ** (1.0 - beta_star / beta)       # ν_t/ν = k/(ω ν)
