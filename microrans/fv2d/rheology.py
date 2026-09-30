"""Fluides non newtoniens (viscosité fonction du taux de cisaillement, « Generalized
Newtonian fluids ») : ν = ν(γ̇), γ̇ = √(2 S:S) (terme (u_r/r)² inclus en axisymétrique).

[physics.viscosity] model = … (grandeurs cinématiques : viscosités en m²/s, contraintes
divisées par ρ) :
  power_law         ν = K γ̇^(n−1)                                  (Ostwald-de Waele)
  carreau           ν = ν∞ + (ν0 − ν∞) [1 + (λγ̇)^a]^((n−1)/a)       (a = 2 : Carreau ;
                                                                    autre a : Carreau-Yasuda)
  cross             ν = ν∞ + (ν0 − ν∞) / [1 + (m γ̇)^n]               (CrossPowerLaw d'OpenFOAM)
  herschel_bulkley  ν = (τ_y + K γ̇^n) / γ̇                           (seuil d'écoulement τ_y)
  bingham           herschel_bulkley avec n = 1 (K = viscosité plastique)
  casson            ν = (√(τ_y/γ̇) + √ν∞)²                           (sang, chocolat…)
Bornes nu_min / nu_max (toutes lois) : pour les lois non bornées (loi puissance, seuil),
défaut 10⁻³ et 10³ × ν(γ̇_ref), γ̇_ref = U_ref / L_ref. Pour les fluides à seuil, nu_max
joue le rôle de la viscosité de la zone « solide » (régularisation bi-visqueuse, comme le
HerschelBulkley d'OpenFOAM) : la zone non cisaillée est un bouchon très visqueux, pas un
solide rigide ; l'erreur sur la vitesse diminue quand nu_max augmente (voir README).

Stationnaire : ν est sous-relaxé (relax, défaut 0.5). Écoulement entraîné par une force
volumique (sans entrée) : relax_U = 1 conseillé ([solver]), sinon la mise en vitesse est
très lente (même limite qu'en newtonien, voir README).

Les modèles de turbulence RANS supposent un fluide newtonien : la combinaison avec une
viscosité non newtonienne est refusée (non validée).
"""
from __future__ import annotations

import numpy as np

MODELS = ("newtonian", "power_law", "carreau", "cross", "herschel_bulkley", "bingham",
          "casson")
PARAMS = {"power_law": ("K", "n"), "carreau": ("nu0", "nu_inf", "lambda", "n"),
          "cross": ("nu0", "nu_inf", "m", "n"), "herschel_bulkley": ("tau_y", "K", "n"),
          "bingham": ("tau_y", "K"), "casson": ("tau_y", "nu_inf")}


class Rheology:
    def __init__(self, spec: dict, gamma_ref: float = 1.0):
        spec = dict(spec)
        self.model = str(spec.pop("model", "newtonian")).lower()
        if self.model not in MODELS:
            raise ValueError(f"[physics.viscosity] model = '{self.model}' inconnu. "
                             f"Choix : {MODELS}")
        missing = [k for k in PARAMS.get(self.model, ()) if k not in spec]
        if missing:
            raise ValueError(f"[physics.viscosity] {self.model} : paramètres manquants "
                             f"{missing} (attendus : {PARAMS[self.model]}).")
        self.p = {k: float(v) for k, v in spec.items() if k not in ("nu_min", "nu_max",
                                                                    "relax")}
        if self.model == "bingham":
            self.p["n"] = 1.0
        # sous-relaxation de ν en stationnaire (itération de Picard) : sans elle, la loi
        # puissance rhéofluidifiante oscille et diverge au démarrage (ν = nu_max où U = 0)
        self.relax = float(spec.get("relax", 0.5))
        self.gamma_ref = float(gamma_ref)
        nref = float(self._raw(np.array([self.gamma_ref]), np)[0])
        self.nu_ref = nref
        self.nu_min = float(spec.get("nu_min", 1e-3 * nref))
        self.nu_max = float(spec.get("nu_max", 1e3 * nref))
        if not 0.0 <= self.nu_min < self.nu_max:
            raise ValueError("[physics.viscosity] : il faut 0 ≤ nu_min < nu_max.")

    @property
    def newtonian(self) -> bool:
        return self.model == "newtonian"

    def _raw(self, g, xp):
        p = self.p
        g = xp.maximum(g, 1e-300)
        m = self.model
        if m == "power_law":
            return p["K"] * g ** (p["n"] - 1.0)
        if m == "carreau":
            a = p.get("a", 2.0)
            return p["nu_inf"] + (p["nu0"] - p["nu_inf"]) * (1.0 + (p["lambda"] * g) ** a) ** (
                (p["n"] - 1.0) / a)
        if m == "cross":
            return p["nu_inf"] + (p["nu0"] - p["nu_inf"]) / (1.0 + (p["m"] * g) ** p["n"])
        if m in ("herschel_bulkley", "bingham"):
            return (p["tau_y"] + p["K"] * g ** p["n"]) / g
        if m == "casson":
            return (xp.sqrt(p["tau_y"] / g) + np.sqrt(p["nu_inf"])) ** 2
        return xp.full_like(g, p.get("nu", 0.0))

    def __call__(self, gamma_dot, xp=np):
        return xp.clip(self._raw(gamma_dot, xp), self.nu_min, self.nu_max)

    def describe(self) -> str:
        args = ", ".join(f"{k}={v:g}" for k, v in self.p.items())
        return f"{self.model} ({args}; ν ∈ [{self.nu_min:.3g}, {self.nu_max:.3g}])"


def reference_viscosity(spec: dict, gamma_ref: float) -> float:
    """Viscosité au taux de cisaillement de référence (sert de ν pour Re, Pr, l'échelle
    des résidus) quand [physics] nu n'est pas donné."""
    return Rheology(spec, gamma_ref).nu_ref
