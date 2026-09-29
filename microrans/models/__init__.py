"""Registre des modèles de turbulence disponibles."""
from __future__ import annotations

from .base import Laminar, TurbulenceModel
from .k_epsilon import LaunderSharmaKE
from .k_omega import WilcoxKOmega2006
from .k_omega_sst import MenterSST
from .spalart_allmaras import SpalartAllmaras

MODELS: dict[str, type[TurbulenceModel]] = {
    "laminar": Laminar,
    "sa": SpalartAllmaras,
    "ke": LaunderSharmaKE,
    "kw": WilcoxKOmega2006,
    "sst": MenterSST,
}

ALIASES = {
    "spalart-allmaras": "sa",
    "k-eps": "ke", "k-epsilon": "ke", "launder-sharma": "ke",
    "k-omega": "kw", "kw2006": "kw", "wilcox": "kw",
    "k-omega-sst": "sst", "menter": "sst",
}

TURBULENT_MODELS = ("sa", "ke", "kw", "sst")


def canonical_name(name: str) -> str:
    key = name.strip().lower()
    key = ALIASES.get(key, key)
    if key not in MODELS:
        raise ValueError(f"Modèle inconnu '{name}'. Choix : {', '.join(MODELS)} "
                         f"(alias : {', '.join(ALIASES)})")
    return key


def get_model(name: str, grid, nu: float, **options) -> TurbulenceModel:
    return MODELS[canonical_name(name)](grid, nu, **options)


__all__ = ["MODELS", "TURBULENT_MODELS", "TurbulenceModel", "get_model", "canonical_name",
           "Laminar", "SpalartAllmaras", "LaunderSharmaKE", "WilcoxKOmega2006", "MenterSST"]
