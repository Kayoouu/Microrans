"""Choix du matériel de calcul du solveur 2D : CPU (NumPy/SciPy) ou GPU (CuPy, CUDA).

Le code des volumes finis est écrit en opérations vectorielles sur des tableaux ; il
s'exécute tel quel sur GPU en remplaçant le module de tableaux `xp` (numpy → cupy) et le
module creux (scipy.sparse → cupyx.scipy.sparse). Le maillage et la préparation (conditions
aux limites, hiérarchie AMG) restent sur CPU ; seuls les champs et opérateurs sont
transférés.

État honnête : le chemin GPU est vérifié en test par un « faux GPU » qui, comme CuPy,
refuse tout mélange implicite CPU/GPU (tests/fake_device.py), mais il n'a PAS été exécuté
sur une vraie carte graphique dans ce dépôt. Gain attendu seulement pour de grands maillages
(≳ 10⁵ cellules) : en dessous, le coût de lancement des noyaux GPU domine.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

BACKENDS = ("cpu", "gpu")


class Backend:
    def __init__(self, name, xp, sparse):
        self.name, self.xp, self.sparse = name, xp, sparse

    @property
    def is_gpu(self) -> bool:
        return self.name != "cpu"

    def asarray(self, a, dtype=None):
        return self.xp.asarray(a, dtype=dtype) if dtype is not None else self.xp.asarray(a)

    def to_host(self, a):
        if isinstance(a, dict):
            return {k: self.to_host(v) for k, v in a.items()}
        if self.xp is np or a is None or isinstance(a, (int, float)):
            return a
        get = getattr(a, "get", None)
        return get() if callable(get) else np.asarray(a)

    def __repr__(self):
        return f"Backend({self.name})"


CPU = Backend("cpu", np, sp)
_custom: dict[str, Backend] = {}


def register(name: str, backend: Backend):
    """Enregistre un backend supplémentaire (utilisé par les tests : faux GPU)."""
    _custom[name] = backend


def get_backend(name: str | Backend | None = "cpu") -> Backend:
    if isinstance(name, Backend):
        return name
    name = (name or "cpu").lower()
    if name == "cpu":
        return CPU
    if name in _custom:
        return _custom[name]
    if name in ("gpu", "cuda", "cupy"):
        try:
            import cupy
            import cupyx.scipy.sparse as csp
            if cupy.cuda.runtime.getDeviceCount() < 1:
                raise RuntimeError("aucun GPU CUDA détecté")
        except Exception as exc:                      # ImportError, CUDARuntimeError...
            raise RuntimeError(
                "Backend GPU indisponible : il faut une carte NVIDIA avec CUDA et CuPy "
                "(pip install cupy-cuda12x). Détail : " + str(exc)) from exc
        return Backend("gpu", cupy, csp)
    raise ValueError(f"Backend inconnu '{name}'. Choix : {BACKENDS}")


def available_backends() -> list[str]:
    out = ["cpu"]
    try:
        get_backend("gpu")
        out.append("gpu")
    except RuntimeError:
        pass
    return out
