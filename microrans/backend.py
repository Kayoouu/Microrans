"""Choix du matériel de calcul du solveur 2D.

Le code des volumes finis est écrit en opérations vectorielles sur des tableaux ; il
s'exécute tel quel sur un autre matériel en remplaçant le module de tableaux `xp` (numpy →
cupy, dpnp) et le module creux. Le maillage et la préparation (conditions aux limites,
hiérarchie AMG) restent sur CPU ; seuls les champs et opérateurs sont transférés.

Backends :
- cpu   : NumPy / SciPy (défaut) ;
- cuda  : cartes NVIDIA, CuPy (alias historique : gpu) ;
- rocm  : cartes AMD, CuPy compilé pour ROCm (Linux) — même code que cuda ;
- intel : cartes et puces graphiques Intel (Arc, Iris Xe, UHD) et processeurs, via dpnp /
          oneAPI (SYCL sur Level Zero ou OpenCL). Variante intel:opencl:gpu, intel:cpu…
          (sélecteur ONEAPI_DEVICE_SELECTOR). Le solveur calcule en double précision :
          refusé si le matériel n'a pas le FP64 (message explicite).

État honnête : cuda/rocm vérifiés par un « faux GPU » (tests/fake_device.py) mais jamais
exécutés sur une vraie carte dans ce dépôt ; intel vérifié sur processeur via le runtime
OpenCL CPU d'Intel (mêmes résultats que NumPy), jamais sur une vraie carte Intel. Gain
attendu seulement pour de grands maillages (≳ 10⁵ cellules) : en dessous, le coût de
lancement des noyaux domine. Une puce intégrée partage la mémoire (et son débit) avec le
processeur : gain limité.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

BACKENDS = ("cpu", "cuda", "rocm", "intel")


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
        return to_host_array(a)

    def __repr__(self):
        return f"Backend({self.name})"


def to_host_array(a):
    """Copie CPU (NumPy) d'un tableau de n'importe quel backend."""
    if isinstance(a, np.ndarray) or a is None or isinstance(a, (int, float)):
        return a
    if type(a).__module__.startswith("dpnp"):
        import dpnp
        return dpnp.asnumpy(a)
    get = getattr(a, "get", None)
    return get() if callable(get) else np.asarray(a)


CPU = Backend("cpu", np, sp)
_custom: dict[str, Backend] = {}


def register(name: str, backend: Backend):
    """Enregistre un backend supplémentaire (utilisé par les tests : faux GPU)."""
    _custom[name] = backend


def frozen() -> bool:
    """Exécutable autonome (PyInstaller) : pas de bibliothèque de carte graphique."""
    import sys
    return bool(getattr(sys, "frozen", False))


_FROZEN_HINT = ("l'exécutable calcule uniquement sur le processeur (CPU) : les bibliothèques "
                "de cartes graphiques (CuPy, dpnp) n'y sont pas incluses. Choisir « CPU », ou "
                "installer la version Python (README § 2 et § 5.3) pour essayer une carte.")


def available(name: str) -> bool:
    """La bibliothèque du backend est-elle installée ? (test rapide, sans l'importer ;
    ne garantit pas qu'une carte compatible soit présente)."""
    import importlib.util
    name = name.lower()
    if name == "cpu" or name in _custom:
        return True
    if frozen():
        return False
    mod = "dpnp" if name.startswith(("intel", "dpnp")) else "cupy"
    return importlib.util.find_spec(mod) is not None


def get_backend(name: str | Backend | None = "cpu") -> Backend:
    if isinstance(name, Backend):
        return name
    name = (name or "cpu").lower()
    if name == "cpu":
        return CPU
    if name in _custom:
        return _custom[name]
    if name not in BACKENDS + ("gpu", "cupy") and not name.startswith(("intel", "dpnp")):
        raise ValueError(f"Backend inconnu '{name}'. Choix : {', '.join(BACKENDS)}")
    if frozen():
        raise RuntimeError(f"Backend {name} indisponible : {_FROZEN_HINT}")
    if name in ("gpu", "cuda", "cupy", "rocm"):
        amd = name == "rocm"
        try:
            import cupy
            import cupyx.scipy.sparse as csp
            if cupy.cuda.runtime.getDeviceCount() < 1:
                raise RuntimeError("aucune carte détectée")
        except Exception as exc:                      # ImportError, CUDARuntimeError...
            hint = ("une carte AMD, ROCm (Linux) et CuPy pour ROCm (pip install "
                    "cupy-rocm-6-0)" if amd else
                    "une carte NVIDIA avec CUDA et CuPy (pip install cupy-cuda12x)")
            raise RuntimeError(f"Backend {name} indisponible : il faut {hint}. Détail : "
                               + str(exc)) from exc
        return Backend("rocm" if amd else "gpu", cupy, csp)
    if name.startswith("intel") or name.startswith("dpnp"):
        return _intel_backend(name)
    raise ValueError(f"Backend inconnu '{name}'. Choix : {BACKENDS}")


def _intel_backend(name: str) -> Backend:
    """dpnp (oneAPI). « intel » : carte graphique Intel (Level Zero, sinon OpenCL) ;
    « intel:<sélecteur> » : sélecteur oneAPI explicite (opencl:gpu, level_zero:gpu, cpu…)."""
    import os
    import sys
    sel = name.split(":", 1)[1] if ":" in name else None
    if sel:
        # lu par le runtime SYCL à son initialisation : sans effet si dpctl est déjà chargé
        if "dpctl" in sys.modules:
            import warnings
            warnings.warn("Sélecteur de matériel ignoré : oneAPI déjà initialisé dans ce "
                          "processus (définir ONEAPI_DEVICE_SELECTOR avant le lancement).")
        os.environ["ONEAPI_DEVICE_SELECTOR"] = sel if ":" in sel else f"*:{sel}"
    try:
        import dpctl
        import dpnp
    except ImportError as exc:
        raise RuntimeError("Backend intel indisponible : installer dpnp (pip install dpnp, "
                           "~2.5 Go avec oneMKL) et le pilote graphique Intel à jour.") from exc
    try:
        dev = dpctl.select_default_device() if sel else dpctl.select_gpu_device()
    except Exception as exc:                          # noqa: BLE001 — aucun matériel SYCL
        raise RuntimeError(f"Backend {name} : aucun matériel oneAPI trouvé ({exc}). "
                           f"Appareils visibles : {list_devices()}") from exc
    if not dev.has_aspect_fp64:
        raise RuntimeError(
            f"Backend {name} : « {dev.name} » ne calcule pas en double précision (FP64), "
            f"indispensable au solveur. Utiliser backend = \"cpu\".")
    from .sparse_generic import GenericSparseModule
    be = Backend(f"intel ({dev.name})", dpnp, GenericSparseModule(dpnp))
    be.device = dev
    return be


def list_devices() -> list[dict]:
    """Matériels de calcul visibles (CUDA/ROCm via CuPy, oneAPI via dpctl, OpenCL via
    PyOpenCL), avec la double précision : sert au diagnostic « microrans devices »."""
    out = [{"backend": "cpu", "name": "processeur (NumPy)", "fp64": True}]
    try:
        import cupy
        for i in range(cupy.cuda.runtime.getDeviceCount()):
            p = cupy.cuda.runtime.getDeviceProperties(i)
            out.append({"backend": "cuda/rocm", "name": p["name"].decode(), "fp64": True})
    except Exception:                                 # noqa: BLE001
        pass
    try:
        import dpctl
        for d in dpctl.get_devices():
            out.append({"backend": f"intel ({d.backend.name}:{d.device_type.name})",
                        "name": d.name, "fp64": bool(d.has_aspect_fp64)})
    except Exception:                                 # noqa: BLE001
        pass
    try:
        import pyopencl as cl
        for p in cl.get_platforms():
            for d in p.get_devices():
                out.append({"backend": f"opencl ({p.name})", "name": d.name.strip(),
                            "fp64": bool(d.double_fp_config)})
    except Exception:                                 # noqa: BLE001
        pass
    return out


def available_backends() -> list[str]:
    out = ["cpu"]
    for name in ("cuda", "intel"):
        try:
            get_backend(name)
            out.append(name)
        except RuntimeError:
            pass
    return out
