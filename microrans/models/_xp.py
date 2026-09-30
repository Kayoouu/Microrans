"""`anp` : fonctions de NumPy aiguillées vers le module du tableau reçu (NumPy, CuPy, dpnp,
faux GPU des tests). CuPy accepte np.maximum(tableau_cupy) (protocole __array_ufunc__),
dpnp non : les modèles de turbulence, communs au 1D (NumPy) et au 2D (tout matériel),
appellent donc anp.maximum, anp.where… au lieu de np.maximum, np.where…"""
from __future__ import annotations

import numpy as np


def _module(values):
    for a in values:
        t = type(a)
        if t is np.ndarray or isinstance(a, (int, float, np.generic)) or a is None:
            continue
        if hasattr(t, "__xp_module__") or t.__module__.startswith(("cupy", "dpnp")):
            from ..linalg import array_module
            return array_module(a)
    return np


class _AnyArray:
    def __getattr__(self, name):
        np_f = getattr(np, name)

        def f(*args, **kw):
            mod = _module(args + tuple(kw.values()))
            return (np_f if mod is np else getattr(mod, name))(*args, **kw)
        f.__name__ = name
        return f


anp = _AnyArray()
