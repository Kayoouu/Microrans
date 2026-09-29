"""« Faux GPU » pour tester le chemin CuPy sans carte graphique.

Les tableaux `DeviceArray` enveloppent des tableaux NumPy mais se comportent comme CuPy
sur les points qui cassent un portage GPU :
- conversion implicite vers NumPy interdite (np.asarray(x) → TypeError) ;
- mélange avec un tableau NumPy (hors scalaires) interdit dans les calculs ;
- indexation par un tableau NumPy interdite.
Les fonctions de création du module (zeros, asarray...) acceptent des données CPU, comme
cupy.asarray (transfert explicite).
"""
from __future__ import annotations

import types

import numpy as np
import numpy.lib.mixins
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CREATION = {"asarray", "array", "zeros", "ones", "full", "empty", "arange", "zeros_like",
            "ones_like", "full_like", "empty_like", "linspace", "eye"}


class MixingError(TypeError):
    pass


def _host(x):
    if isinstance(x, DeviceArray):
        return x._a
    if isinstance(x, (list, tuple)):
        return type(x)(_host(v) for v in x)
    if isinstance(x, dict):
        return {k: _host(v) for k, v in x.items()}
    return x


def _check(x, where):
    if isinstance(x, np.ndarray):
        raise MixingError(f"tableau CPU (numpy {x.shape}) mélangé à des tableaux GPU dans "
                          f"{where}")
    if isinstance(x, (list, tuple)):
        for v in x:
            _check(v, where)
    if isinstance(x, dict):
        for v in x.values():
            _check(v, where)


def _wrap(x):
    if isinstance(x, np.ndarray):
        return DeviceArray(x)
    if isinstance(x, tuple):
        return tuple(_wrap(v) for v in x)
    if isinstance(x, list):
        return [_wrap(v) for v in x]
    if isinstance(x, (np.generic,)):
        return DeviceArray(np.asarray(x))       # CuPy : réductions → tableaux 0-d
    return x


class DeviceArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    __slots__ = ("_a",)

    def __init__(self, a):
        self._a = a

    # --- protocole NumPy -------------------------------------------------------
    def __array__(self, dtype=None, copy=None):
        raise MixingError("conversion implicite GPU → NumPy (utiliser .get())")

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kw):
        _check(inputs, ufunc.__name__)
        if out is not None:
            _check(out, ufunc.__name__)
            kw["out"] = tuple(_host(o) for o in out)
        res = getattr(ufunc, method)(*_host(list(inputs)), **kw)
        if out is not None:
            return out[0] if len(out) == 1 else out
        return _wrap(res)

    def __array_function__(self, func, types_, args, kwargs):
        _check(args, func.__name__)
        _check(kwargs, func.__name__)
        return _wrap(func(*_host(list(args)), **_host(kwargs)))

    # --- interface tableau -----------------------------------------------------
    def _idx(self, key):
        if isinstance(key, tuple):
            for k in key:
                if isinstance(k, np.ndarray):
                    raise MixingError("indexation d'un tableau GPU par un tableau NumPy")
            return tuple(_host(k) for k in key)
        if isinstance(key, np.ndarray):
            raise MixingError("indexation d'un tableau GPU par un tableau NumPy")
        return _host(key)

    def __getitem__(self, key):
        return _wrap(self._a[self._idx(key)])

    def __setitem__(self, key, value):
        _check(value, "affectation")
        self._a[self._idx(key)] = _host(value)

    def __len__(self):
        return len(self._a)

    def __iter__(self):
        return (_wrap(v) for v in self._a)

    def __float__(self):
        return float(self._a)

    def __int__(self):
        return int(self._a)

    def __bool__(self):
        return bool(self._a)

    def __index__(self):
        return int(self._a)

    def __format__(self, spec):
        return format(self._a.item() if self._a.ndim == 0 else self._a, spec)

    def __repr__(self):
        return f"DeviceArray({self._a!r})"

    def get(self):
        return self._a.copy()

    def __getattr__(self, name):
        if name.startswith("__"):                # pas de __array_interface__ & co
            raise AttributeError(name)
        attr = getattr(self._a, name)
        if callable(attr):
            def method(*a, **k):
                _check(a, name)
                return _wrap(attr(*_host(list(a)), **_host(k)))
            return method
        return _wrap(attr)


# --------------------------------------------------------------------------- module creux
class FakeCSR:
    def __init__(self, arg, shape=None):
        if isinstance(arg, FakeCSR):
            self._m = arg._m.copy()
        elif isinstance(arg, tuple):
            _check(arg, "csr_matrix")
            self._m = sp.csr_matrix(tuple(_host(list(arg))) if not isinstance(arg[1], tuple)
                                    else (_host(arg[0]), tuple(_host(list(arg[1])))),
                                    shape=shape)
        elif sp.issparse(arg):
            raise MixingError("matrice creuse CPU passée au module GPU")
        else:
            raise TypeError(type(arg))

    @property
    def shape(self):
        return self._m.shape

    @property
    def nnz(self):
        return self._m.nnz

    data = property(lambda self: DeviceArray(self._m.data))
    indices = property(lambda self: DeviceArray(self._m.indices))
    indptr = property(lambda self: DeviceArray(self._m.indptr))
    has_sorted_indices = property(lambda self: self._m.has_sorted_indices)

    def sort_indices(self):
        self._m.sort_indices()

    def diagonal(self):
        return DeviceArray(self._m.diagonal())

    def toarray(self):
        return DeviceArray(self._m.toarray())

    def tocsr(self):
        return self

    def __matmul__(self, x):
        if not isinstance(x, DeviceArray):
            raise MixingError("produit matrice GPU × vecteur CPU")
        return DeviceArray(self._m @ x._a)


def _coo(arg, shape=None):
    raise NotImplementedError("coo_matrix non utilisé sur GPU")


def _spsolve(A, b):
    if not isinstance(b, DeviceArray):
        raise MixingError("spsolve : second membre CPU")
    return DeviceArray(spla.spsolve(A._m.tocsc(), b._a))


# --------------------------------------------------------------------------- module xp
class _FakeXP(types.ModuleType):
    def __getattr__(self, name):
        if name == "linalg":
            return _LINALG
        attr = getattr(np, name)
        if isinstance(attr, type) or not callable(attr):
            return attr
        creation = name in CREATION

        def fn(*a, **k):
            if not creation:
                _check(a, name)
                _check(k, name)
            return _wrap(attr(*_host(list(a)), **_host(k)))
        fn.__name__ = name
        return fn


fakecupy = _FakeXP("fakecupy")
fake_sparse = types.SimpleNamespace(csr_matrix=FakeCSR, coo_matrix=_coo)
fakecupy.__sparse_module__ = fake_sparse
fakecupy.__spsolve__ = _spsolve
_LINALG = types.SimpleNamespace(
    norm=lambda x, *a, **k: _wrap(np.linalg.norm(_host(x), *a, **k)) if isinstance(
        x, DeviceArray) else (_ for _ in ()).throw(MixingError("norm sur tableau CPU")),
    inv=lambda x: _wrap(np.linalg.inv(_host(x))))
DeviceArray.__xp_module__ = fakecupy


def register_fake_gpu():
    from microrans.backend import Backend, register
    be = Backend("fakegpu", fakecupy, fake_sparse)
    register("fakegpu", be)
    return be
