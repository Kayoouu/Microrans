"""Noyaux multi-cœurs (Numba, facultatif) des opérateurs volumes finis les plus coûteux.

Désactivé par défaut ([solver] numba = false) : les opérateurs de fvm.py restent en NumPy
pur. Avec [solver] numba = true (Numba installé, CPU), les sommes « faces → cellules » (np.bincount) deviennent
des boucles centrées sur les cellules (chaque cellule additionne ses propres faces : pas
d'écriture concurrente, résultat indépendant du nombre de fils), l'interpolation aux faces
et le gradient de Green-Gauss sont fusionnés (sans tableaux intermédiaires). Seul l'ordre
des additions change par rapport à np.bincount : écarts d'arrondi (~1e-16 relatif),
résultats identiques d'une exécution à l'autre et quel que soit le nombre de fils.

MESURÉ sur la machine de développement (machine virtuelle, 4 cœurs annoncés), de bout en
bout : 1 fil ×0.97 à ×1.18 selon le cas (cavité 16 384 et 147 456 cellules ×1.13, NACA
0012 SA ×1.18, plaque plane ×1.11, cylindre instationnaire ×0.97) ; 2 et 4 fils PLUS LENTS
(×0.76 et ×0.21 sur la cavité) : chaque itération enchaîne des centaines de petites boucles
parallèles et la synchronisation des fils domine sur cette machine. Gain jugé trop faible
pour être actif par défaut. Sur un vrai processeur le résultat peut différer : le mesurer
avec « microrans bench --numba » ([solver] threads = nombre de fils).
"""
from __future__ import annotations

import os

import numpy as np

try:                                                    # pragma: no cover - dépend du poste
    import numba
    from numba import njit, prange
    AVAILABLE = True
except ImportError:                                     # pragma: no cover
    AVAILABLE = False

    def njit(*a, **k):
        return (lambda f: f) if not a or not callable(a[0]) else a[0]
    prange = range


def set_threads(n: int | None) -> int:
    """Nombre de fils (0 ou None : tous les cœurs). Renvoie le nombre effectif."""
    if not AVAILABLE:
        return 1
    n = int(n or 0)
    top = numba.config.NUMBA_NUM_THREADS
    numba.set_num_threads(top if n <= 0 else max(1, min(n, top)))
    return numba.get_num_threads()


def wanted(enabled: bool = False) -> bool:
    """Chemin Numba utilisé ? (installé et demandé : [solver] numba = true, ou variable
    d'environnement MICRORANS_NUMBA=1 ; MICRORANS_NUMBA=0 l'interdit)."""
    env = os.environ.get("MICRORANS_NUMBA")
    return AVAILABLE and (env == "1" or (bool(enabled) and env != "0"))


def segments(owner_of_entry: np.ndarray, n: int):
    """CSR « cellule → entrées » : (ptr, idx) tels que les entrées de la cellule c soient
    idx[ptr[c]:ptr[c+1]] (ordre croissant conservé)."""
    order = np.argsort(owner_of_entry, kind="stable").astype(np.int64)
    ptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(np.bincount(owner_of_entry, minlength=n), out=ptr[1:])
    return ptr, order


# ------------------------------------------------------------------------------ noyaux
@njit(parallel=True, cache=True)
def interp1(w, phi, P, N, out):
    for f in prange(P.shape[0]):
        out[f] = w[f] * phi[P[f]] + (1.0 - w[f]) * phi[N[f]]


@njit(parallel=True, cache=True)
def interp2(w, phi, P, N, out):
    m = phi.shape[1]
    for f in prange(P.shape[0]):
        a, b, wf = P[f], N[f], w[f]
        for j in range(m):
            out[f, j] = wf * phi[a, j] + (1.0 - wf) * phi[b, j]


@njit(parallel=True, cache=True)
def segsum(ptr, idx, vals, out):
    """out[c] = Σ vals[idx[k]], k ∈ [ptr[c], ptr[c+1])."""
    for c in prange(ptr.shape[0] - 1):
        s = 0.0
        for k in range(ptr[c], ptr[c + 1]):
            s += vals[idx[k]]
        out[c] = s


@njit(parallel=True, cache=True)
def facesum(ptr, faces, sign, fi, fb, ni, out):
    """Σ des flux sortants par cellule (faces internes signées, puis frontières)."""
    for c in prange(ptr.shape[0] - 1):
        s = 0.0
        for k in range(ptr[c], ptr[c + 1]):
            f = faces[k]
            if f < ni:
                s += sign[k] * fi[f]
            else:
                s += fb[f - ni]
        out[c] = s


@njit(parallel=True, cache=True)
def green_gauss(ptr, faces, sign, w, P, N, phi, phib, Si, Sb, V, ni, out):
    """Gradient de Green-Gauss fusionné (interpolation linéaire aux faces internes), 2D ou
    3D (nombre de colonnes de out)."""
    nd = out.shape[1]
    for c in prange(ptr.shape[0] - 1):
        for j in range(nd):
            out[c, j] = 0.0
        for k in range(ptr[c], ptr[c + 1]):
            f = faces[k]
            if f < ni:
                s = sign[k] * (w[f] * phi[P[f]] + (1.0 - w[f]) * phi[N[f]])
                for j in range(nd):
                    out[c, j] += s * Si[f, j]
            else:
                b = f - ni
                for j in range(nd):
                    out[c, j] += phib[b] * Sb[b, j]
        for j in range(nd):
            out[c, j] /= V[c]


@njit(parallel=True, cache=True)
def csr_matvec(indptr, indices, data, x, y):
    for i in prange(indptr.shape[0] - 1):
        s = 0.0
        for k in range(indptr[i], indptr[i + 1]):
            s += data[k] * x[indices[k]]
        y[i] = s
