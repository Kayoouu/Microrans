"""Solveurs linéaires creux pour les volumes finis (inspirés d'OpenFOAM et PETSc).

- `AggregationAMG` : multigrille algébrique par agrégation de paires (comme le GAMG
  « faceAreaPair » d'OpenFOAM). La hiérarchie (qui agrège qui) est construite UNE fois à
  partir du graphe du maillage ; à chaque résolution seuls les opérateurs grossiers
  A_c = PᵀAP sont recalculés (sommes creuses, peu coûteuses). Lisseur de Chebyshev
  diagonal (entièrement vectorisé, donc exécutable tel quel sur GPU).
- `pcg` / `pbicgstab` : Krylov préconditionnés, écrits en opérations vectorielles pour
  fonctionner aussi bien avec NumPy/SciPy qu'avec CuPy (même code).

Tolérance : `rtol` relative au résidu INITIAL (‖b − A x0‖₂), comme `relTol` d'OpenFOAM.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


# ------------------------------------------------------------------------ backend
def array_module(x):
    """numpy, cupy (ou module de test) selon le type du tableau."""
    mod = getattr(type(x), "__xp_module__", None)
    if mod is not None:
        return mod
    if type(x).__module__.startswith("cupy"):
        import cupy
        return cupy
    return np


def sparse_module(xp):
    if xp is np:
        return sp
    mod = getattr(xp, "__sparse_module__", None)
    if mod is not None:
        return mod
    import cupyx.scipy.sparse as csp
    return csp


# ------------------------------------------------------------------------ agrégation
def pairwise_aggregation(n: int, i, j, w, rounds: int = 30):
    """Appariement glouton « poignée de main » vectorisé.

    Chaque nœud libre propose à son voisin libre le plus fortement couplé ; les propositions
    mutuelles forment une paire. Les nœuds restés seuls rejoignent l'agrégat de leur voisin
    le plus fort (ou restent singletons). Retourne (agg, nc) avec agg[fin] = indice grossier.
    """
    i = np.asarray(i)
    j = np.asarray(j)
    w = np.asarray(w, dtype=float)
    # départage déterministe des poids égaux (maillages réguliers), symétrique par arête
    h = ((i * 2654435761 + j * 40503) % 1000003) / 1000003.0
    w = w * (1.0 + 1e-6 * h)
    src = np.concatenate([i, j])
    dst = np.concatenate([j, i])
    ww = np.concatenate([w, w])
    match = np.full(n, -1)
    for _ in range(rounds):
        free = (match[src] < 0) & (match[dst] < 0) & (src != dst)
        if not np.any(free):
            break
        s, d, wf = src[free], dst[free], ww[free]
        order = np.lexsort((-wf, s))              # par nœud, poids décroissant
        s, d = s[order], d[order]
        first = np.ones(len(s), dtype=bool)
        first[1:] = s[1:] != s[:-1]
        best = np.full(n, -1)
        best[s[first]] = d[first]
        cand = np.nonzero(best >= 0)[0]
        mutual = cand[best[best[cand]] == cand]
        mutual = mutual[mutual < best[mutual]]    # chaque paire une seule fois
        if len(mutual) == 0:
            break
        match[mutual] = best[mutual]
        match[best[mutual]] = mutual
    agg = np.full(n, -1)
    paired = np.nonzero((match >= 0) & (np.arange(n) < match))[0]
    agg[paired] = np.arange(len(paired))
    agg[match[paired]] = agg[paired]
    nc = len(paired)
    # nœuds isolés : rejoindre l'agrégat du voisin apparié le plus fort
    lone = agg < 0
    if np.any(lone):
        sel = lone[src] & ~lone[dst]
        s, d, wf = src[sel], dst[sel], ww[sel]
        order = np.lexsort((-wf, s))
        s, d = s[order], d[order]
        first = np.ones(len(s), dtype=bool)
        first[1:] = s[1:] != s[:-1]
        agg[s[first]] = agg[d[first]]
        rest = np.nonzero(agg < 0)[0]
        agg[rest] = nc + np.arange(len(rest))
        nc += len(rest)
    return agg, nc


def galerkin(A, agg, nc):
    """A_c = Pᵀ A P pour un prolongateur constant par morceaux (P[i, agg[i]] = 1)."""
    A = A.tocoo()
    xp = array_module(A.data)
    spm = sparse_module(xp)
    aggx = xp.asarray(agg)
    Ac = spm.coo_matrix((A.data, (aggx[A.row], aggx[A.col])), shape=(nc, nc)).tocsr()
    Ac.sum_duplicates()
    return Ac


class AggregationAMG:
    """Hiérarchie d'agrégation fixe (dépend du graphe), opérateurs recalculés par matrice.

    smoother : 'gs' (Gauss-Seidel symétrique, noyau C de pyamg, CPU), 'chebyshev'
    (vectorisé, CPU ou GPU) ou 'auto' (gs si possible).
    """

    def __init__(self, n: int, i, j, w, coarse_size: int = 300, max_levels: int = 25,
                 passes: int = 2, smoother: str = "auto"):
        self.n = n
        self.smoother = smoother
        self.levels = []                     # liste de (agg, nc)
        cur_n, ci, cj, cw = n, np.asarray(i), np.asarray(j), np.asarray(w, dtype=float)
        while cur_n > coarse_size and len(self.levels) < max_levels:
            agg = np.arange(cur_n)
            nc = cur_n
            for _ in range(passes):              # 2 appariements successifs → agrégats de ~4
                a2, nc2 = pairwise_aggregation(nc, agg[ci], agg[cj], cw)
                agg = a2[agg]
                nc = nc2
            if nc > 0.8 * cur_n:
                break                            # agrégation inefficace : on s'arrête
            self.levels.append((agg, nc))
            # graphe grossier : somme des poids des arêtes entre agrégats
            keep = agg[ci] != agg[cj]
            a, b = agg[ci[keep]], agg[cj[keep]]
            lo, hi = np.minimum(a, b), np.maximum(a, b)
            key = lo * nc + hi
            uk, inv = np.unique(key, return_inverse=True)
            cw = np.bincount(inv, cw[keep], len(uk))
            ci, cj = uk // nc, uk % nc
            cur_n = nc
        self._xp_levels = None

    @classmethod
    def from_matrix(cls, A, **kw):
        A = A.tocoo()
        off = A.row < A.col
        return cls(A.shape[0], A.row[off], A.col[off], np.abs(A.data[off]), **kw)

    def setup(self, A, scale_correction: bool = True):
        """Préconditionneur V-cycle pour la matrice A (même graphe que la hiérarchie).

        Les structures creuses des niveaux grossiers sont calculées une fois pour une
        structure fine donnée ; ensuite A_c = PᵀAP n'est qu'une somme indexée (bincount).
        """
        xp = array_module(A.data)
        spm = sparse_module(xp)
        A = spm.csr_matrix(A)
        if not A.has_sorted_indices:
            A.sort_indices()
        key = (xp is np, A.shape[0], A.nnz, int(A.indices[:64].sum()), int(A.indptr[-1]))
        if getattr(self, "_struct_key", None) != key:
            self._build_structure(A, xp, spm)
            self._struct_key = key
        smoother = self.smoother
        if smoother == "auto":
            smoother = "gs" if (xp is np and _pyamg_gs() is not None) else "chebyshev"
        ops = []
        data = A.data
        Al = A
        for (agg, nc), (pos, nnz, indices, indptr) in zip(self._xp_levels[1], self._maps):
            ops.append(_level_data(Al, xp, smoother == "chebyshev"))
            data = xp.bincount(pos, weights=data, minlength=nnz)
            Al = spm.csr_matrix((data, indices, indptr), shape=(nc, nc))
        if xp is np:
            lu = sp.linalg.splu(Al.tocsc())
            csolve = lu.solve
        else:                                    # GPU : inverse dense du niveau grossier
            inv = xp.linalg.inv(Al.toarray())
            csolve = lambda r: inv @ r            # noqa: E731
        return AMGPrecond(ops, self._xp_levels[1], csolve, xp, smoother, scale_correction)

    def _build_structure(self, A, xp, spm):
        self._xp_levels = (xp, [(xp.asarray(agg), nc) for agg, nc in self.levels])
        self._maps = []
        n = A.shape[0]
        rows = xp.repeat(xp.arange(n), xp.diff(A.indptr))
        cols = A.indices
        for agg, nc in self._xp_levels[1]:
            key = agg[rows].astype(xp.int64) * nc + agg[cols]
            uk, pos = xp.unique(key, return_inverse=True)
            indices = (uk % nc).astype(xp.int32)
            indptr = xp.concatenate([xp.zeros(1, dtype=xp.int64),
                                     xp.cumsum(xp.bincount(uk // nc, minlength=nc))])
            self._maps.append((pos.ravel(), len(uk), indices, indptr.astype(xp.int32)))
            rows, cols = uk // nc, indices


def _pyamg_gs():
    """Noyau C++ de Gauss-Seidel de pyamg (appel direct, sans surcoût de vérification)."""
    try:
        from pyamg import amg_core
    except ImportError:
        return None
    return amg_core.gauss_seidel


def _level_data(A, xp, need_lmax=True):
    d = A.diagonal()
    dinv = xp.where(xp.abs(d) > 1e-300, 1.0 / d, 1.0)
    if not need_lmax:
        return A, dinv, 2.0
    # borne de Gershgorin de λmax(D⁻¹A) (≈ 2 pour une M-matrice)
    n = A.shape[0]
    rows = xp.repeat(xp.arange(n), xp.diff(A.indptr))
    rs = xp.bincount(rows, weights=xp.abs(A.data), minlength=n)
    lmax = float(xp.max(rs * xp.abs(dinv)))
    return A, dinv, lmax


class AMGPrecond:
    """V-cycle : lissage avant/après, correction grossière mise à l'échelle.

    La mise à l'échelle e ← (eᵀr / eᵀAe) e de la correction grossière (« scaleCorrection »
    d'OpenFOAM) compense l'agrégation non lissée ; elle rend le préconditionneur non
    linéaire, d'où l'usage du gradient conjugué FLEXIBLE (`fcg`).
    """

    def __init__(self, ops, levels, csolve, xp, smoother="chebyshev", scale=True,
                 degree: int = 2):
        self.ops, self.levels, self.csolve, self.xp = ops, levels, csolve, xp
        self.smoother, self.scale, self.degree = smoother, scale, degree
        self._gs = _pyamg_gs() if smoother == "gs" else None

    def _smooth(self, A, dinv, lmax, b, x):
        if self._gs is not None:           # balayage symétrique avant + arrière
            x = x.copy()
            n = A.shape[0]
            self._gs(A.indptr, A.indices, A.data, x, b, 0, n, 1)
            self._gs(A.indptr, A.indices, A.data, x, b, n - 1, -1, -1)
            return x
        # Chebyshev sur [λmax/8, λmax] (Saad, Iterative Methods, alg. 12.1)
        a, bb = lmax / 8.0, lmax * 1.02
        theta, delta = 0.5 * (bb + a), 0.5 * (bb - a)
        sigma = theta / delta
        rho = 1.0 / sigma
        r = dinv * (b - A @ x)
        d = r / theta
        for _ in range(self.degree):
            x = x + d
            r = r - dinv * (A @ d)
            rho_new = 1.0 / (2.0 * sigma - rho)
            d = rho_new * rho * d + 2.0 * rho_new / delta * r
            rho = rho_new
        return x

    def _vcycle(self, lvl, b):
        if lvl == len(self.ops):
            return self.csolve(b)
        xp = self.xp
        A, dinv, lmax = self.ops[lvl]
        agg, nc = self.levels[lvl]
        x = self._smooth(A, dinv, lmax, b, xp.zeros_like(b))
        r = b - A @ x
        e = self._vcycle(lvl + 1, xp.bincount(agg, weights=r, minlength=nc))[agg]
        if self.scale:
            eAe = float(e @ (A @ e))
            if eAe > 0.0:
                e = e * (float(e @ r) / eAe)
        return self._smooth(A, dinv, lmax, b, x + e)

    def __call__(self, r):
        return self._vcycle(0, r)


class JacobiPrecond:
    def __init__(self, A):
        xp = array_module(A.data)
        d = A.diagonal()
        self.dinv = xp.where(xp.abs(d) > 1e-300, 1.0 / d, 1.0)

    def __call__(self, r):
        return self.dinv * r


# ------------------------------------------------------------------------ Krylov
def pcg(A, b, x0, M, rtol=1e-6, maxiter=500):
    """Gradient conjugué préconditionné. Retourne (x, itérations, convergé)."""
    xp = array_module(b)
    x = x0.copy()
    r = b - A @ x
    r0 = float(xp.linalg.norm(r))
    if r0 == 0.0:
        return x, 0, True
    z = M(r)
    p = z.copy()
    rz = float(r @ z)
    for k in range(1, maxiter + 1):
        Ap = A @ p
        pAp = float(p @ Ap)
        if pAp == 0.0:
            return x, k, False
        alpha = rz / pAp
        x += alpha * p
        r -= alpha * Ap
        if float(xp.linalg.norm(r)) <= rtol * r0:
            return x, k, True
        z = M(r)
        rz_new = float(r @ z)
        p = z + (rz_new / rz) * p
        rz = rz_new
    return x, maxiter, False


def fcg(A, b, x0, M, rtol=1e-6, maxiter=500):
    """Gradient conjugué flexible (préconditionneur variable, ex. AMG mis à l'échelle)."""
    xp = array_module(b)
    x = x0.copy()
    r = b - A @ x
    r0 = float(xp.linalg.norm(r))
    if r0 == 0.0:
        return x, 0, True
    p = M(r)
    Ap = A @ p
    for k in range(1, maxiter + 1):
        pAp = float(p @ Ap)
        if pAp <= 0.0:
            return x, k, False
        alpha = float(p @ r) / pAp
        x += alpha * p
        r -= alpha * Ap
        if float(xp.linalg.norm(r)) <= rtol * r0:
            return x, k, True
        z = M(r)
        Az = A @ z
        beta = -float(Az @ p) / pAp
        p = z + beta * p
        Ap = Az + beta * Ap
    return x, maxiter, False


def pbicgstab(A, b, x0, M, rtol=1e-6, maxiter=500):
    """BiCGStab préconditionné à droite. Retourne (x, itérations, convergé)."""
    xp = array_module(b)
    x = x0.copy()
    r = b - A @ x
    r0n = float(xp.linalg.norm(r))
    if r0n == 0.0:
        return x, 0, True
    rhat = r.copy()
    rho = alpha = omega = 1.0
    v = xp.zeros_like(b)
    p = xp.zeros_like(b)
    for k in range(1, maxiter + 1):
        rho_new = float(rhat @ r)
        if rho_new == 0.0:
            return x, k, False
        beta = (rho_new / rho) * (alpha / omega)
        p = r + beta * (p - omega * v)
        ph = M(p)
        v = A @ ph
        den = float(rhat @ v)
        if den == 0.0:
            return x, k, False
        alpha = rho_new / den
        s = r - alpha * v
        if float(xp.linalg.norm(s)) <= rtol * r0n:
            x += alpha * ph
            return x, k, True
        sh = M(s)
        t = A @ sh
        tt = float(t @ t)
        omega = float(t @ s) / tt if tt > 0 else 0.0
        x += alpha * ph + omega * sh
        r = s - omega * t
        if float(xp.linalg.norm(r)) <= rtol * r0n:
            return x, k, True
        if omega == 0.0:
            return x, k, False
        rho = rho_new
    return x, maxiter, False


# ------------------------------------------------------------------------ façade
SOLVERS = ("auto", "direct", "amg", "bicgstab", "pyamg")


class LinearSolver:
    """Choix du solveur par équation, à la manière de `fvSolution` d'OpenFOAM.

    - direct   : LU creuse (SuperLU ; ordonnancement MMD(AᵀA+A) si symétrique) ;
    - amg      : symétrique → CG flexible + AMG agrégé ; sinon BiCGStab + AMG (sans mise
                 à l'échelle, pour garder un préconditionneur linéaire) ;
    - bicgstab : BiCGStab + Jacobi (idéal pour la quantité de mouvement, relTol ~0.1) ;
    - pyamg    : AMG agrégation lissée de pyamg (si installé), reconstruit à chaque appel ;
    - auto     : symétrique → amg (direct si n < `direct_max` ET tolérance serrée < 1e-3) ;
                 non symétrique → bicgstab.
    La hiérarchie AMG est construite une seule fois à partir du graphe (i, j, w) du maillage.
    """

    def __init__(self, n: int, i, j, w, direct_max: int = 20000):
        self.n, self._graph = n, (i, j, w)
        self.direct_max = direct_max
        self._amg = None
        self.stats: dict = {}

    @property
    def amg(self) -> AggregationAMG:
        if self._amg is None:
            i, j, w = self._graph
            self._amg = AggregationAMG(self.n, i, j, w)
        return self._amg

    def resolve(self, method: str, symmetric: bool, rtol: float = 1e-6) -> str:
        if method == "auto":
            if symmetric:
                return "direct" if (self.n < self.direct_max and rtol < 1e-3) else "amg"
            return "bicgstab"
        if method not in SOLVERS:
            raise ValueError(f"Solveur linéaire inconnu '{method}'. Choix : {SOLVERS}")
        return method

    def solve(self, A, b, x0=None, method: str = "auto", rtol: float = 1e-6,
              maxiter: int = 1000, symmetric: bool = False, tag: str = ""):
        xp = array_module(b)
        method = self.resolve(method, symmetric, rtol)
        x0 = xp.zeros_like(b) if x0 is None else x0
        its = 0
        if method == "direct":
            x = _direct(A, b, symmetric, xp)
        elif method == "pyamg" and xp is not np:
            method = "amg"
            M = self.amg.setup(A, scale_correction=symmetric)
            x, its, _ = (fcg if symmetric else pbicgstab)(A, b, x0, M, rtol, maxiter)
        elif method == "pyamg":
            try:
                import pyamg
                ml = pyamg.smoothed_aggregation_solver(A.tocsr())
                r0 = b - A @ x0
                res: list = []
                x = x0 + ml.solve(r0, tol=rtol, accel="cg" if symmetric else "bicgstab",
                                  maxiter=maxiter, residuals=res)
                its = len(res)
            except ImportError:
                x = _direct(A, b, symmetric, xp)
        else:
            if method == "amg":
                M = self.amg.setup(A, scale_correction=symmetric)
                krylov = fcg if symmetric else pbicgstab
            else:
                M = JacobiPrecond(A)
                krylov = pbicgstab
            x, its, ok = krylov(A, b, x0, M, rtol, maxiter)
            if not ok or not bool(xp.all(xp.isfinite(x))):
                x = _direct(A, b, symmetric, xp)          # filet de sécurité
                method += "→direct"
        st = self.stats.setdefault(tag or method, [0, 0])
        st[0] += 1
        st[1] += its
        return x


def _direct(A, b, symmetric, xp):
    if xp is np:
        lu = sp.linalg.splu(A.tocsc(), permc_spec="MMD_AT_PLUS_A" if symmetric else "COLAMD")
        return lu.solve(b)
    solve = getattr(xp, "__spsolve__", None)
    if solve is None:
        from cupyx.scipy.sparse.linalg import spsolve as solve
    return solve(A.tocsr(), b)
