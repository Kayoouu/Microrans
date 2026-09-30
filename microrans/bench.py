"""Mesure CPU contre un autre backend (« microrans bench ») sur des cavités entraînées de
taille croissante : temps par itération SIMPLE et écart des résultats. Sert à décider,
sur SA machine, si une carte graphique accélère réellement le calcul (le gain dépend du
matériel, du pilote et de la taille du maillage : il ne se devine pas).
"""
from __future__ import annotations

import time
import warnings

import numpy as np


def _cavity(n, backend):
    from .fv2d import Settings, Solver2D
    from .mesh2d import rectangle_mesh
    m = rectangle_mesh(0, 1, 0, 1, n, n, names={"left": "w", "right": "w", "bottom": "w",
                                                   "top": "lid"})
    # mêmes solveurs linéaires des deux côtés (AMG + Chebyshev, pas de LU) : comparaison
    # du matériel, pas des algorithmes
    st = Settings(backend=backend, solver_p="amg", solver_U="bicgstab")
    return Solver2D(m, 0.01, {"lid": {"type": "wall", "U": [1.0, 0.0]},
                              "w": {"type": "wall"}}, settings=st)


def run_bench(backend: str, sizes=(64, 128, 256), iters: int = 20, verbose=True):
    rows = []
    warnings.simplefilter("ignore")
    _cavity(8, backend).run_steady(max_iter=3, tol=1e-30)   # compilation des noyaux (JIT)
    for n in sizes:
        res = {}
        for bn in ("cpu", backend):
            s = _cavity(n, bn)
            if bn == "cpu":
                s.fvm.lin.amg.smoother = "chebyshev"              # même lisseur partout
            t = time.perf_counter()
            s.run_steady(max_iter=iters, tol=1e-30)
            res[bn] = (time.perf_counter() - t, s.to_cpu())
        (tc, sc), (tb, sb) = res["cpu"], res[backend]
        diff = float(np.max(np.abs(sc.U - sb.U)))
        row = {"cells": n * n, "cpu_ms_per_it": 1e3 * tc / iters,
               "backend_ms_per_it": 1e3 * tb / iters, "speedup": tc / tb, "max_diff_U": diff}
        rows.append(row)
        if verbose:
            print(f"{n * n:>8d} cellules : CPU {row['cpu_ms_per_it']:8.1f} ms/it   "
                  f"{backend} {row['backend_ms_per_it']:8.1f} ms/it   gain ×{row['speedup']:.2f}"
                  f"   écart {diff:.1e}")
    return rows


def run_numba_bench(sizes=(128, 256), threads=(1, 2, 4), iters: int = 20, verbose=True):
    """NumPy contre noyaux Numba ([solver] numba = true) avec 1, 2, 4… fils, sur des
    cavités : temps par itération SIMPLE (compilation exclue) et écart des résultats."""
    from .fv2d import Settings, Solver2D
    from .fv2d import kernels
    from .mesh2d import rectangle_mesh
    if not kernels.AVAILABLE:
        raise RuntimeError("Numba n'est pas installé (pip install numba).")
    warnings.simplefilter("ignore")
    rows = []
    for n in sizes:
        m = rectangle_mesh(0, 1, 0, 1, n, n, names={"left": "w", "right": "w", "bottom": "w",
                                                       "top": "lid"})
        res = {}
        for label, st in [("NumPy", Settings())] + [
                (f"Numba {t} fil(s)", Settings(numba=True, threads=t)) for t in threads]:
            s = Solver2D(m, 0.01, {"lid": {"type": "wall", "U": [1.0, 0.0]},
                                   "w": {"type": "wall"}}, settings=st)
            s.run_steady(max_iter=2, tol=1e-30)             # compilation hors mesure
            t0 = time.perf_counter()
            s.run_steady(max_iter=iters, tol=1e-30)
            res[label] = (1e3 * (time.perf_counter() - t0) / iters, s)
        ref_t, ref_s = res["NumPy"]
        for label, (ms, s) in res.items():
            row = {"cells": n * n, "mode": label, "ms_per_it": ms, "speedup": ref_t / ms,
                   "max_diff_U": float(np.max(np.abs(s.U - ref_s.U)))}
            rows.append(row)
            if verbose:
                print(f"{n * n:>8d} cellules  {label:<16s} {ms:8.1f} ms/it   gain ×{ref_t / ms:.2f}"
                      f"   écart {row['max_diff_U']:.1e}")
    return rows
