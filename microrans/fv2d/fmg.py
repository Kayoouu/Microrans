"""Démarrage multigrille (« Full Multigrid initialization » de Fluent, grid sequencing).

Le cas est d'abord résolu sur des maillages 2, 4… fois plus grossiers, reconstruits à partir
des paramètres du cas (nombre de mailles ÷ 2 par direction, tailles × 2), puis chaque
solution est interpolée sur le maillage suivant (fv2d/restart.py). SIMPLE propage
l'information d'une maille par itération : sur un maillage grossier, l'écoulement
d'ensemble s'établit en beaucoup moins d'itérations, chacune 4 fois moins chère.

[solver] fmg_levels = n (0 : désactivé) ; fmg_tol : tolérance des niveaux grossiers.
Stationnaire seulement. Mesures (README § 7) : cavité 128² 2.5× plus rapide, cylindre 1.4×,
aucun gain sur la plaque plane et le profil turbulents (convergence dominée par les
équations de turbulence). Maillages lus dans un fichier ou préréglés : non pris en charge
(pas de version grossière à reconstruire).
"""
from __future__ import annotations

import copy
from pathlib import Path


def _half(n, minimum=1):
    return max(minimum, int(round(n / 2)))


def coarsen_config(cfg: dict, level: int) -> dict | None:
    """Copie du cas avec un maillage 2**level fois plus grossier, ou None si le type de
    maillage ne s'y prête pas."""
    c = copy.deepcopy(cfg)
    m = c.get("mesh", {})
    kind = m.get("type", "unstructured").lower()
    if "preset" in m or kind == "file":
        return None
    for _ in range(level):
        if isinstance(m.get("extrude"), dict):            # 3D extrudé : couches ÷ 2 aussi
            m["extrude"]["nz"] = _half(m["extrude"].get("nz", 1))
        if kind == "rectangle":
            m["nx"], m["ny"] = _half(m["nx"], 2), _half(m["ny"], 2)
        elif kind == "box":
            m["nx"], m["ny"], m["nz"] = (_half(m["nx"], 2), _half(m["ny"], 2),
                                         _half(m["nz"], 2))
        elif kind == "blocks":
            for b in m["blocks"]:
                b["cells"] = [_half(n) for n in b["cells"]]
        elif kind == "ogrid":
            na = _half(m.get("n_around", 128), 8)
            m["n_around"] = na + (na % 2)            # pair : la coupe à l'axe reste possible
            m["n_radial"] = _half(m.get("n_radial", 64), 4)
            m["first_height"] = 2.0 * m.get("first_height", 1e-3)
        elif kind in ("unstructured", "hybrid"):
            m["h_max"] = 2.0 * m.get("h_max", 1.0)
            m["h_surface"] = 2.0 * m.get("h_surface", 0.05)
            for r in m.get("refinements", []):
                if "h" in r:
                    r["h"] = 2.0 * r["h"]
            if kind == "hybrid":
                lay = m.setdefault("layers", {})
                lay["n"] = _half(lay.get("n", 10), 2)
                lay["first_height"] = 2.0 * lay.get("first_height", 1e-3)
        else:
            return None
    return c


def fmg_initialize(cfg: dict, levels: int, base_dir=".", out_dir=".", verbose=False,
                   callback=None):
    """Résout les niveaux grossiers (du plus grossier au moins grossier) ; renvoie
    (fichier de reprise du dernier niveau, à interpoler sur le maillage fin, ou None ;
    résumé par niveau)."""
    from .case import build_solver
    from .restart import save_checkpoint

    sc = cfg.get("solver", {})
    tol = float(sc.get("fmg_tol", max(float(sc.get("tol", 1e-6)), 1e-4)))
    prev, report = None, []
    stop = [False]

    def cb(s, n):
        if callback and callback(s, n):
            stop[0] = True
            return True
        return None

    for k in range(levels, 0, -1):
        c = coarsen_config(cfg, k)
        if c is None:
            if verbose:
                print("Démarrage multigrille : maillage lu dans un fichier ou préréglé, "
                      "pas de niveau grossier (ignoré).")
            return None, report
        c.setdefault("initial", {}).pop("restart", None)
        if prev is not None:
            c["initial"]["restart"] = str(prev)
        solver = build_solver(c, base_dir, verbose=False)
        ok = solver.run_steady(tol=tol, callback=cb)
        prev = Path(out_dir) / f"fmg_niveau{k}.npz"
        save_checkpoint(solver, prev)
        report.append({"level": k, "n_cells": solver.mesh.n_cells,
                       "iterations": solver.iterations, "converged": bool(ok),
                       "wall_time_s": round(solver.wall_time, 2)})
        if verbose:
            print(f"  multigrille niveau {k} : {solver.mesh.n_cells} cellules, "
                  f"{solver.iterations} itérations ({solver.wall_time:.1f} s)")
        if stop[0]:
            break
    return prev, report
