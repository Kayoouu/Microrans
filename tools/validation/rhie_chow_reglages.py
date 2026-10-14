"""C24 (audit 3) : la solution stationnaire convergée dépend-elle des réglages qui ne
devraient changer que le chemin de convergence (`relax_U`, `pseudo_cfl`) ?

Chaque exemple est calculé avec ses réglages, puis avec relax_U = 0.7, 0.5 et
pseudo_cfl = 20, jusqu'à ce que ses efforts ne bougent plus (variation relative < 1e-8
sur 200 itérations) ou 8 000 itérations. Écart relatif des efforts (norme du vecteur, et
composantes x / y) par rapport au calcul par défaut.

    python tools/validation/rhie_chow_reglages.py <dossier> [exemple …]
    (variable d'environnement VARIANTES : par ex. "relax_U=0.5 pseudo_cfl=20")
"""
import copy
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from microrans.cli import examples_dir  # noqa: E402
from microrans.fv2d.case import run_case  # noqa: E402
from microrans.mesh2d.builder import load_config  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
CASES = sys.argv[2:] or ["cylindre_re20", "cavite_re100", "naca0012_sa", "plaque_plane_sa"]
VARIANTS = ["défaut"] + os.environ.get(
    "VARIANTES", "relax_U=0.7 relax_U=0.5 pseudo_cfl=20").split()


def forces(solver):
    return {k: np.asarray(f["total"], dtype=float).copy() for k, f in solver.forces().items()}


def run(name, variant):
    cfg = copy.deepcopy(load_config(examples_dir() / f"{name}.toml"))
    sc = cfg["solver"]
    sc.update(max_iter=200, monitor_tol=None, tol=1e-30, log_every=100000)
    sc.pop("fmg_levels", None)
    if variant != "défaut":
        k, v = variant.split("=")
        sc[k] = float(v)
    t0 = time.perf_counter()
    _, solver = run_case(cfg, base_dir=examples_dir(),
                         out_dir=OUT / f"{name}_{variant.replace('=', '')}", verbose=False,
                         plot=False, return_solver=True)
    f0 = forces(solver)
    change = np.inf
    while solver.iterations_total < 8000:
        solver.run_steady(max_iter=200, tol=1e-30)
        f1 = forces(solver)
        change = max(float(np.linalg.norm(f1[k] - f0[k]) / max(np.linalg.norm(f1[k]), 1e-300))
                     for k in f1)
        f0 = f1
        if change < 1e-8:
            break
    last = solver.history[-1]
    r = {"iterations": solver.iterations_total, "variation_200it": change,
         "residu_max": max(float(v) for k, v in last.items() if k != "iteration"
                           and not k.startswith("probe")),
         "relax_U": solver.settings.relax_U, "forces": {k: v.tolist() for k, v in f0.items()},
         "temps_s": round(time.perf_counter() - t0, 1)}
    print(f"== {name} {variant} : {json.dumps(r, ensure_ascii=False)}", flush=True)
    return r


res = {}
for name in CASES:
    res[name] = {v: run(name, v) for v in VARIANTS}
    ref = res[name]["défaut"]["forces"]
    for v, r in res[name].items():
        r["ecart"] = {k: {"norme": float(np.linalg.norm(np.subtract(f, ref[k])) /
                                         np.linalg.norm(ref[k])),
                          "x": float(f[0] / ref[k][0] - 1) if ref[k][0] else None,
                          "y": float(f[1] / ref[k][1] - 1) if abs(ref[k][1]) > 1e-9 * np.linalg.norm(ref[k]) else None}
                      for k, f in r["forces"].items()}
        txt = []
        for k, e in r["ecart"].items():
            dx = "-" if e["x"] is None else f"{e['x']:+.3%}"
            dy = "-" if e["y"] is None else f"{e['y']:+.3%}"
            txt.append(f"{k}: |ΔF|/|F| {e['norme']:.2e} Δx {dx} Δy {dy}")
        print(f"   {name:16s} {v:14s} it {r['iterations']:5d}  var {r['variation_200it']:.1e}  "
              + "  ".join(txt), flush=True)
    (OUT / f"rhie_chow_{'_'.join(CASES)}.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
