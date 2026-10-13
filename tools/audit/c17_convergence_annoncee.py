"""Campagne 17 : « Convergé » veut-il dire fini ? Chaque exemple stationnaire est lancé avec
ses propres réglages (comme l'utilisateur), puis prolongé sans arrêt sur les efforts ni sur
les résidus, du même nombre d'itérations (au moins 500, au plus 3 000). Relevés : état
annoncé et raison de l'arrêt (résidus / efforts stabilisés), résidus en fin de calcul et en
fin de prolongation, changement des efforts (|ΔF| / |F| par paroi) et des champs
(max |Δ| / étendue du champ ; pression : / max(étendue, ½ max|U|²)) pendant la
prolongation.

Lecture : un exemple annoncé « convergé » dont les efforts bougent encore nettement
pendant la prolongation est une convergence annoncée trompeuse (résultat faux silencieux).

    python tools/audit/c17_convergence_annoncee.py <dossier> [exemple[:clé=valeur,…] …]

Variante de réglages : `naca0012_sa:solver.pseudo_cfl=20` (valeurs lues en TOML).
"""
import copy
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from microrans.cli import examples_dir  # noqa: E402
from microrans.fv2d.case import run_case  # noqa: E402
from microrans.mesh2d.builder import load_config  # noqa: E402
from microrans.tomlio import loads  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
only = {}
for a in sys.argv[2:]:
    name, _, sets = a.partition(":")
    only.setdefault(name, []).append([kv.split("=", 1) for kv in sets.split(",") if kv])


def host(a):
    return np.asarray(a.get() if hasattr(a, "get") else a, dtype=float)


def snapshot(solver):
    forces = {k: host(f["total"]).copy() for k, f in solver.forces().items()}
    fields = {k: host(v).copy() for k, v in solver.fields().items()}
    return forces, fields


def residuals(solver):
    last = solver.history[-1] if solver.history else {}
    return {k: float(v) for k, v in last.items()
            if k != "iteration" and not k.startswith("probe")}


def variants():
    for f in sorted(examples_dir().glob("*.toml")):
        for sets in (only.get(f.stem, []) if only else [[]]):
            cfg = load_config(f)
            for k, v in sets:
                *path, last = k.split(".")
                d = cfg
                for p in path:
                    d = d.setdefault(p, {})
                d[last] = loads(f"v = {v}")["v"]
            name = f.stem + "".join(f" {k}={v}" for k, v in sets)
            yield f, name, cfg


res = []
for f, name, cfg in variants():
    sc = cfg.get("solver")
    if not sc or sc.get("mode", "steady") != "steady" or "mesh" not in cfg:
        continue
    r = {"ex": name}
    t0 = time.perf_counter()
    try:
        s, solver = run_case(copy.deepcopy(cfg), base_dir=f.parent, out_dir=OUT / name.replace(" ", "_"),
                             verbose=False, plot=False, return_solver=True)
        r.update(annonce="convergé" if s["converged"] else "NON convergé",
                 raison=("efforts" if solver.converged == "forces" else
                         "résidus" if solver.converged else "max_iter"),
                 iterations=int(s["iterations"]), residus=residuals(solver),
                 monitor_tol=getattr(solver.settings, "monitor_tol", None),
                 tol=getattr(solver.settings, "tol", None))
        f1, u1 = snapshot(solver)
        n = int(min(max(solver.iterations, 500), 3000))
        solver.settings.monitor_tol = None
        solver.run_steady(max_iter=n, tol=1e-300)
        f2, u2 = snapshot(solver)
        r["prolongation"] = n
        r["residus_apres"] = residuals(solver)
        r["dF_sur_F"] = {k: float(np.linalg.norm(f2[k] - f1[k]) /
                                  max(np.linalg.norm(f1[k]), 1e-300)) for k in f1}
        r["F_avant"] = {k: f1[k].tolist() for k in f1}
        r["F_apres"] = {k: f2[k].tolist() for k in f2}
        r["d_champs"] = {}
        # pression : rapportée à au moins ½ max|U|² (canal périodique : p presque uniforme,
        # son étendue n'est que du bruit)
        q = 0.5 * float(np.max(np.linalg.norm(u1["U"].reshape(len(u1["U"]), -1), axis=1))) ** 2 \
            if "U" in u1 else 0.0
        for k in u1:
            a, b = u1[k], u2[k]
            if a.shape != b.shape or not a.size:
                continue
            span = max(float(np.ptp(a)), q if k == "p" else 0.0)
            r["d_champs"][k] = float(np.max(np.abs(b - a)) / span) if span > 0 else 0.0
    except Exception as e:                        # relevé, pas d'arrêt de la campagne
        r["erreur"] = f"{type(e).__name__}: {e}"[:300]
    r["s"] = round(time.perf_counter() - t0, 1)
    res.append(r)
    print(json.dumps(r, ensure_ascii=False), flush=True)

(OUT / "c17.json").write_text(json.dumps(res, indent=1, ensure_ascii=False),
                              encoding="utf-8")
print("fin")
