"""NACA 0012, α = 4°, Re = 1e6, Spalart-Allmaras (exemple naca0012_sa) : maillage en O
contre maillage en C, trois maillages chacun (rapport 2), C_l, C_d (pression, frottement)
et GCI (Celik et al. 2008).

Niveau k : 64·2^k mailles sur le profil, 32·2^k dans la direction normale, 16·2^k le long
du sillage (en C), première maille 4e-5 / 2^k (niveau 1 = l'exemple, y⁺ ≈ 1). Mêmes
réglages du solveur pour les deux maillages : pas de pseudo-temps local (pseudo_cfl = 20 ;
sans lui, les efforts du maillage en C dérivent encore après 3 000 itérations), arrêt quand
les efforts varient de moins de 1e-5 sur 100 itérations (comme l'exemple ; 1e-6 demande
9 550 itérations au niveau 0 en C, contre 1 972 en O).
Sillage proche : cellules traversées par le segment x = 1.1 (0.1 corde derrière le bord de
fuite), |y| ≤ 0.15, et vitesse minimale |U| de ces cellules.

    python tools/validation/naca_o_c.py <dossier> [niveaux, défaut 2 1 0]
    (variable d'environnement KINDS = "ogrid" ou "cgrid" : un seul des deux maillages)
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
from microrans.fv2d.sampling import locate  # noqa: E402
from microrans.gci import gci  # noqa: E402
from microrans.mesh2d.builder import load_config  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
LEVELS = [int(a) for a in sys.argv[2:]] or [2, 1, 0]


def run(kind, level):
    c = copy.deepcopy(load_config(examples_dir() / "naca0012_sa.toml"))
    k = 2 ** level
    c["mesh"].update(type=kind, n_around=64 * k, n_radial=32 * k, first_height=4e-5 / k)
    if kind == "cgrid":
        c["mesh"]["n_wake"] = 16 * k
    # nonorth_limit = 0.5 : réglage des mesures publiées (15 octobre 2026, ancien défaut) ;
    # sans limiteur (défaut depuis le 16), le pseudo-temps converge bien plus lentement sur
    # ce profil (niveau 1 : plus de 3 000 itérations au lieu de 2 063)
    c["solver"].update(max_iter=30000, tol=1e-9, monitor_tol=1e-5, pseudo_cfl=20.0,
                       log_every=1000, nonorth_limit=0.5)
    t0 = time.perf_counter()
    s, solver = run_case(c, out_dir=OUT / f"{kind}_{level}", verbose=False, plot=False,
                         return_solver=True)
    a = s["airfoil"]
    r = {"cellules": solver.mesh.n_cells, "Cl": a["Cl"], "Cd": a["Cd"],
         "Cd_pression": a.get("Cd_pressure"), "Cd_frottement": a.get("Cd_viscous"),
         "Cm": a.get("Cm"), "iterations": s.get("iterations"), "converge": s.get("converged"),
         "temps_s": round(time.perf_counter() - t0, 1)}
    y = np.linspace(-0.15, 0.15, 20001)
    cells = locate(solver.mesh, np.column_stack([np.full_like(y, 1.1), y]))
    cells = np.unique(cells[cells >= 0])
    r["sillage_cellules"] = int(len(cells))
    r["sillage_U_min"] = float(np.linalg.norm(np.asarray(solver.U)[cells], axis=1).min())
    print(f"== {kind} niveau {level} : {r}", flush=True)
    return r


KINDS = os.environ.get("KINDS", "ogrid cgrid").split()
res = {kind: {lv: run(kind, lv) for lv in LEVELS} for kind in KINDS}
if len(LEVELS) == 3:
    for kind, r in res.items():
        cells = [r[lv]["cellules"] for lv in LEVELS]
        for q in ("Cl", "Cd", "Cd_pression", "Cd_frottement"):
            vals = [r[lv][q] for lv in LEVELS]
            if None in vals:
                continue
            try:
                r[f"gci_{q}"] = gci(vals, cells, 2)
            except ValueError as e:
                r[f"gci_{q}"] = str(e)
            print(f"== {kind} {q} : {vals} → {r[f'gci_{q}']}", flush=True)
(OUT / f"naca_o_c_{'_'.join(KINDS)}.json").write_text(
    json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
