"""NACA 0012, α = 4°, Re = 1e6, Spalart-Allmaras (exemple naca0012_sa) : effet de la taille
du domaine (rayon R du champ lointain) sur C_l et C_d, maillage en O et maillage en C.

Maillage de l'exemple (128 mailles sur le profil, 1re maille 2e-5, 32 mailles de sillage en
C) ; R = 30, 100, 300 avec 64, 70, 76 mailles radiales (même raison de progression, 1.22).
Réglages du solveur de l'exemple sauf relax_turb = 0.5 (le maillage en C converge mieux :
résidu de ν̃ 4.8e-5 à 3 000 itérations contre 8.8e-5 avec 0.8), sans arrêt sur les efforts :
3 000 itérations, puis 3 000 de plus ; l'écart entre les deux relevés mesure l'erreur
d'itération.

    python tools/validation/naca_domaine.py <dossier> [R …]
    (variable d'environnement KINDS = "ogrid" ou "cgrid" : un seul des deux maillages)
"""
import copy
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from microrans.cli import examples_dir  # noqa: E402
from microrans.fv2d.case import run_case  # noqa: E402
from microrans.mesh2d.builder import load_config  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
N_RADIAL = {30.0: 64, 100.0: 70, 300.0: 76}
RADII = [float(a) for a in sys.argv[2:]] or list(N_RADIAL)
KINDS = os.environ.get("KINDS", "ogrid cgrid").split()


def run(kind, R):
    c = copy.deepcopy(load_config(examples_dir() / "naca0012_sa.toml"))
    c["mesh"].update(type=kind, farfield_radius=R, n_radial=N_RADIAL[R])
    c["solver"].update(max_iter=3000, monitor_tol=None, tol=1e-12, relax_turb=0.5,
                       log_every=100000)
    t0 = time.perf_counter()
    s, solver = run_case(c, out_dir=OUT / f"{kind}_R{R:g}", verbose=False, plot=False,
                         return_solver=True)
    a = s["airfoil"]
    first = {"Cd": a["Cd"], "Cl": a["Cl"], "Cd_pression": a["Cd_pressure"]}
    f1 = solver.forces()["airfoil"]["total"].copy()
    solver.run_steady(max_iter=3000)
    f2 = solver.forces()["airfoil"]["total"]
    # écoulement amont selon x (l'incidence est celle du profil) : C_d ∝ F_x, C_l ∝ F_y
    second = {"Cd": a["Cd"] * float(f2[0] / f1[0]), "Cl": a["Cl"] * float(f2[1] / f1[1])}
    last = solver.history[-1]
    r = {"cellules": solver.mesh.n_cells, "it3000": first, "it6000": second,
         "residus_6000": {k: float(v) for k, v in last.items() if k != "iteration"},
         "temps_s": round(time.perf_counter() - t0, 1)}
    print(f"== {kind} R = {R:g} : {r}", flush=True)
    return r


res = {kind: {R: run(kind, R) for R in RADII} for kind in KINDS}
(OUT / f"naca_domaine_{'_'.join(KINDS)}.json").write_text(
    json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
