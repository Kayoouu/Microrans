"""Convection naturelle en cavité carrée chauffée latéralement (de Vahl Davis 1983,
Pr = 0.71) : Nusselt moyen sur la paroi chaude, u maximal sur la verticale x = 0.5,
v maximal sur l'horizontale y = 0.5 (vitesses en α/L). Tableau du README (§ 6, thermique).
u_max, v_max : valeurs aux centres des mailles, moyenne des deux colonnes (rangées)
encadrant le plan médian, maximum sur les rangées (colonnes) — méthode du tableau
d'origine, retrouvée sur l'ancien code (Ra = 1e5 : 34.845 / 68.409 contre 34.840 / 68.405
publiés ; l'échantillonnage interpolé le long des lignes donne 35.09 / 69.46).

Même mise en données que l'exemple convection_naturelle_ra1e5 (ν = Pr, β g = Ra·Pr,
ΔT = 1) ; maillage n × n, uniforme ou resserré aux parois (rapport de taille g entre le
centre et la paroi, comme l'exemple : g = 3). Résidus à 1e-9.

    python tools/validation/vahl_davis.py <dossier> [Ra:n:g …]   (défaut : le tableau)
"""
import copy
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from microrans.cli import examples_dir  # noqa: E402
from microrans.fv2d.case import run_case  # noqa: E402
from microrans.mesh2d.builder import load_config  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
CASES = sys.argv[2:] or ["1e3:48:3", "1e4:48:3", "1e5:48:3", "1e6:96:3"]
REF = {1e3: (1.118, 3.649, 3.697), 1e4: (2.243, 16.178, 19.617),
       1e5: (4.519, 34.73, 68.59), 1e6: (8.800, 64.63, 219.36)}

res = {}
for case in CASES:
    ra, n, g = case.split(":")
    ra, n, g = float(ra), int(n), float(g)
    cfg = copy.deepcopy(load_config(examples_dir() / "convection_naturelle_ra1e5.toml"))
    cfg["mesh"].update(nx=n, ny=n)
    grading = [[[0.5, 0.5, g], [0.5, 0.5, 1.0 / g]]] * 2 if g != 1.0 else [1.0, 1.0]
    cfg["mesh"]["grading"] = grading
    cfg["energy"]["beta"] = ra * 0.71
    cfg["physics"]["reference_velocity"] = max(1.0, 0.25 * ra ** 0.5)
    cfg["solver"].update(max_iter=20000, tol=1e-9, monitor_tol=None)
    s, solver = run_case(cfg, out_dir=OUT / case.replace(":", "_"), verbose=False,
                         plot=False, return_solver=True)
    _, q = solver.wall_heat_flux("chaud")
    nu = float(np.sum(q * solver.fvm.magSb[solver.patch_slices["chaud"]]))
    cc, U = np.asarray(solver.mesh.cell_centers), np.asarray(solver.U)

    def mid_max(axis, comp):
        # moyenne des deux files de cellules de part et d'autre du plan médian
        x = np.unique(np.round(cc[:, axis], 12))
        lo, hi = np.isclose(cc[:, axis], x[x < 0.5].max()), np.isclose(cc[:, axis], x[x > 0.5].min())
        other = 1 - axis
        a = U[lo][np.argsort(cc[lo, other]), comp]
        b = U[hi][np.argsort(cc[hi, other]), comp]
        return float(np.abs(0.5 * (a + b)).max())
    umax, vmax = mid_max(0, 0), mid_max(1, 1)
    ref = REF.get(ra)
    r = {"Nu": nu, "u_max": umax, "v_max": vmax, "iterations": s["iterations"],
         "converge": s["converged"]}
    if ref:
        r["ecarts"] = [f"{a / b - 1:+.2%}" for a, b in zip((nu, umax, vmax), ref)]
    res[case] = r
    print(f"== Ra {ra:g}, {n}², g = {g:g} : {json.dumps(r, ensure_ascii=False)}", flush=True)
(OUT / "vahl_davis.json").write_text(json.dumps(res, indent=1, ensure_ascii=False),
                                     encoding="utf-8")
