"""Limiteur de la correction non orthogonale ([solver] nonorth_limit, ψ) : convergence,
stabilité et efforts sur les cas non orthogonaux, pour plusieurs valeurs de ψ.

Cas (maillages des exemples) : NACA 0012 SA en O (naca0012_sa, 76°) et en C (même cas,
[mesh] type = "cgrid") ; cylindre Re = 20 laminaire sur le maillage hybride de
mesh_cylindre_hybride (55°) et sur des triangles seuls (même fichier, type unstructured) ;
deux profils de mesh_naca_multi (63°), laminaire Re = 100, linearUpwind (divergeait à
l'itération 28 lors de l'audit 1, avant C24 : docs/audit_utilisateur.md, Observations).

Chaque calcul va jusqu'à max_iter sans arrêt sur les efforts (tol = 1e-12) ; relevés :
première itération où tous les résidus passent sous 1e-5 et 1e-6, résidus finaux, efforts
finaux, divergence éventuelle.

    python tools/validation/nonorth_limiteur.py <dossier> [ψ …]    (défaut : 0.5 1.0)
    (variable d'environnement CASES : sous-ensemble, ex. CASES="naca_o cyl_hybride")
"""
import copy
import csv
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from microrans.cli import examples_dir  # noqa: E402
from microrans.fv2d.case import run_case  # noqa: E402
from microrans.mesh2d.builder import build_mesh, load_config  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
LIMITS = [float(a) for a in sys.argv[2:]] or [0.5, 1.0]
EX = examples_dir()


def naca(kind):
    c = copy.deepcopy(load_config(EX / "naca0012_sa.toml"))
    if kind == "cgrid":
        c["mesh"].update(type="cgrid", n_wake=32)
    c["solver"].update(max_iter=3000)
    return c


def cylinder(kind):
    c = copy.deepcopy(load_config(EX / "mesh_cylindre_hybride.toml"))
    c["mesh"]["type"] = kind
    if kind == "unstructured":
        c["mesh"].pop("layers", None)
    c["physics"] = {"reynolds": 20, "reference_velocity": 1.0, "reference_length": 1.0,
                    "model": "laminar"}
    c["initial"] = {"U": [1.0, 0.0]}
    c["boundary"] = {"inlet": {"type": "inlet", "U": [1.0, 0.0]},
                     "outlet": {"type": "outlet"}, "bottom": {"type": "symmetry"},
                     "top": {"type": "symmetry"}, "cylinder": {"type": "wall"}}
    c["solver"] = {"mode": "steady", "max_iter": 2000}
    c["output"] = {"forces": ["cylinder"]}
    return c


def multi():
    c = copy.deepcopy(load_config(EX / "mesh_naca_multi.toml"))
    c["physics"] = {"reynolds": 100, "reference_velocity": 1.0, "reference_length": 1.0,
                    "model": "laminar"}
    c["initial"] = {"U": [1.0, 0.0]}
    c["boundary"] = {"inlet": {"type": "inlet", "U": [1.0, 0.0]},
                     "outlet": {"type": "outlet"}, "bottom": {"type": "symmetry"},
                     "top": {"type": "symmetry"}, "main": {"type": "wall"},
                     "flap": {"type": "wall"}}
    c["solver"] = {"mode": "steady", "max_iter": 1500}
    c["output"] = {"forces": ["main", "flap"]}
    return c


CASES = {"naca_o": lambda: naca("ogrid"), "naca_c": lambda: naca("cgrid"),
         "cyl_hybride": lambda: cylinder("hybrid"),
         "cyl_triangles": lambda: cylinder("unstructured"), "multi": multi}
names = os.environ.get("CASES", " ".join(CASES)).split()

res = {}
for name in names:
    base = CASES[name]()
    t0 = time.perf_counter()
    mesh = build_mesh(base, base_dir=str(EX))
    q = mesh.quality()
    print(f"== {name} : {mesh.n_cells} cellules, non-orthogonalité max "
          f"{q['non_orthogonality_max_deg']:.1f}°, maillage {time.perf_counter() - t0:.0f} s",
          flush=True)
    for lim in LIMITS:
        c = copy.deepcopy(base)
        c["solver"].update(tol=1e-12, monitor_tol=None, nonorth_limit=lim, log_every=500)
        out = OUT / f"{name}_psi{lim:g}"
        t0 = time.perf_counter()
        r = {"psi": lim}
        try:
            s = run_case(c, base_dir=str(EX), out_dir=out, verbose=False, plot=False,
                         mesh=mesh)
            r["iterations"] = s.get("iterations")
            r["efforts"] = {k: {q: v[q] for q in ("Cd", "Cl") if q in v}
                            for k, v in s.items() if isinstance(v, dict) and "Cd" in v}
        except Exception as e:  # divergence : message du solveur
            r["erreur"] = f"{type(e).__name__}: {str(e)[:160]}"
        r["temps_s"] = round(time.perf_counter() - t0, 1)
        hist = out / "history.csv"
        if hist.exists():
            rows = list(csv.DictReader(hist.open(encoding="utf-8")))
            keys = [k for k in rows[0] if k not in ("iteration", "continuity")]
            mx = [max(float(row[k]) for k in keys) for row in rows]
            for thr in (1e-5, 1e-6):
                r[f"it_sous_{thr:g}"] = next((int(rows[i]["iteration"]) for i, m in
                                              enumerate(mx) if m < thr), None)
            r["residus_fin"] = {k: float(rows[-1][k]) for k in keys}
            r["residu_max_fin"] = mx[-1]
            r["residu_max_min"] = min(mx)
        res.setdefault(name, []).append(r)
        print(f"   ψ = {lim:g} : {json.dumps(r, ensure_ascii=False)}", flush=True)
(OUT / f"nonorth_limiteur_{'_'.join(names)}.json").write_text(
    json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
