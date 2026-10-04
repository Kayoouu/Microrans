"""Campagne 7 : reprise exacte (README : « identique au bit près »). Pour chaque exemple :
A = 2N itérations d'un coup ; B = N itérations puis `--continue` pour N de plus (même
dossier). Les champs de fields.vtk et les efforts de summary.json doivent être égaux au
bit près. Instationnaire : 6 pas d'un coup contre 3 + 3."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from microrans.mesh2d.io import read_vtk  # noqa: E402
from microrans.tomlio import loads  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
only = sys.argv[2:]
N = 15
CASES = ["cavite_re100", "naca0012_sa", "plaque_plane_transition_t3a", "plaque_plane_loi_de_paroi",
         "convection_naturelle_ra1e5", "melange_deux_courants", "sang_carreau_artere",
         "filtre_poreux_conduite", "eolienne_disque_actuateur", "sphere_re100_axisym",
         "tuyau_turbulent_sst", "conduite_carree_3d", "canal_turbulent_3d",
         "cavite_cubique_re100_3d", "compressible_rampe_mach2",
         "compressible_naca0012_transsonique", "compressible_plaque_laminaire",
         "cylindre_re100_urans", "compressible_tube_sod"]
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")


def run(name, out, sets, cont=False):
    cmd = [sys.executable, "-m", "microrans", "run2d", name, "-o", str(out), "--no-plot", "-q",
           "--set", *sets] + (["--continue"] if cont else [])
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    return p.returncode, p.stderr.strip()[-300:]


def forces(path):
    s = json.loads(path.read_text())
    return {f"{k}.{c}": v[c] for k, v in s.items() if isinstance(v, dict)
            for c in ("Cd", "Cl", "Cm") if c in v}


for name in CASES:
    if only and name not in only:
        continue
    cfg = loads((REPO / "microrans/examples" / f"{name}.toml").read_text(encoding="utf-8"))
    s = cfg.get("solver", {})
    unsteady = str(s.get("mode", "steady")).lower() in ("transient", "unsteady") or "t_end" in s
    if unsteady and "dt" not in s:
        print(name, "ignoré : pas de dt fixe (pas de temps par CFL, découpage 3 + 3 impossible)")
        continue
    if unsteady:
        dt = float(s["dt"])
        t0 = float(s.get("t_start", 0.0))
        full, half, more = [f"solver.t_end={t0 + 6 * dt!r}"], [f"solver.t_end={t0 + 3 * dt!r}"], \
            [f"solver.t_end={t0 + 6 * dt!r}"]
    else:
        full, half, more = [f"solver.max_iter={2 * N}", "solver.tol=1e-30"], \
            [f"solver.max_iter={N}", "solver.tol=1e-30"], [f"solver.max_iter={N}", "solver.tol=1e-30"]
    a, b = OUT / name / "A", OUT / name / "B"
    r = {"ex": name}
    r["A"] = run(name, a, full)
    r["B1"] = run(name, b, half)
    r["B2"] = run(name, b, more, cont=True)
    try:
        fa, fb = read_vtk(a / "fields.vtk")["cell_data"], read_vtk(b / "fields.vtk")["cell_data"]
        r["champs_differents"] = {k: float(np.max(np.abs(fa[k] - fb[k])))
                                  for k in fa if k in fb and not np.array_equal(fa[k], fb[k])}
        r["champs_manquants"] = sorted(set(fa) ^ set(fb))
        ca, cb = forces(a / "summary.json"), forces(b / "summary.json")
        r["efforts_differents"] = {k: (ca[k], cb.get(k)) for k in ca if ca[k] != cb.get(k)}
        sa, sb = json.loads((a / "summary.json").read_text()), json.loads((b / "summary.json").read_text())
        r["iterations"] = (sa.get("iterations"), sb.get("iterations"), sb.get("restart", {}).get("mode"))
    except Exception as e:
        r["err"] = repr(e)
    with open(OUT / "c7.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(name, "champs≠", list(r.get("champs_differents", {}))[:6], "efforts≠",
          len(r.get("efforts_differents", {})), r.get("iterations"), r.get("err", ""), flush=True)
print("fin")
