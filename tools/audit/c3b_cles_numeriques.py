"""Campagne 3b (argument « contexte » : 3c, chaque clé dans un cas où elle sert) : chaque clé numérique de [solver] / [flow] / [physics] reçoit un texte, puis
-1 : erreur claire avant calcul, avertissement, accepté en silence, ou erreur interne ?"""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from microrans.tomlio import dumps, loads  # noqa: E402

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)


def ex(n):
    return loads(
        (REPO / "microrans/examples" / f"{n}.toml").read_text(encoding="utf-8")
    )


INC_SOLVER = [
    "relax_U",
    "relax_p",
    "relax_turb",
    "relax_T",
    "relax_scalar",
    "pseudo_dt",
    "pseudo_cfl",
    "max_co",
    "max_dt",
    "cn_theta",
    "ddt_phi_coeff",
    "threads",
    "n_outer",
    "n_corr",
    "n_nonorth",
    "nonorth_limit",
    "max_iter",
    "tol",
    "monitor_tol",
    "monitor_window",
    "dt",
    "t_end",
    "fmg_levels",
]
COMP_SOLVER = [
    "venkat_k",
    "limiter_freeze",
    "entropy_fix",
    "cfl",
    "cfl_max",
    "cfl_growth",
    "cfl_cuts",
    "linear_sweeps",
    "linear_iter",
    "linear_tol",
    "first_order_iter",
    "viscous_factor",
    "max_iter",
    "tol",
    "order",
]
cases = []
if len(sys.argv) > 2 and sys.argv[2] == "contexte":
    cyl = ex("cylindre_re100_urans")
    cyl["solver"]["t_end"] = cyl["solver"].get("dt", 0.05) * 3
    for k in ["cn_theta", "ddt_phi_coeff", "dt", "max_co", "max_dt", "nonorth_limit"]:
        c = copy.deepcopy(cyl)
        if k in ("max_co", "max_dt"):
            c["solver"]["adjust_dt"] = True
        if k == "cn_theta":
            c["solver"]["time_scheme"] = "crankNicolson"
        cases.append(("inst", "solver", k, c))
    naca = ex("compressible_naca0012_transsonique")
    naca["solver"]["max_iter"] = 3
    naca["mesh"]["n_around"] = 48
    naca["mesh"]["n_radial"] = 16
    for k in ["cfl_growth", "cfl_cuts", "linear_sweeps", "linear_iter", "linear_tol"]:
        cases.append(("impl", "solver", k, naca))
    pl = ex("compressible_plaque_laminaire")
    pl["solver"]["max_iter"] = 3
    cases.append(("ns", "solver", "viscous_factor", pl))
cav = ex("cavite_re100")
cav["solver"]["max_iter"] = 3
ramp = ex("compressible_rampe_mach2")
ramp["solver"]["max_iter"] = 3
for k in [] if len(sys.argv) > 2 else INC_SOLVER:
    cases.append(("inc", "solver", k, cav))
for k in [] if len(sys.argv) > 2 else COMP_SOLVER:
    cases.append(("comp", "solver", k, ramp))
for k in (
    []
    if len(sys.argv) > 2
    else ["mach", "reynolds", "alpha", "T", "p", "gamma", "prandtl", "mu"]
):
    cases.append(("comp", "flow", k, ramp))
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")
for kind, sec, k, base in cases:
    for val in ("abc", -1):
        c = copy.deepcopy(base)
        c.setdefault(sec, {})[k] = val
        if sec == "solver" and k == "max_iter" and val == "abc":
            pass
        w = OUT / f"{kind}_{sec}_{k}_{val}"
        w.mkdir(exist_ok=True)
        c.setdefault("output", {})["directory"] = str(w / "res")
        (w / "cas.toml").write_text(dumps(c), encoding="utf-8")
        p = subprocess.run(
            [sys.executable, "-m", "microrans", "run2d", "cas.toml", "--no-plot", "-q"],
            cwd=w,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        err = p.stderr.strip().splitlines() or [""]
        line = next(
            (ln for ln in err if ln.startswith(("Erreur", "ATTENTION", "Remarque"))),
            err[-1],
        )
        cls = (
            "INTERNE"
            if "interne" in p.stderr
            else "ERREUR"
            if p.returncode == 2
            else "AVERT"
            if "ATTENTION" in p.stderr or "Remarque" in p.stderr
            else "SILENCE"
        )
        r = {
            "kind": kind,
            "sec": sec,
            "key": k,
            "val": val,
            "rc": p.returncode,
            "cls": cls,
            "msg": line[:260],
        }
        with open(OUT / "c3b.jsonl", "a") as fh:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        print(
            f"{cls:8s} {kind} {sec}.{k}={val!r} rc={p.returncode} : {line[:150]}",
            flush=True,
        )
print("fin")
