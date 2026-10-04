"""Campagne 12 : cas dégénérés et paramètres limites (vitesse nulle partout, une seule
maille, Reynolds très grand, plages de polaire ou de balayage vides, inversées ou de pas
nul) : réponse claire, résultat sensé, ou valeur absurde / blocage ?"""

import json
import math
import os
import resource
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")


def cap():
    resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))


def run(tag, args, timeout=180):
    try:
        p = subprocess.run([sys.executable, "-m", "microrans", *args], cwd=OUT, env=env,
                           capture_output=True, text=True, preexec_fn=cap, timeout=timeout)
        rc, so, se = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        rc, so, se = "blocage", "", ""
    lines = [ln for ln in (se + "\n" + so).splitlines() if ln.strip()]
    msg = next((ln for ln in lines if ln.startswith(("Erreur", "ATTENTION"))), lines[-1] if lines else "")
    out = {"essai": tag, "rc": rc, "message": msg[:250]}
    summ = OUT / f"o_{tag}" / "summary.json"
    if summ.is_file():
        s = json.loads(summ.read_text())
        big = {f"{k}.{c}": v[c] for k, v in s.items() if isinstance(v, dict)
               for c in ("Cd", "Cl", "Cm") if c in v and (v[c] is None or not math.isfinite(v[c])
                                                          or abs(v[c]) > 1e6)}
        out["coefficients_absurdes"] = big
        out["U_ref"] = s.get("reference_velocity")
    print(json.dumps(out, ensure_ascii=False), flush=True)
    with open(OUT / "c12.jsonl", "a") as fh:
        fh.write(json.dumps(out, ensure_ascii=False) + "\n")


S5 = ["--no-plot", "--set", "solver.max_iter=20"]
run("vitesse_nulle", ["run2d", "cavite_re100", "-o", "o_vitesse_nulle", *S5,
                      "boundary.lid.U=[0.0,0.0]"])
run("une_maille", ["run2d", "cavite_re100", "-o", "o_une_maille", *S5, "mesh.nx=1", "mesh.ny=1"])
run("pave_une_maille", ["run2d", "cavite_cubique_re100_3d", "-o", "o_pave_une_maille", *S5,
                        "mesh.nx=1", "mesh.ny=1", "mesh.nz=1"])
run("reynolds_1e8", ["run2d", "cavite_re100", "-o", "o_reynolds_1e8", "--no-plot", "--set",
                     "solver.max_iter=300", "physics.nu=1e-8"])
run("polaire_inversee", ["polar", "naca0012_sa", "--alpha", "10", "0", "2", "-o", "o_pi", "--no-plot",
                         "--set", "solver.max_iter=3"])
run("polaire_pas_nul", ["polar", "naca0012_sa", "--alpha", "0", "4", "0", "-o", "o_pn", "--no-plot",
                        "--set", "solver.max_iter=3"], timeout=60)
run("polaire_un_point", ["polar", "naca0012_sa", "--alpha", "2", "2", "1", "-o", "o_p1", "--no-plot",
                         "--set", "solver.max_iter=3"])
run("balayage_pas_negatif", ["sweep", "cavite_re100", "--param", "physics.nu", "--range", "0.01",
                             "0.03", "-0.01", "-o", "o_bn", "--no-plot", "--set", "solver.max_iter=3"],
    timeout=60)
run("balayage_cle_inconnue", ["sweep", "cavite_re100", "--param", "physics.nuu", "--values", "0.01",
                              "0.02", "-o", "o_bk", "--no-plot", "--set", "solver.max_iter=3"])
run("balayage_maillage", ["sweep", "cavite_re100", "--param", "mesh.nx", "--values", "8", "16",
                          "-o", "o_bm", "--no-plot", "--set", "solver.max_iter=3"])
run("instationnaire_t_end_0", ["run2d", "cylindre_re100_urans", "-o", "o_t0", "--no-plot", "--set",
                               "solver.t_end=0.0"])
run("dt_minuscule", ["run2d", "cylindre_re100_urans", "-o", "o_dtmin", "--no-plot", "--set",
                     "solver.dt=1e-12", "solver.t_end=3e-12"])
print("fin")
