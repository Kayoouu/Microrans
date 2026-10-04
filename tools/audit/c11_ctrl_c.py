"""Campagne 11 : Ctrl-C pendant un calcul en ligne de commande (cavité stationnaire,
cylindre instationnaire, NACA compressible). Le README promet un checkpoint « à l'arrêt
demandé » : qu'est-ce qui est réellement écrit, peut-on reprendre ?"""

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")
LONG = ["solver.max_iter=100000", "solver.tol=1e-30"]
CASES = (("cavite_re100", LONG), ("cylindre_re100_urans", []),
         ("compressible_naca0012_transsonique", LONG))

res = []
for name, sets in CASES:
    out = OUT / f"o_{name}"
    cmd = [sys.executable, "-m", "microrans", "run2d", name, "-o", str(out), "--no-plot"]
    p = subprocess.Popen(cmd + (["--set", *sets] if sets else []), cwd=OUT, env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    time.sleep(12)
    t0 = time.time()
    p.send_signal(signal.SIGINT)
    so, se = p.communicate(timeout=120)
    files = sorted(os.listdir(out)) if out.is_dir() else []
    r = {"cas": name, "rc": p.returncode, "arret_s": round(time.time() - t0, 1),
         "fichiers": files, "fin_sortie": (so + se).strip().splitlines()[-1:]}
    if "checkpoint.npz" in files:
        q = subprocess.run(cmd + ["-q", "--continue", "--set", "solver.max_iter=5"], cwd=OUT,
                           env=env, capture_output=True, text=True)
        r["reprise_rc"] = q.returncode
    res.append(r)
    print(json.dumps(r, ensure_ascii=False))
(OUT / "c11.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
