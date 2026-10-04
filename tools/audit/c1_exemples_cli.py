"""Campagne 1 : la commande d'en-tête de chaque exemple, copiée telle quelle, lancée dans
un dossier vide, un exemple à la fois. Durée, pic mémoire, code de retour, sorties,
summary.json (convergence, NaN, avertissements)."""

import json
import math
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
only = sys.argv[2:]
SUB = "rans urans mesh run2d polar sweep examples gui schemes verify devices bench".split()


def commands(f):
    cmds = []
    for ln in f.read_text(encoding="utf-8").splitlines():
        if not ln.startswith("#"):
            break
        for m in re.finditer(r"microrans\s+(\w+)([^()]*)", ln):
            if m.group(1) in SUB:
                c = ("microrans " + m.group(1) + m.group(2)).split("   ")[0].strip()
                cmds.append(c)
    return cmds


def bad_numbers(x, path=""):
    out = []
    if isinstance(x, dict):
        for k, v in x.items():
            out += bad_numbers(v, f"{path}.{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            out += bad_numbers(v, f"{path}[{i}]")
    elif isinstance(x, float) and not math.isfinite(x):
        out.append(path)
    return out


env = dict(
    os.environ, PYTHONPATH=str(REPO), QT_QPA_PLATFORM="offscreen", MPLBACKEND="Agg"
)
env.pop("MICRORANS_RESULTS", None)
res = []
for f in sorted((REPO / "microrans/examples").glob("*.toml")):
    if only and f.stem not in only:
        continue
    cmds = commands(f)
    if not cmds:
        res.append({"ex": f.stem, "erreur": "aucune commande dans l'en-tête"})
        continue
    cmd = cmds[0]
    work = OUT / f.stem
    work.mkdir(exist_ok=True)
    # MICRORANS_EXE=…/dist/microrans/microrans : même campagne avec l'exécutable
    base = (
        [os.environ["MICRORANS_EXE"]]
        if os.environ.get("MICRORANS_EXE")
        else [sys.executable, "-m", "microrans"]
    )
    argv = base + shlex.split(cmd)[1:]
    t0 = time.perf_counter()
    with open(work / "stdout.txt", "w") as so, open(work / "stderr.txt", "w") as se:
        p = subprocess.Popen(argv, cwd=work, env=env, stdout=so, stderr=se)
        _, status, ru = os.wait4(p.pid, 0)
    dt = time.perf_counter() - t0
    rc = os.waitstatus_to_exitcode(status)
    r = {
        "ex": f.stem,
        "cmd": cmd,
        "rc": rc,
        "s": round(dt, 1),
        "Mo": ru.ru_maxrss // 1024,
    }
    files = sorted(
        str(p.relative_to(work))
        for p in work.rglob("*")
        if p.is_file() and p.name not in ("stdout.txt", "stderr.txt")
    )
    r["n_fichiers"] = len(files)
    sums = [p for p in work.rglob("summary.json")]
    if sums:
        try:
            s = json.loads(sums[0].read_text(), parse_constant=lambda c: float(c))
            r["converged"] = s.get("converged")
            r["warnings"] = s.get("warnings")
            r["non_finis"] = bad_numbers(s)
        except Exception as e:
            r["summary_err"] = repr(e)
    err = (work / "stderr.txt").read_text()
    r["stderr"] = err[-600:]
    out = (work / "stdout.txt").read_text()
    r["stdout_fin"] = out[-900:]
    res.append(r)
    print(
        json.dumps(
            {k: r[k] for k in ("ex", "rc", "s", "Mo", "n_fichiers") if k in r},
            ensure_ascii=False,
        ),
        flush=True,
    )
    with open(OUT / "c1.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")
print("fin")
