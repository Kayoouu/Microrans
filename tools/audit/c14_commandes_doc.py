"""Campagne 14 : chaque commande « microrans … » des blocs de code du README et des guides
(tutoriel, dépannage), copiée telle quelle et lancée dans un dossier vide, dans l'ordre du
texte (une commande peut dépendre de la précédente, ex. --continue) ; code de retour,
durée, dernière ligne utile."""

import json
import os
import re
import resource
import shlex
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
DOCS = [REPO / "README.md", REPO / "docs/tutoriel.md", REPO / "docs/depannage.md"]
SKIP = ("microrans gui", "microrans schemes")    # interface ; étude longue qui écrit dans docs
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")


def cap():
    resource.setrlimit(resource.RLIMIT_AS, (6 << 30, 6 << 30))


for doc in DOCS:
    text = doc.read_text(encoding="utf-8")
    blocks = re.findall(r"```(?:bash|sh|console)?\n(.*?)```", text, re.S)
    cmds = []
    for b in blocks:
        for ln in b.splitlines():
            ln = ln.strip()
            if ln.startswith("$ "):
                ln = ln[2:]
            if ln.startswith("microrans "):
                cmds.append(ln.split("  #")[0].split(" #")[0].strip())
    work = OUT / doc.stem
    work.mkdir(exist_ok=True)
    for c in cmds:
        if c.startswith(SKIP):
            continue
        t0 = time.perf_counter()
        try:
            p = subprocess.run([sys.executable, "-m", "microrans", *shlex.split(c)[1:]], cwd=work,
                               env=env, capture_output=True, text=True, timeout=900,
                               preexec_fn=cap)
            rc, so, se = p.returncode, p.stdout, p.stderr
        except subprocess.TimeoutExpired:
            rc, so, se = "> 900 s", "", ""
        lines = [ln for ln in (se + "\n" + so).splitlines() if ln.strip()]
        msg = next((ln for ln in lines if ln.startswith(("Erreur", "ATTENTION"))),
                   lines[-1] if lines else "")
        r = {"doc": doc.name, "cmd": c, "rc": rc, "s": round(time.perf_counter() - t0, 1),
             "message": msg[:250]}
        with open(OUT / "c14.jsonl", "a") as fh:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        print(f"[{doc.name}] rc={rc} {r['s']:7.1f} s  {c[:90]}  ⇒ {msg[:120]}", flush=True)
print("fin")
