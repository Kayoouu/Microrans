"""Campagne 9b : un cas calculé sur son propre maillage exporté (.msh, .su2) puis réimporté
([mesh] type = "file", sans patch_types) donne-t-il le même résultat ? (20 itérations)."""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from microrans.mesh2d.builder import build_mesh  # noqa: E402
from microrans.mesh2d.io import write_mesh  # noqa: E402
from microrans.tomlio import dumps, loads  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
EX = REPO / "microrans/examples"
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")


def run(cfg, w):
    w.mkdir(parents=True, exist_ok=True)
    cfg = copy.deepcopy(cfg)
    cfg["solver"].update(max_iter=20, tol=1e-30)
    cfg.setdefault("output", {})["directory"] = str(w / "res")
    (w / "cas.toml").write_text(dumps(cfg), encoding="utf-8")
    p = subprocess.run([sys.executable, "-m", "microrans", "run2d", str(w / "cas.toml"),
                        "--no-plot", "-q"], env=env, capture_output=True, text=True)
    if p.returncode not in (0, 1):
        return {"erreur": p.stderr.strip()[-300:]}
    s = json.loads((w / "res" / "summary.json").read_text())
    return {f"{k}.{c}": v[c] for k, v in s.items() if isinstance(v, dict)
            for c in ("Cd", "Cl", "yplus_max") if c in v}


for name in ("naca0012_sa", "cylindre_re20", "plaque_plane_sa", "plaque_plane_loi_de_paroi"):
    cfg = loads((EX / f"{name}.toml").read_text(encoding="utf-8"))
    m = build_mesh(cfg, base_dir=EX)
    r = {"ex": name, "origine": run(cfg, OUT / name / "origine")}
    for ext in ("msh", "su2"):
        write_mesh(m, OUT / name / f"m.{ext}")
        c = copy.deepcopy(cfg)
        c["mesh"] = {"type": "file", "path": str(OUT / name / f"m.{ext}")}
        c.pop("bodies", None)
        c.pop("domain", None)
        r[ext] = run(c, OUT / name / ext)
        o = r["origine"]
        r[f"ecart_{ext}"] = {k: (o[k], r[ext].get(k)) for k in o
                             if r[ext].get(k) is None or abs(o[k] - r[ext][k]) > 1e-9 * max(abs(o[k]), 1e-12)}
    with open(OUT / "c9b.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(name, "msh :", r["ecart_msh"] or "identique", "| su2 :", r["ecart_su2"] or "identique",
          r["msh"].get("erreur", ""), flush=True)
print("fin")
