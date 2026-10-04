"""Campagne 8 : README § 8, limite 17 : « une couche entre deux plans de symétrie redonne
exactement le 2D ». Chaque exemple 2D éligible est calculé N itérations en 2D, puis extrudé
sur une couche (z de 0 à 0.1, faces back / front en symétrie, vecteurs complétés par 0) ;
efforts et vitesse moyenne comparés (écart relatif)."""

import copy
import json
import math
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from microrans.tomlio import dumps, loads  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
args = [a for a in sys.argv[2:] if not a.startswith("N=")]
N = next((int(a[2:]) for a in sys.argv[2:] if a.startswith("N=")), 20)
CASES = ["cavite_re100", "cylindre_re20", "naca0012_sa", "plaque_plane_sa",
         "plaque_plane_loi_de_paroi", "plaque_plane_transition_t3a",
         "convection_naturelle_ra1e5", "melange_deux_courants", "cylindre_re100_urans"]
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")
DZ = 0.1


def to3d(cfg):
    c = copy.deepcopy(cfg)
    m = c["mesh"]
    m["extrude"] = {"z0": 0.0, "z1": DZ, "nz": 1,
                    "patch_types": {"back": "symmetry", "front": "symmetry"}}
    c.setdefault("boundary", {}).update(back={"type": "symmetry"}, front={"type": "symmetry"})
    for spec in c["boundary"].values():
        if isinstance(spec.get("U"), list) and len(spec["U"]) == 2:
            spec["U"] = spec["U"] + [0.0]
    for sec, key in (("initial", "U"), ("physics", "body_force"), ("energy", "gravity")):
        v = c.get(sec, {}).get(key)
        if isinstance(v, list) and len(v) == 2:
            c[sec][key] = v + [0.0]
    o = c.setdefault("output", {})
    o.pop("animate", None)
    if o.get("probes"):
        o["probes"] = [list(p) + [DZ / 2] for p in o["probes"]]
    for ln in o.get("lines", []) or []:
        ln["start"], ln["end"] = list(ln["start"]) + [DZ / 2], list(ln["end"]) + [DZ / 2]
    if isinstance(o.get("moment_center"), list) and len(o["moment_center"]) == 2:
        o["moment_center"] = o["moment_center"] + [0.0]
    return c


def short(cfg):
    s = cfg.setdefault("solver", {})
    if str(s.get("mode", "steady")).lower() in ("transient", "unsteady"):
        s["t_end"] = float(s.get("t_start", 0.0)) + 3 * float(s["dt"])
    else:
        s["max_iter"], s["tol"] = N, 1e-30


def run(cfg, w, ex_dir):
    w.mkdir(parents=True, exist_ok=True)
    cfg = copy.deepcopy(cfg)
    cfg.setdefault("output", {})["directory"] = str(w / "res")
    path = ex_dir / f"_audit_{w.parent.name}_{w.name}.toml"   # même dossier (fichiers lus)
    path.write_text(dumps(cfg), encoding="utf-8")
    try:
        p = subprocess.run([sys.executable, "-m", "microrans", "run2d", str(path), "--no-plot",
                            "-q"], env=env, capture_output=True, text=True)
    finally:
        path.unlink()
    s = json.loads((w / "res" / "summary.json").read_text()) if p.returncode in (0, 1) else {}
    return p.returncode, p.stderr.strip()[-400:], s


for name in CASES:
    if args and name not in args:
        continue
    ex_dir = REPO / "microrans/examples"
    cfg = loads((ex_dir / f"{name}.toml").read_text(encoding="utf-8"))
    short(cfg)
    r = {"ex": name, "N": N}
    rc2, e2, s2 = run(cfg, OUT / name / "2d", ex_dir)
    rc3, e3, s3 = run(to3d(cfg), OUT / name / "3d", ex_dir)
    r["rc"], r["err"] = (rc2, rc3), (e2, e3)
    worst = []
    for k, v in s2.items():
        if isinstance(v, dict) and "Cd" in v and k in s3:
            for c in ("Cd", "Cl", "Cm", "yplus_max", "Nu_mean"):
                if c in v and c in s3[k]:
                    a, b = v[c], s3[k][c]
                    rel = abs(a - b) / max(abs(a), 1e-12)
                    worst.append((rel, f"{k}.{c}", a, b))
    for i in range(2):
        a, b = s2.get("U_mean", [math.nan] * 2)[i], s3.get("U_mean", [math.nan] * 3)[i]
        worst.append((abs(a - b) / max(abs(a), 1e-12), f"U_mean.{i}", a, b))
    worst.sort(reverse=True)
    r["pire"] = worst[:5]
    r["iterations"] = (s2.get("iterations"), s3.get("iterations"))
    with open(OUT / "c8.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(name, r["rc"], "pire écart relatif", f"{worst[0][0]:.2e}" if worst else None,
          worst[0][1] if worst else "", flush=True)
print("fin")
