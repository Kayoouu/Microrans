"""Campagne 6 : l'interface calcule-t-elle exactement ce que calcule la ligne de commande ?
Chaque exemple : N itérations (instationnaire : 3 pas) en ligne de commande, puis le même
réglage dans l'interface (ouvrir l'exemple, changer N, Lancer). Les résumés doivent être
identiques (écart relatif < 1e-12) ; tout écart = l'interface a changé le cas."""

import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
only = sys.argv[2:]
os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["MICRORANS_RESULTS"] = str(OUT / "gui")
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from microrans.cli import examples_dir  # noqa: E402
from microrans.gui.app import MainWindow  # noqa: E402
from microrans.tomlio import loads  # noqa: E402

N = 20
IGNORE = {"wall_time_s", "checkpoint", "time_s", "restart", "file", "directory", "backend"}
dialogs = []
for kind in ("warning", "information", "critical", "question"):
    setattr(QMessageBox, kind, staticmethod(
        lambda *a, k=kind: dialogs.append((k, a[1], str(a[2])[:300])) or QMessageBox.Yes))


def flat(d, p=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            if k not in IGNORE:
                out.update(flat(v, f"{p}{k}."))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(flat(v, f"{p}{i}."))
    else:
        out[p[:-1]] = d
    return out


def overrides(cfg):
    s = cfg.get("solver", {})
    mode = str(s.get("mode", "steady")).lower()
    if mode in ("transient", "unsteady") or (cfg.get("physics", {}).get("compressible")
                                             and "t_end" in s):
        dt = float(s.get("dt", 0.0) or 0.0)
        if dt <= 0:
            return None
        return {"solver.t_end": float(s.get("t_start", 0.0)) + 3 * dt}
    return {"solver.max_iter": N}


app = QApplication.instance() or QApplication([])
win = MainWindow()
win.quiet = False
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")
for f in sorted(examples_dir().glob("*.toml")):
    if (only and f.stem not in only) or f.stem.startswith("mesh_"):
        continue
    cfg = loads(f.read_text(encoding="utf-8"))
    ov = overrides(cfg)
    r = {"ex": f.stem, "overrides": ov}
    if ov is None:
        r["note"] = "pas de dt fixe : ignoré"
    else:
        cli_out = OUT / "cli" / f.stem
        sets = [f"{k}={v!r}" for k, v in ov.items()]
        p = subprocess.run([sys.executable, "-m", "microrans", "run2d", f.stem, "-o",
                            str(cli_out), "--no-plot", "-q", "--set", *sets],
                           env=env, capture_output=True, text=True)
        r["cli_rc"] = p.returncode
        dialogs.clear()
        win.errors.clear()
        win.open_case(f)
        for k, v in ov.items():
            sec, key = k.split(".")
            win.cfg.setdefault(sec, {})[key] = v
        win.load_cfg(win.cfg)
        win.generate_mesh()
        while win.thread is not None:
            app.processEvents()
            time.sleep(0.02)
        win.summary = None
        win.run_2d()
        while win.thread is not None:
            app.processEvents()
            time.sleep(0.02)
        r["gui_dialogs"] = list(dialogs)
        r["gui_errors"] = [e[:300] for e in win.errors]
        try:
            a = flat(json.loads((cli_out / "summary.json").read_text()))
            b = flat(json.loads((win.out_dir() / "summary.json").read_text()))
            diffs = []
            for k in sorted(set(a) | set(b)):
                x, y = a.get(k), b.get(k)
                if isinstance(x, (int, float)) and isinstance(y, (int, float)) and not (
                        isinstance(x, bool) or isinstance(y, bool)):
                    if x != y and not math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-300):
                        diffs.append((k, x, y))
                elif x != y:
                    diffs.append((k, x, y))
            r["diffs"] = diffs[:30]
            r["n_diffs"] = len(diffs)
        except Exception as e:
            r["err"] = repr(e)
    with open(OUT / "c6.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False, default=str) + "\n")
    print(f.stem, "écarts", r.get("n_diffs"), r.get("err", ""), r.get("note", ""), flush=True)
print("fin")
