"""Campagne 2 : chaque exemple dans l'interface (hors écran) : ouverture, enregistrement
immédiat (le cas change-t-il ?), maillage, calcul court, tous les champs, toutes les
grandeurs pariétales, un profil, la vue du maillage, coupes x / y / z en 3D, « Continuer ».
Boîtes de dialogue et exceptions enregistrées."""

import copy
import json
import os
import sys
import time
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
os.environ["QT_QPA_PLATFORM"] = "offscreen"
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
os.environ["MICRORANS_RESULTS"] = str(OUT / "res")
only = sys.argv[2:]
from PySide6.QtWidgets import QApplication, QMessageBox, QFileDialog  # noqa: E402
from microrans.gui.app import MainWindow  # noqa: E402
from microrans.cli import examples_dir  # noqa: E402
from microrans.tomlio import loads  # noqa: E402
from microrans.fv2d.validate import check_case  # noqa: E402

dialogs, excs = [], []


def rec(kind):
    def f(*a, **k):
        dialogs.append(
            (
                kind,
                str(a[1]) if len(a) > 1 else "",
                str(a[2])[:400] if len(a) > 2 else "",
            )
        )
        return QMessageBox.Yes if kind == "question" else QMessageBox.Ok

    return f


for kind in ("warning", "information", "critical", "question"):
    setattr(QMessageBox, kind, staticmethod(rec(kind)))
QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: ("", ""))
QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: ("", ""))
sys.excepthook = lambda t, v, tb: excs.append(
    "".join(traceback.format_exception(t, v, tb))[-1500:]
)

app = QApplication.instance() or QApplication([])
win = MainWindow()
win.quiet = False


def wait(timeout=600):
    t0 = time.time()
    while win.thread is not None and time.time() - t0 < timeout:
        app.processEvents()
        time.sleep(0.02)
    app.processEvents()
    return round(time.time() - t0, 1)


def diff(a, b, p=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a:
                out.append(f"+ {p}{k} = {b[k]!r}"[:200])
            elif k not in b:
                out.append(f"- {p}{k}")
            else:
                out += diff(a[k], b[k], f"{p}{k}.")
    elif a != b and not (
        isinstance(a, (int, float))
        and isinstance(b, (int, float))
        and float(a) == float(b)
    ):
        out.append(f"~ {p[:-1]} : {a!r} -> {b!r}"[:200])
    return out


def axes_texts():
    t = []
    for ax in win.canvas.fig.axes:
        t.append(ax.get_title())
        t += [x.get_text() for x in ax.texts]
    return [x for x in t if x]


for f in sorted(examples_dir().glob("*.toml")):
    if only and f.stem not in only:
        continue
    dialogs.clear()
    excs.clear()
    win.errors.clear()
    r = {"ex": f.stem}
    try:
        orig = loads(f.read_text(encoding="utf-8"))
        win.open_case(f)
        app.processEvents()
        r["dlg_ouverture"] = list(dialogs)
        dialogs.clear()
        w = OUT / f.stem
        w.mkdir(exist_ok=True)
        win.case_path = w / "enregistre.toml"
        win.save_case()
        saved = loads((w / "enregistre.toml").read_text(encoding="utf-8"))
        r["diff_enregistrement"] = diff(orig, saved)
        try:
            r["warn_enregistre"] = check_case(copy.deepcopy(saved))
        except Exception as e:
            r["err_enregistre"] = str(e)[:400]
        # calcul court
        cfg = win.cfg
        s = cfg.setdefault("solver", {})
        comp = bool(cfg.get("physics", {}).get("compressible"))
        if str(s.get("mode", "steady")).lower() in ("transient", "unsteady") or (
            comp and s.get("t_end")
        ):
            dt = float(s.get("dt", 1e-3) or 1e-3)
            s["t_end"] = float(s.get("t_start", 0.0)) + 3 * dt
            s["max_steps"] = 3
        s["max_iter"] = min(int(s.get("max_iter", 20) or 20), 15)
        cfg.setdefault("output", {})["directory"] = str(w / "res")
        win.load_cfg(cfg)
        r["dlg_load"] = list(dialogs)
        dialogs.clear()
        t0 = time.time()
        win.generate_mesh()
        wait()
        r["t_maillage"] = round(time.time() - t0, 1)
        r["dlg_maillage"] = list(dialogs)
        dialogs.clear()
        r["cellules"] = getattr(win.mesh, "n_cells", None)
        if win.cfg.get("boundary"):
            t0 = time.time()
            win.summary = None
            win.run_2d()
            wait()
            r["t_calcul"] = round(time.time() - t0, 1)
            r["dlg_calcul"] = list(dialogs)
            dialogs.clear()
            r["iterations"] = (win.summary or {}).get("iterations")
            if win.summary is not None:
                nf = win.field_combo.count()
                for i in range(nf):
                    win.field_combo.setCurrentIndex(i)
                    win.plot_field()
                    bad = [
                        t
                        for t in axes_texts()
                        if "hors" in t or "rreur" in t or "impossible" in t
                    ]
                    if bad:
                        r.setdefault("champs_textes", []).append(
                            (win.field_combo.itemData(i), bad)
                        )
                r["n_champs"] = nf
                nw = win.wall_q.count()
                for i in range(nw):
                    win.wall_q.setCurrentIndex(i)
                    win.plot_wall()
                r["n_parietal"] = nw
                m = win.mesh
                lo, hi = m.points.min(0), m.points.max(0)
                a = lo + 0.1 * (hi - lo)
                b = hi - 0.1 * (hi - lo)
                win.line_start.set_value(tuple(float(x) for x in a))
                win.line_end.set_value(tuple(float(x) for x in b))
                win.plot_line()
                r["profil"] = axes_texts()[:3]
                if getattr(m, "dim", 2) == 3:
                    for axis in ("x", "y", "z"):
                        win.slice_axis.setCurrentIndex(win.slice_axis.findData(axis))
                        win.slice_value.set_value(None)
                        win.vec_check.setChecked(True)
                        win.field_combo.setCurrentIndex(0)
                        win.plot_field()
                        r.setdefault("coupes", []).append(axes_texts()[:2])
                    win.vec_check.setChecked(False)
                    win.slice_axis.setCurrentIndex(0)
                win.draw_mesh()
                r["dlg_trace"] = list(dialogs)
                dialogs.clear()
                # Continuer (stationnaire seulement)
                if str(win.cfg["solver"].get("mode", "steady")).lower() == "steady":
                    it0 = win.summary.get("iterations")
                    win.continue_2d()
                    wait()
                    r["continuer"] = [
                        it0,
                        (win.summary or {}).get("iterations"),
                        (win.summary or {}).get("restart", {}).get("mode"),
                    ]
                    r["dlg_continuer"] = list(dialogs)
                    dialogs.clear()
    except Exception:
        r["exception_script"] = traceback.format_exc()[-1500:]
    r["erreurs"] = list(win.errors)
    r["exceptions_qt"] = list(excs)
    with open(OUT / "c2.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False, default=str) + "\n")
    print(
        f.stem,
        "erreurs",
        len(win.errors),
        "exc",
        len(excs),
        "dlg",
        sum(len(v) for k, v in r.items() if k.startswith("dlg")),
        flush=True,
    )
print("fin")
