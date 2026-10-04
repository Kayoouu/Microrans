"""Campagne 15 : bouton « Arrêter » de l'interface (pendant un maillage, un calcul, un
balayage), reprise après arrêt, durée d'un calcul complet comparée à la ligne de commande,
mémoire de l'interface au fil d'un long calcul instationnaire."""

import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["MICRORANS_RESULTS"] = str(OUT / "gui")
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from microrans.cli import examples_dir  # noqa: E402
from microrans.gui.app import MainWindow  # noqa: E402
from microrans.gui.widgets import set_combo  # noqa: E402

dialogs = []
for kind in ("warning", "information", "critical", "question"):
    setattr(QMessageBox, kind, staticmethod(
        lambda *a, k=kind: dialogs.append((k, a[1], str(a[2])[:300])) or QMessageBox.Yes))
app = QApplication.instance() or QApplication([])
win = MainWindow()
win.quiet = False
res = {}


def rss():
    with open("/proc/self/status") as f:
        return next(int(ln.split()[1]) // 1024 for ln in f if ln.startswith("VmRSS"))


def wait(stop_after=None, timeout=900):
    t0 = time.time()
    stopped_at = None
    while win.thread is not None and time.time() - t0 < timeout:
        app.processEvents()
        time.sleep(0.02)
        if stop_after is not None and stopped_at is None and time.time() - t0 > stop_after:
            win.stop()
            stopped_at = time.time()
    app.processEvents()
    return round(time.time() - t0, 1), (round(time.time() - stopped_at, 1) if stopped_at else None)


def ex(name, **solver):
    win.open_case(examples_dir() / f"{name}.toml")
    win.cfg.setdefault("solver", {}).update(solver)
    win.load_cfg(win.cfg)


# 1. arrêt pendant un maillage (hybride, ~25 s)
ex("mesh_cylindre_hybride")
dialogs.clear()
win.generate_mesh()
total, after = wait(stop_after=2.0)
res["arret_maillage"] = {"duree_totale_s": total, "delai_apres_arret_s": after,
                         "maillage": getattr(win.mesh, "n_cells", None), "dialogues": list(dialogs)}
# 2. arrêt pendant un calcul, puis Continuer
ex("naca0012_sa", max_iter=3000)
win.generate_mesh()
wait()
dialogs.clear()
win.summary = None
win.run_2d()
total, after = wait(stop_after=6.0)
it_stop = (win.summary or {}).get("iterations")
ck = (win.out_dir() / "checkpoint.npz").is_file()
win.cfg["solver"]["max_iter"] = 20
win.load_cfg(win.cfg)
win.continue_2d()
wait()
res["arret_calcul"] = {"delai_apres_arret_s": after, "iterations_a_l_arret": it_stop,
                       "checkpoint": ck, "apres_continuer": (win.summary or {}).get("iterations"),
                       "mode": (win.summary or {}).get("restart", {}).get("mode"),
                       "dialogues": list(dialogs)}
# 3. arrêt pendant un balayage
ex("cavite_re100", max_iter=3000, tol=1e-12)
win.generate_mesh()
wait()
set_combo(win.sweep_param, "physics.nu")
win.sweep_values.setText("0.01, 0.005, 0.0025, 0.00125")
win.sweep_jobs.setValue(1)
dialogs.clear()
win.run_sweep_gui()
total, after = wait(stop_after=4.0)
res["arret_balayage"] = {"delai_apres_arret_s": after, "points": len(win.sweep_rows),
                         "dialogues": list(dialogs)}
# 4. durée d'un calcul complet : interface contre ligne de commande (cavité, défaut)
ex("cavite_re100")
win.generate_mesh()
wait()
t0 = time.time()
win.run_2d()
wait()
t_gui = round(time.time() - t0, 1)
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")
t0 = time.time()
subprocess.run([sys.executable, "-m", "microrans", "run2d", "cavite_re100", "-o", str(OUT / "cli"),
                "-q"], env=env, capture_output=True)
res["duree_cavite"] = {"interface_s": t_gui, "ligne_de_commande_s": round(time.time() - t0, 1),
                       "iterations": (win.summary or {}).get("iterations")}
# 5. mémoire au fil d'un calcul instationnaire (cylindre, 400 pas)
ex("cylindre_re100_urans", t_end=20.0)
win.cfg["output"].pop("animate", None)
win.load_cfg(win.cfg)
win.generate_mesh()
wait()
m0 = rss()
win.run_2d()
wait()
res["memoire_instationnaire"] = {"avant_Mo": m0, "apres_Mo": rss(),
                                 "pic_processus_Mo": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024}
(OUT / "c15.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
print(json.dumps(res, ensure_ascii=False, indent=1))
