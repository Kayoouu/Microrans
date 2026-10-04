"""Campagne 4 : combinaisons 3D dans l'interface (hors écran), actions d'un utilisateur :
extrusion de chaque exemple 2D (case à cocher), pavé 3D × modèles, instationnaire ×
schémas, extrémités en z, passage 3D ↔ 2D après maillage, sondes, balayage, Continuer,
export du maillage. Boîtes de dialogue, erreurs et exceptions enregistrées."""

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
from microrans.gui.widgets import set_combo  # noqa: E402
from microrans.cli import examples_dir  # noqa: E402

dialogs, excs, save_target = [], [], [""]


def rec(kind):
    def f(*a, **k):
        dialogs.append(
            (
                kind,
                str(a[1]) if len(a) > 1 else "",
                str(a[2])[:500] if len(a) > 2 else "",
            )
        )
        return QMessageBox.Yes if kind == "question" else QMessageBox.Ok

    return f


for kind in ("warning", "information", "critical", "question"):
    setattr(QMessageBox, kind, staticmethod(rec(kind)))
QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (save_target[0], ""))
sys.excepthook = lambda t, v, tb: excs.append(
    "".join(traceback.format_exception(t, v, tb))[-1500:]
)
app = QApplication.instance() or QApplication([])
win = MainWindow()
win.quiet = False


def wait(timeout=900):
    t0 = time.time()
    while win.thread is not None and time.time() - t0 < timeout:
        app.processEvents()
        time.sleep(0.02)
    app.processEvents()


def texts():
    out = []
    for ax in win.canvas.fig.axes:
        out.append(ax.get_title())
        out += [t.get_text() for t in ax.texts]
    return [t for t in out if t]


def short(cfg, n=5):
    s = cfg.setdefault("solver", {})
    if str(s.get("mode", "steady")).lower() == "transient":
        dt = float(s.get("dt", 0.01) or 0.01)
        s["t_end"] = float(s.get("t_start", 0)) + 3 * dt
    s["max_iter"] = n


def mesh_run_plot(r, slices=True):
    t0 = time.time()
    win.generate_mesh()
    wait()
    r["t_maillage"] = round(time.time() - t0, 1)
    r["cellules"] = getattr(win.mesh, "n_cells", None)
    r["dim_maillage"] = getattr(win.mesh, "dim", 2)
    t0 = time.time()
    win.summary = None
    win.run_2d()
    wait()
    r["t_calcul"] = round(time.time() - t0, 1)
    r["iterations"] = (win.summary or {}).get("iterations")
    r["dim_resume"] = (win.summary or {}).get("dimension")
    if win.summary is None:
        return
    for i in range(win.field_combo.count()):
        win.field_combo.setCurrentIndex(i)
        win.plot_field()
    if slices and getattr(win.mesh, "dim", 2) == 3:
        for axis in ("x", "y", "z"):
            win.slice_axis.setCurrentIndex(win.slice_axis.findData(axis))
            win.slice_value.set_value(None)
            win.field_combo.setCurrentIndex(0)
            win.vec_check.setChecked(True)
            win.plot_field()
            r.setdefault("coupes", []).append(texts()[:1])
        win.vec_check.setChecked(False)
        win.slice_axis.setCurrentIndex(0)
    for i in range(win.wall_q.count()):
        win.wall_q.setCurrentIndex(i)
        win.plot_wall()
    win.draw_mesh()


def extrude(z1, nz=2, ends=None):
    win.extrude_on.setChecked(True)
    win.ext_z0.set_value(0.0)
    win.ext_z1.set_value(z1)
    win.ext_nz.setValue(nz)
    if ends:
        set_combo(win.ext_ends, ends)
    win._store_forms()


SC = {}


def scenario(f):
    SC[f.__name__] = f
    return f


# --- extrusion de chaque exemple 2D (case à cocher), extrémités par défaut (périodiques)
EX2D = [
    "cavite_re100",
    "cylindre_re20",
    "naca0012_sa",
    "plaque_plane_sa",
    "convection_naturelle_ra1e5",
    "melange_deux_courants",
    "sang_carreau_artere",
    "cylindre_re100_urans",
    "tuyau_turbulent_sst",
    "filtre_poreux_conduite",
    "plaque_plane_transition_t3a",
    "compressible_rampe_mach2",
    "mesh_cylindre_hybride",
]
for name in EX2D:

    def make(name=name):
        def f(r):
            win.open_case(examples_dir() / f"{name}.toml")
            short(win.cfg)
            win.load_cfg(win.cfg)
            extrude(0.2, 2)
            r["U_cl"] = {
                k: v.get("U")
                for k, v in win.cfg.get("boundary", {}).items()
                if isinstance(v, dict) and "U" in v
            }
            r["probes"] = win.cfg.get("output", {}).get("probes")
            mesh_run_plot(r)

        f.__name__ = f"extr_{name}"
        return f

    scenario(make())


@scenario
def extr_symetrie(r):
    win.open_case(examples_dir() / "cylindre_re20.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    extrude(0.2, 1, "symmetry")
    mesh_run_plot(r)
    r["cl_back"] = win.cfg["boundary"].get("back")


@scenario
def extr_parois(r):
    win.open_case(examples_dir() / "cavite_re100.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    extrude(0.5, 4, "wall")
    mesh_run_plot(r)


@scenario
def extr_a_regler(r):
    win.open_case(examples_dir() / "cavite_re100.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    extrude(0.5, 4, "patch")
    mesh_run_plot(r)


# --- pavé 3D : modèles
for model in ("laminar", "sa", "ke", "kw", "sst", "sst_gamma"):

    def make(model=model):
        def f(r):
            win.open_case(examples_dir() / "conduite_carree_3d.toml")
            win.cfg["mesh"].update(ny=8, nz=8)
            short(win.cfg)
            win.load_cfg(win.cfg)
            set_combo(win.model_combo, model)
            win._store_forms()
            mesh_run_plot(r, slices=False)

        f.__name__ = f"box_{model}"
        return f

    scenario(make())


@scenario
def box_sst_loi_paroi(r):
    win.open_case(examples_dir() / "canal_turbulent_3d.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    set_combo(win.model_combo, "sst")
    win._store_forms()
    set_combo(win.wall_treat_combo, "wall_function")
    win._store_forms()
    mesh_run_plot(r, slices=False)


for scheme in ("euler", "backward", "crankNicolson", "rk3"):

    def make(scheme=scheme):
        def f(r):
            win.open_case(examples_dir() / "cavite_cubique_re100_3d.toml")
            win.cfg["mesh"].update(nx=8, ny=8, nz=8)
            win.load_cfg(win.cfg)
            set_combo(win.mode_combo, "transient")
            win._store_forms()
            win.cfg["solver"].update(time_scheme=scheme, dt=0.01, t_end=0.03)
            win.load_cfg(win.cfg)
            mesh_run_plot(r, slices=False)

        f.__name__ = f"box_inst_{scheme}"
        return f

    scenario(make())


@scenario
def box_thermique(r):
    win.open_case(examples_dir() / "cavite_cubique_re100_3d.toml")
    win.cfg["mesh"].update(nx=8, ny=8, nz=8)
    short(win.cfg)
    win.load_cfg(win.cfg)
    win.energy_on.setChecked(True)
    win._store_forms()
    r["energy"] = win.cfg.get("energy")
    mesh_run_plot(r, slices=False)


@scenario
def box_depuis_cylindre(r):
    """exemple 2D avec corps → type « Pavé 3D »"""
    win.open_case(examples_dir() / "cylindre_re20.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    set_combo(win.mesh_type, "box")
    win._store_forms()
    r["mesh"] = {k: v for k, v in win.cfg["mesh"].items() if k != "names"}
    mesh_run_plot(r, slices=False)
    r["bc"] = {k: v.get("type") for k, v in win.cfg.get("boundary", {}).items()}


@scenario
def decoche_apres_maillage(r):
    """3D maillé puis case décochée, « Lancer » sans remailler"""
    win.open_case(examples_dir() / "cavite_re100.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    extrude(0.2, 2)
    win.generate_mesh()
    wait()
    win.extrude_on.setChecked(False)
    win._store_forms()
    win.summary = None
    win.run_2d()
    wait()
    r["iterations"] = (win.summary or {}).get("iterations")
    r["dim_resume"] = (win.summary or {}).get("dimension")
    r["dim_maillage"] = getattr(win.mesh, "dim", 2)


@scenario
def coche_apres_maillage(r):
    """2D maillé puis case cochée, « Lancer » sans remailler"""
    win.open_case(examples_dir() / "cavite_re100.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    win.generate_mesh()
    wait()
    extrude(0.2, 2)
    win.summary = None
    win.run_2d()
    wait()
    r["iterations"] = (win.summary or {}).get("iterations")
    r["dim_resume"] = (win.summary or {}).get("dimension")
    r["dim_maillage"] = getattr(win.mesh, "dim", 2)


@scenario
def change_nz_sans_remailler(r):
    win.open_case(examples_dir() / "cavite_cubique_re100_3d.toml")
    win.cfg["mesh"].update(nx=6, ny=6, nz=6)
    short(win.cfg)
    win.load_cfg(win.cfg)
    win.generate_mesh()
    wait()
    win.cfg["mesh"]["nz"] = 3
    win.load_cfg(win.cfg)
    win.summary = None
    win.run_2d()
    wait()
    r["cellules_calcul"] = (win.summary or {}).get("n_cells") or (
        win.summary or {}
    ).get("cells")
    r["cellules_maillage"] = getattr(win.mesh, "n_cells", None)


@scenario
def sondes_saisies(r):
    win.open_case(examples_dir() / "cavite_cubique_re100_3d.toml")
    win.cfg["mesh"].update(nx=6, ny=6, nz=6)
    short(win.cfg)
    win.load_cfg(win.cfg)
    win.probes_edit.setText("0.5 0.5 0.5 ; 0,25 0,75 0,5")
    win.probes_edit.editingFinished.emit()
    win._store_forms()
    r["probes_cfg"] = win.cfg.get("output", {}).get("probes")
    mesh_run_plot(r, slices=False)
    r["probes_resume"] = {k: v for k, v in (win.summary or {}).items() if "probe" in k}


@scenario
def sondes_2d_puis_extrusion(r):
    win.open_case(examples_dir() / "cylindre_re100_urans.toml")
    short(win.cfg)
    win.load_cfg(win.cfg)
    r["probes_avant"] = win.cfg.get("output", {}).get("probes")
    extrude(0.2, 1, "symmetry")
    r["probes_apres"] = win.cfg.get("output", {}).get("probes")
    mesh_run_plot(r, slices=False)


@scenario
def balayage_3d(r):
    win.open_case(examples_dir() / "conduite_carree_3d.toml")
    win.cfg["mesh"].update(ny=6, nz=6)
    short(win.cfg)
    win.load_cfg(win.cfg)
    win.generate_mesh()
    wait()
    set_combo(win.sweep_param, "physics.nu")
    win.sweep_values.setText("0.01, 0.02")
    win.sweep_jobs.setValue(1)
    win.run_sweep_gui()
    wait()
    r["lignes"] = len(win.sweep_rows)


@scenario
def continuer_3d(r):
    win.open_case(examples_dir() / "conduite_carree_3d.toml")
    win.cfg["mesh"].update(ny=6, nz=6)
    short(win.cfg)
    win.load_cfg(win.cfg)
    mesh_run_plot(r, slices=False)
    win.continue_2d()
    wait()
    r["continuer"] = [
        (win.summary or {}).get("iterations"),
        (win.summary or {}).get("restart", {}).get("mode"),
    ]


@scenario
def export_maillage_3d(r):
    win.open_case(examples_dir() / "conduite_carree_3d.toml")
    win.cfg["mesh"].update(ny=6, nz=6)
    win.load_cfg(win.cfg)
    win.generate_mesh()
    wait()
    for ext in ("vtk", "msh", "su2"):
        save_target[0] = str(OUT / f"export.{ext}")
        win.export_mesh()
        app.processEvents()
        r[f"export_{ext}"] = Path(save_target[0]).is_file()


@scenario
def coupe_valeurs(r):
    win.open_case(examples_dir() / "cavite_cubique_re100_3d.toml")
    win.cfg["mesh"].update(nx=6, ny=6, nz=6)
    short(win.cfg)
    win.load_cfg(win.cfg)
    mesh_run_plot(r, slices=False)
    for axis, val in (
        ("x", 0.0),
        ("x", 1.0),
        ("y", 1.0),
        ("z", 0.5),
        ("z", -3.0),
        ("x", 1e-9),
    ):
        win.slice_axis.setCurrentIndex(win.slice_axis.findData(axis))
        win.slice_value.set_value(val)
        win.field_combo.setCurrentIndex(0)
        win.plot_field()
        r.setdefault("coupes", []).append((axis, val, texts()[:2]))


for name in SC:
    if only and name not in only:
        continue
    dialogs.clear()
    excs.clear()
    win.errors.clear()
    r = {"sc": name}
    try:
        SC[name](r)
    except Exception:
        r["exception_script"] = traceback.format_exc()[-1500:]
    r["dialogues"] = list(dialogs)
    r["erreurs"] = [e[:600] for e in win.errors]
    r["exceptions_qt"] = list(excs)
    with open(OUT / "c4.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False, default=str) + "\n")
    print(
        name,
        "it",
        r.get("iterations"),
        "err",
        len(win.errors),
        "exc",
        len(excs) + ("exception_script" in r),
        "dlg",
        len(dialogs),
        flush=True,
    )
print("fin")
