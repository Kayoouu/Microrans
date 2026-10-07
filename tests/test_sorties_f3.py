"""Lot F3 (audit 2) : figures, sorties, interface. Un test par point d'audit."""
import copy
import json
import os

import numpy as np
import pytest

from microrans.cli import examples_dir, main
from microrans.mesh2d.builder import load_config


def _ex(name):
    return copy.deepcopy(load_config(examples_dir() / f"{name}.toml"))


def _run(tmp_path, name, *sets):
    out = tmp_path / name
    main(["run2d", name, "-o", str(out), "--no-plot", "-q", "--set", "solver.max_iter=20",
          *sets])
    return out, json.loads((out / "summary.json").read_text(encoding="utf-8"))


def test_wall_csv_integrates_to_the_summary_forces(tmp_path):
    """L7 : les CSV pariétaux n'avaient ni aire ni normale (efforts impossibles à refaire
    hors de microrans) ; le frottement du résumé (vecteur complet) n'était pas celui du
    CSV (composante tangentielle : 0.2 % d'écart sur le cylindre Re 20)."""
    for name, coef in (("cylindre_re20", lambda s: 1 / (0.5 * 1.0 * 1.0)),
                       ("sphere_re100_axisym", lambda s: 1 / (0.5 * np.pi / 4)),
                       ("conduite_carree_3d", lambda s: 1 / (0.5 * s["reference_area"])),
                       ("compressible_plaque_laminaire",
                        lambda s: 1 / (0.5 * s["freestream"]["rho"]
                                       * s["freestream"]["speed"] ** 2))):
        out, s = _run(tmp_path, name)
        p_inf = s["freestream"]["p"] if "freestream" in s else 0.0
        for f in sorted(out.glob("wall_*.csv")):
            patch = f.stem[5:]
            head = f.read_text(encoding="utf-8").splitlines()[0].split(",")
            dim = 3 if "z" in head else 2
            for col in ("p", "area", *(f"n{c}" for c in "xyz"[:dim]),
                        *(f"tau_{c}" for c in "xyz"[:dim])):
                assert col in head, (name, col)
            d = np.genfromtxt(f, delimiter=",", names=True)
            n = np.column_stack([d[f"n{c}"] for c in "xyz"[:dim]])
            tv = np.column_stack([d[f"tau_{c}"] for c in "xyz"[:dim]])
            A = d["area"][:, None]
            Fp = ((d["p"] - p_inf)[:, None] * n * A).sum(axis=0)
            Fv = (tv * A).sum(axis=0)
            k = coef(s)
            assert Fp[0] * k == pytest.approx(s[patch]["Cd_pressure"], rel=1e-10, abs=1e-14)
            assert Fv[0] * k == pytest.approx(s[patch]["Cd_viscous"], rel=1e-10), name
            if dim == 2:                           # tau_w = composante tangentielle
                t = np.column_stack([-n[:, 1], n[:, 0]])
                assert np.allclose(d["tau_w"], np.sum(tv * t, axis=1), rtol=1e-9,
                                   atol=1e-12 * np.abs(tv).max())


def test_convergence_axis_ignores_zero_residuals(tmp_path, monkeypatch):
    """P2 (trouvé en mesurant) : un résidu nul (Uy à la 1re itération de la cavité) était
    tracé à 1e-300 : axe de 1e-314 à 1e14, courbes écrasées en haut de la figure."""
    from matplotlib.figure import Figure
    seen = {}
    orig = Figure.savefig

    def savefig(self, fname, *a, **k):
        if str(fname).endswith("convergence.png"):
            seen["ylim"] = self.axes[0].get_ylim()
        return orig(self, fname, *a, **k)
    monkeypatch.setattr(Figure, "savefig", savefig)
    main(["run2d", "cavite_re100", "-o", str(tmp_path / "c"), "-q", "--set", "mesh.nx=8",
          "mesh.ny=8", "solver.max_iter=5"])
    assert 1e-30 < seen["ylim"][0] < seen["ylim"][1] < 1e3


def test_3d_figures_plane_axis_names_and_width(tmp_path, monkeypatch):
    """L5 : figures 3D toujours dans le plan z médian (conduite : bande de 2 mailles le long
    de l'axe), sans noms d'axes, avec une grande marge blanche (domaine haut)."""
    import matplotlib
    matplotlib.use("Agg")
    from microrans.fv2d.case import run_case
    from microrans.fv2d.validate import check_case
    from microrans.mesh2d import plot as mplot
    titles, sizes, labels = [], [], []
    orig = mplot.plot_field

    def spy(mesh, values, path=None, **kw):
        ax = orig(mesh, values, ax=None, **kw)        # sans fichier : figure gardée
        titles.append(ax.get_title())
        sizes.append(tuple(ax.figure.get_size_inches()))
        labels.append((ax.get_xlabel(), ax.get_ylabel()))
        import matplotlib.pyplot as plt
        plt.close(ax.figure)
        return ax
    c = _ex("conduite_carree_3d")
    c["mesh"].update(ny=8, nz=8)
    c["solver"]["max_iter"] = 3
    c["output"]["slice_axis"] = "x"                   # celui de l'exemple
    import microrans.fv2d.post as post
    monkeypatch.setattr(post, "plot_field", spy)
    run_case(c, out_dir=tmp_path / "x", verbose=False)
    assert titles[0] == "|U| (plan x = 0.25)" and labels[0] == ("y", "z")
    assert "vorticité ω_x" in titles[2]
    titles.clear(), sizes.clear(), labels.clear()
    c["output"].update(slice_axis="z")
    run_case(c, out_dir=tmp_path / "z", verbose=False)
    assert titles[0] == "|U| (plan z = 0)" and labels[0] == ("x", "y")
    assert sizes[0][0] < 6.0                           # domaine 0.5 × 1 : avant 10 × 9
    c["output"].update(slice_axis="w", slice_value="a")
    with pytest.raises(ValueError, match=r"(?s)slice_axis = 'w' inconnu.*slice_value = "
                                         r"'a'"):
        check_case(c)
    assert _ex("conduite_carree_3d")["output"]["slice_axis"] == "x"
    c2 = _ex("cavite_re100")
    c2["output"]["slice_axis"] = "x"
    assert any("sans effet en 2D" in w for w in check_case(c2))


@pytest.fixture(scope="module")
def win():
    pytest.importorskip("PySide6.QtWidgets", exc_type=ImportError)
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication

    from microrans.gui.app import MainWindow
    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    w.quiet = True
    yield w
    w.close()
    app.processEvents()


def test_live_convergence_plot_does_not_hold_the_solver(win, monkeypatch):
    """P2 : interface 2 à 6 fois plus lente que la ligne de commande (cavité : 34 à 41 s
    contre 6.6 s) : un tracé complet de la convergence à chaque progrès (toutes les 0.25 s)
    dans le fil de l'interface, qui garde le verrou de Python pendant le tracé."""
    draws = []
    orig = win.canvas.draw
    monkeypatch.setattr(win.canvas, "draw", lambda *a: draws.append(a) or orig(*a))
    win.history, win._it0, win.solver = [], None, None
    win._run_meta = (True, 100, 1.0)
    win._plot_t = win._plot_cost = 0.0
    for i in range(1, 21):                         # 20 progrès rapprochés
        win._on_progress({"iteration": i, "Ux": 1.0 / i, "Uy": 0.0 if i == 1 else 0.5 / i,
                          "p": 2.0 / i})
    assert len(win.history) == 20
    assert len(draws) <= 2                         # avant : 20 tracés complets
    win._plot_history()                            # dernier état ; résidu nul non tracé
    lo, hi = win.canvas.fig.axes[0].get_ylim()
    assert 1e-3 < lo < hi < 10


def _wait(win, timeout=180):
    import time

    from PySide6.QtWidgets import QApplication
    t0 = time.time()
    while win.thread is not None and time.time() - t0 < timeout:
        QApplication.processEvents()
        time.sleep(0.02)


def test_stop_button_during_meshing(win):
    """U20 : « Arrêter » sans effet pendant le maillage (hybride : mené à son terme 21.8 s
    après la demande)."""
    import time

    from PySide6.QtWidgets import QApplication
    win.open_case(examples_dir() / "mesh_cylindre_hybride.toml")
    before = win.mesh
    win.errors.clear()
    win.generate_mesh()
    t0 = time.time()
    while time.time() - t0 < 1.0:
        QApplication.processEvents()
        time.sleep(0.02)
    assert win.thread is not None                  # maillage en cours (≈ 20 s en tout)
    win.stop()
    t1 = time.time()
    _wait(win)
    delay = time.time() - t1
    assert delay < 3.0, delay
    assert win.mesh is before and not win.errors
    assert win.run_btn.isEnabled() and not win.stop_btn.isEnabled()
    from microrans.stop import check_stop
    check_stop()                                   # hors de l'interface : sans effet


def test_extrusion_checkbox_fits_outputs_and_refuses_2d_options(win):
    """U15 : « Extruder » laissait les sondes et lignes de profil à 2 composantes (refus au
    lancement, après le maillage), acceptait un cas axisymétrique (refus au lancement
    seulement) et n'était pas proposé pour le maillage multi-blocs."""
    from microrans.fv2d.validate import check_case
    win.open_case(examples_dir() / "melange_deux_courants.toml")
    win.probes_edit.setText("5 0.5")
    win.extrude_on.setChecked(True)
    oc = win.cfg["output"]
    assert oc["lines"][0]["start"] == [9.9, 0.0, 0.5] and oc["probes"] == [[5.0, 0.5, 0.5]]
    assert win.probes_edit.text() == "5 0.5 0.5"
    check_case(copy.deepcopy(win.cfg))               # avant : refus (start et end à 3 …)
    win.extrude_on.setChecked(False)
    assert win.cfg["output"]["lines"][0]["end"] == [9.9, 1.0]
    win.errors.clear()
    win.open_case(examples_dir() / "sphere_re100_axisym.toml")
    win.extrude_on.setChecked(True)
    assert not win.extrude_on.isChecked() and "extrude" not in win.cfg["mesh"]
    assert any("Extrusion 3D impossible : axisymétrique" in e for e in win.errors)
    win.errors.clear()
    win.open_case(examples_dir() / "plaque_plane_sa.toml")
    assert win.cfg["mesh"]["type"] == "blocks" and not win.box_extrude.isHidden()
    win.extrude_on.setChecked(True)
    assert win.cfg["mesh"]["extrude"]["nz"] >= 1 and win._dim() == 3
    win.extrude_on.setChecked(False)


def test_3d_results_page_details(win):
    """U19 : profil par défaut (0, 0, 0) → (1, 0, 0) même hors du domaine ; « Ouvrir
    l'animation » actif en 3D ; « Zoom sur les corps » actif mais sans effet en coupe x / y ;
    vue 3D du maillage : légende sur le dessin ; types « patch » / « cyclic » affichés."""
    from PySide6.QtWidgets import QPushButton
    win.open_case(examples_dir() / "conduite_carree_3d.toml")
    win.cfg["mesh"].update(ny=8, nz=8)
    win.load_cfg(win.cfg)
    win.generate_mesh()
    _wait(win)
    assert win.mesh is not None and win.mesh.dim == 3
    assert win.line_start.value() == [0.0, 0.0, 0.0]
    assert win.line_end.value() == [0.5, 0.0, 0.0]          # x ∈ [0, 0.5] (avant : 1)
    anim = [b for b in win.findChildren(QPushButton) if b.text() == "Ouvrir l'animation"]
    assert anim and not anim[0].isEnabled()
    win.slice_axis.setCurrentIndex(win.slice_axis.findData("x"))
    assert not win.zoom_check.isEnabled()
    win.slice_axis.setCurrentIndex(win.slice_axis.findData("z"))
    assert win.zoom_check.isEnabled()
    info = win.mesh_info.text()
    assert "inlet (périodique)" in info and "bottom (paroi)" in info and "cyclic" not in info
    win.draw_mesh()
    win.canvas.canvas.draw()
    ax3d = next(a for a in win.canvas.fig.axes if a.name == "3d")
    leg = ax3d.get_legend()
    fig_h = win.canvas.fig.bbox.height
    top = leg.get_window_extent().y1 / fig_h
    assert top < ax3d.get_position().y0 + 0.02               # légende sous la vue
    win.line_end.set_value((0.3, 0.1, 0.0))                  # ligne saisie : gardée
    win._default_profile_line(win.mesh)
    assert win.line_end.value() == [0.3, 0.1, 0.0]


def test_save_as_elsewhere_keeps_relative_files(win, tmp_path, monkeypatch):
    """U18 : « Enregistrer sous » dans un autre dossier : contour profil_volet.dat (et
    maillage importé, fichier de reprise) introuvables ensuite."""
    from PySide6.QtWidgets import QFileDialog

    from microrans.mesh2d.geometry import shape_from_dict
    win.open_case(examples_dir() / "mesh_naca_multi.toml")
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    target = [tmp_path / "a" / "cas.toml"]
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target[0]), "")))
    win.save_case_as()
    saved = load_config(tmp_path / "a" / "cas.toml")
    body = next(b for b in saved["bodies"] if b.get("type") == "file")
    assert body["path"] == "profil_volet.dat"
    assert (tmp_path / "a" / "profil_volet.dat").read_bytes() == \
        (examples_dir() / "profil_volet.dat").read_bytes()
    shape_from_dict(body, tmp_path / "a")         # avant : fichier de contour introuvable
    (tmp_path / "a" / "res").mkdir()
    (tmp_path / "a" / "res" / "checkpoint.npz").write_bytes(b"x")
    win.restart_path.setText("res/checkpoint.npz")
    target[0] = tmp_path / "b" / "cas2.toml"
    win.save_case_as()
    saved = load_config(tmp_path / "b" / "cas2.toml")
    assert saved["initial"]["restart"] == str((tmp_path / "a" / "res" /
                                               "checkpoint.npz").resolve())
    assert (tmp_path / "b" / "profil_volet.dat").is_file()
    win.restart_path.setText("")
