"""Lot G2 (audit 3) : figures. Un test par point d'audit (L8 à L11)."""
import copy
import os

import numpy as np
import pytest

from microrans.cli import examples_dir
from microrans.mesh2d.builder import load_config


def _ex(name):
    return copy.deepcopy(load_config(examples_dir() / f"{name}.toml"))


def _spy_savefig(monkeypatch):
    """{nom du fichier : [libellés des courbes de chaque graphe]} des figures enregistrées."""
    from matplotlib.figure import Figure
    seen = {}
    orig = Figure.savefig

    def savefig(self, fname, *a, **k):
        seen[os.path.basename(str(fname))] = [[ln.get_label() for ln in ax.get_lines()]
                                              for ax in self.axes]
        return orig(self, fname, *a, **k)
    monkeypatch.setattr(Figure, "savefig", savefig)
    return seen


def test_probes_are_not_plotted_as_residuals(tmp_path, monkeypatch, capsys):
    """L8 : les valeurs des sondes (probe1_Ux…) étaient tracées dans convergence.png, sur
    l'axe « résidu normalisé » (lignes plates vers 1, comme des résidus bloqués ; valeurs
    négatives absentes), et mêlées aux résidus dans la ligne de progrès."""
    from microrans.fv2d.case import run_case
    seen = _spy_savefig(monkeypatch)
    c = _ex("cavite_re100")
    c["mesh"].update(nx=8, ny=8)
    c["solver"].update(max_iter=5, log_every=1)
    c["output"]["probes"] = [[0.5, 0.5], [0.25, 0.75]]
    out = tmp_path / "c"
    run_case(c, base_dir=str(examples_dir()), out_dir=out, verbose=True, plot=True)
    conv = [lab for ax in seen["convergence.png"] for lab in ax]
    assert conv and not [lab for lab in conv if "probe" in lab], conv
    assert (out / "probes.png").exists()
    assert seen["probes.png"][0] == ["sonde 1", "sonde 2"]        # un graphe par grandeur
    assert len(seen["probes.png"]) == 3                          # Ux, Uy, p
    lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.startswith("  it ")]
    assert lines
    for ln in lines:
        res, _, pr = ln.partition("| sondes :")
        assert "probe" not in res and "probe1_Ux=" in pr, ln
    head = (out / "history.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "probe2_p" in head                                     # toujours dans le CSV


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


def test_gui_live_curve_has_no_probes(win):
    """L8, interface : la courbe « Convergence » et le texte de progrès mêlaient aussi les
    sondes aux résidus."""
    win.history, win._it0, win.solver = [], None, None
    win._run_meta = (True, 100, 1.0)
    win._plot_t = win._plot_cost = 0.0
    for i in range(1, 6):
        win._on_progress({"iteration": i, "Ux": 1.0 / i, "Uy": 0.5 / i, "p": 2.0 / i,
                          "continuity": 1e-3 / i, "probe1_Ux": -0.2, "probe1_p": 1.5})
    win._plot_history()
    labels = [ln.get_label() for ln in win.canvas.fig.axes[0].get_lines()]
    assert labels == ["Ux", "Uy", "p", "continuity"], labels
    assert "probe" not in win.run_info.text()


def test_profile_round_off_noise_is_not_drawn_as_a_curve(tmp_path, monkeypatch):
    """L9 : line_*.png traçait le bruit d'arrondi à pleine échelle (conduite carrée et canal
    3D : U_y, U_z ≈ 1e-15 dessinés comme des courbes qui varient fortement, seul « 1e-15 »
    en coin d'axe le disait)."""
    from matplotlib.figure import Figure

    from microrans.fv2d.case import run_case
    seen = {}
    orig = Figure.savefig

    def savefig(self, fname, *a, **k):
        if os.path.basename(str(fname)).startswith("line_"):
            seen.update({ax.get_ylabel(): (ax.get_title(), ax.get_lines()[0].get_ydata())
                         for ax in self.axes})
        return orig(self, fname, *a, **k)
    monkeypatch.setattr(Figure, "savefig", savefig)
    c = _ex("conduite_carree_3d")
    c["mesh"].update(ny=8, nz=8)
    c["solver"]["max_iter"] = 20
    run_case(c, base_dir=str(examples_dir()), out_dir=tmp_path / "c", verbose=False)
    for k in ("Uy", "Uz"):
        title, y = seen[k]
        assert title.startswith("nul : bruit d'arrondi"), (k, title)
        assert not y[np.isfinite(y)].any()
    assert seen["Ux"][0] == ""                         # écoulement principal : inchangé
    assert seen["Ux"][1].max() > 1.0


def test_profile_note_thresholds():
    """L9 : règle du bruit (< 1e-9 de la référence) et des petites grandeurs (< 1e-2)."""
    from microrans.fv2d.sampling import profile_note
    v, note = profile_note("Uy", [1e-16, -2e-16, np.nan], u_ref=1.0, p_ref=0.5)
    assert np.array_equal(v[:2], [0.0, 0.0]) and np.isnan(v[2]) and note.startswith("nul")
    v, note = profile_note("p", [0.3 + 1e-12, 0.3 - 1e-12], u_ref=1.0, p_ref=0.5)
    assert v == pytest.approx([0.3, 0.3]) and note.startswith("constant")
    v, note = profile_note("p", [-5e-14, -6e-14], u_ref=20.0, p_ref=200.0)   # canal 3D
    assert list(v) == [0.0, 0.0] and note.startswith("constant")
    v, note = profile_note("Uz", [4e-4, -1e-4], u_ref=1.0, p_ref=0.5)    # petit, réel
    assert list(v) == [4e-4, -1e-4] and note == "petit : max |Uz| = 4.0e-04 |U|max"
    assert profile_note("Ux", [1.0, 0.5], 1.0, 0.5)[1] is None
    assert profile_note("T", [1e-20, 0.0], 1.0, 0.5)[1] is None          # autre grandeur
    assert profile_note("Uy", [1e-16], None, None)[1] is None            # sans référence


def test_pipe_figures_show_the_whole_radius(tmp_path, monkeypatch):
    """L10 : tuyau axisymétrique (rayon 0 à 1) : figures de champs cadrées sur la paroi
    (r de 0.7 à 1.3, x de −0.2 à 1.0), cœur du tuyau caché, moitié du cadre vide : la
    paroi du tuyau était prise pour un corps."""
    import types

    import microrans.fv2d.post as post
    from microrans.fv2d.case import run_case
    from microrans.mesh2d.builder import build_mesh
    zooms = []
    orig = post.plot_field
    monkeypatch.setattr(post, "plot_field",
                        lambda *a, **k: zooms.append(k.get("zoom")) or orig(*a, **k))
    c = _ex("tuyau_turbulent_sst")
    c["solver"]["max_iter"] = 3
    run_case(c, base_dir=str(examples_dir()), out_dir=tmp_path / "t", verbose=False)
    assert zooms and all(z is None for z in zooms), zooms        # vue globale
    for name, box in (("cylindre_re20", (-1.5, 4.5, -2.0, 2.0)),   # corps : cadrage gardé
                      ("sphere_re100_axisym", (-1.5, 4.5, -1.5, 2.0))):
        c = _ex(name)
        m = build_mesh(c, base_dir=str(examples_dir()))
        m.set_patch_types({k: "wall" for k, v in c["boundary"].items()
                           if v.get("type") == "wall" and k in m.patch_types})
        assert post._zoom(types.SimpleNamespace(mesh=m)) == pytest.approx(box, abs=1e-3)


def test_sod_tube_band_has_profiles_not_unreadable_field_figures(tmp_path, monkeypatch):
    """L11 : tube de Sod (bande de 1 × 400 mailles) : mesh.png, Mach.png, p.png, rho.png,
    T.png étaient des lignes d'un pixel, illisibles ; la bonne figure est line_axe.png.
    Profil compressible : U_y (≈ 1e-14) suit la règle de L9."""
    from matplotlib.figure import Figure

    import microrans.fv2d.post as post
    from microrans.fv2d.case import run_case
    from microrans.mesh2d.builder import build_mesh
    titles = {}
    orig = Figure.savefig

    def savefig(self, fname, *a, **k):
        name = os.path.basename(str(fname))
        titles.update({(name, ax.get_ylabel()): ax.get_title() for ax in self.axes})
        return orig(self, fname, *a, **k)
    monkeypatch.setattr(Figure, "savefig", savefig)
    c = _ex("compressible_tube_sod")
    c["solver"]["t_end"] = 0.2 * c["solver"]["t_end"]
    out = tmp_path / "s"
    with pytest.warns(UserWarning, match=r"Figures de champs omises : domaine en bande "
                                         r"\(400 fois plus long que large\)"):
        run_case(c, base_dir=str(examples_dir()), out_dir=out, verbose=False)
    pngs = sorted(p.name for p in out.glob("*.png"))
    assert pngs == ["line_axe.png", "probes.png"], pngs          # sondes : L8
    for fig in ("line_axe.png", "probes.png"):         # U_y ≈ 1e-14 (profil), 1e-22 (sondes)
        assert titles[fig, "Uy"].startswith("nul : bruit d'arrondi"), titles[fig, "Uy"]
        assert titles[fig, "rho"] == ""
    # domaine ordinaire (sang, 1:20) : figures de champs gardées
    m = build_mesh(_ex("sang_carreau_artere"), base_dir=str(examples_dir()))
    assert not post._thin_band(m) and not post._thin_band(m, mirror=1)


def test_gui_profile_round_off_noise(win, tmp_path):
    """L9, interface : le « Profil » de l'interface traçait aussi le bruit d'arrondi à
    pleine échelle (U_y ≈ 1e-15 de la conduite carrée)."""
    from microrans.fv2d.case import run_case
    c = _ex("conduite_carree_3d")
    c["mesh"].update(ny=8, nz=8)
    c["solver"]["max_iter"] = 20
    summary, solver = run_case(c, base_dir=str(examples_dir()), out_dir=tmp_path / "g",
                               verbose=False, plot=False, return_solver=True)
    win.load_cfg(c)                                   # cas 3D : ligne à 3 composantes
    win._run_done((summary, solver))
    titles = {}
    for key in ("Uy", "Ux"):
        win.field_combo.setCurrentIndex(win.field_combo.findData(key))
        win.plot_line()
        titles[key] = win.canvas.fig.axes[0].get_title()
    assert titles["Uy"].startswith("Profil de") and "nul : bruit d'arrondi" in titles["Uy"]
    assert "—" not in titles["Ux"], titles["Ux"]
