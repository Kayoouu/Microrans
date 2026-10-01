"""Interface : formulaires ↔ cas TOML (viscosité non newtonienne, scalaires, colonne
« scalaires » des conditions limites), hors écran."""
import os

import pytest

pytest.importorskip("PySide6.QtWidgets", exc_type=ImportError)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def win():
    from PySide6.QtWidgets import QApplication

    from microrans.gui.app import MainWindow
    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    yield w
    w.close()
    app.processEvents()


def test_viscosity_and_scalars_roundtrip(win):
    from microrans.gui.widgets import set_combo
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 2, "y0": 0, "y1": 1, "nx": 8, "ny": 4,
                    "names": {"left": "inlet", "right": "outlet", "bottom": "wall",
                              "top": "wall"}},
           "physics": {"nu": 0.01, "viscosity": {"model": "carreau", "nu0": 1.0,
                                                 "nu_inf": 0.01, "lambda": 2.0, "n": 0.4}},
           "scalars": {"c": {"diffusivity": 1e-3, "scheme": "upwind"}},
           "boundary": {"inlet": {"type": "inlet", "flow_rate": 1.0, "profile": "parabolic",
                                  "scalars": {"c": 1.0}},
                        "outlet": {"type": "outlet"}, "wall": {"type": "wall"}}}
    win.load_cfg(cfg)
    assert win.visc_fields["lambda"].isEnabled() and not win.visc_fields["tau_y"].isEnabled()
    # changement de loi : paramètres de l'ancienne loi retirés du cas
    set_combo(win.visc_combo, "power_law")
    win.visc_fields["K"].set_value(0.1)
    win.visc_fields["n"].set_value(0.5)
    win._visc_changed()
    v = win.cfg["physics"]["viscosity"]
    assert v["model"] == "power_law" and v["K"] == 0.1 and "nu0" not in v
    # retour au newtonien : section supprimée
    set_combo(win.visc_combo, "newtonian")
    win._visc_changed()
    assert "viscosity" not in win.cfg["physics"]
    # scalaires : édition du tableau, clés sans colonne (scheme) conservées
    win.scalar_table.item(0, 3).setText("1")
    assert win.cfg["scalars"]["c"] == {"scheme": "upwind", "diffusivity": 1e-3, "source": 1.0}
    win._scalar_add()
    assert list(win.cfg["scalars"]) == ["c", "c1"]
    # conditions limites : colonne scalaires, profil de débit conservé
    win.bc_table.item(0, 8).setText("c=0.5; c1=2")
    b = win.cfg["boundary"]["inlet"]
    assert b["scalars"] == {"c": 0.5, "c1": 2.0} and b["profile"] == "parabolic"


def test_porous_table_roundtrip(win):
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 2, "y0": 0, "y1": 1, "nx": 8, "ny": 4},
           "physics": {"nu": 0.01},
           "porous": [{"name": "filtre", "region": "circle", "center": [1.0, 0.5],
                       "radius": 0.2, "permeability": 0.01, "forchheimer": [2.0, 4.0],
                       "angle": 30.0}],
           "boundary": {}}
    win.load_cfg(cfg)
    assert win.porous_table.item(0, 1).text() == "cercle 1 0.5 0.2"
    assert win.porous_table.item(0, 2).text() == "100"
    win.porous_table.item(0, 1).setText("rect 0 1 0 0.5")
    z = win.cfg["porous"][0]
    assert z == {"name": "filtre", "region": "rectangle", "x0": 0.0, "x1": 1.0, "y0": 0.0,
                 "y1": 0.5, "darcy": 100.0, "forchheimer": [2.0, 4.0], "angle": 30.0}
    win.porous_table.item(0, 1).setText("x > 1.5")
    assert win.cfg["porous"][0]["expression"] == "x > 1.5"


def test_swirl_and_actuator_disk_forms(win):
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 2, "y0": 0, "y1": 1, "nx": 8, "ny": 4,
                    "names": {"left": "inlet", "right": "outlet", "bottom": "axis",
                              "top": "wall"}},
           "physics": {"nu": 0.01, "axisymmetric": True, "swirl": True},
           "actuator_disk": [{"name": "rotor", "x0": 0.9, "x1": 1.1, "radius": 0.5,
                              "thrust_coefficient": 0.4, "mode": "propeller",
                              "distribution": "optimal"}],
           "boundary": {"inlet": {"type": "inlet", "U": [1.0, 0.0], "U_theta": "2*y"},
                        "outlet": {"type": "outlet"}, "axis": {"type": "axis"},
                        "wall": {"type": "wall", "omega": 3.0}}}
    win.load_cfg(cfg)
    assert win.disk_table.item(0, 5).text() == "CT=0.4 helice"
    win.disk_table.item(0, 6).setText("0.1")
    d = win.cfg["actuator_disk"][0]
    assert d["torque"] == 0.1 and d["distribution"] == "optimal" and d["mode"] == "propeller"
    rows = {win.bc_table.item(r, 0).text(): r for r in range(win.bc_table.rowCount())}
    assert win.bc_table.item(rows["wall"], 9).text() == "Ω=3.0"
    win.bc_table.item(rows["wall"], 9).setText("Ω=5")
    assert win.cfg["boundary"]["wall"]["omega"] == 5.0
    assert win.cfg["boundary"]["inlet"]["U_theta"] == "2*y"


def test_transition_example_roundtrip(win):
    """Cas T3A : modèle sst_gamma et convection de la turbulence conservés par les formulaires."""
    from microrans.cli import examples_dir
    from microrans.mesh2d.builder import load_config
    win.load_cfg(load_config(examples_dir() / "plaque_plane_transition_t3a.toml"))
    assert win.model_combo.currentData() == "sst_gamma"
    win._store_forms()
    assert win.cfg["physics"]["model"] == "sst_gamma"
    assert win.cfg["solver"]["convection_turb"] == "linearUpwindLimited"


def test_coupled_algorithm_roundtrip(win):
    from microrans.cli import examples_dir
    from microrans.mesh2d.builder import load_config
    cfg = load_config(examples_dir() / "cavite_re100.toml")
    cfg["solver"]["algorithm"] = "coupled"
    win.load_cfg(cfg)
    win._store_forms()
    assert win.cfg["solver"]["algorithm"] == "coupled"


def test_numbers_typed_with_dot_or_comma_in_french_locale(win):
    """Audit C1 : avec Windows en français, QDoubleValidator supprimait le point en
    silence (« 0.5 » tapé → 5). Le point et la virgule sont acceptés quelle que soit la
    langue du système."""
    from PySide6.QtCore import QLocale
    from PySide6.QtTest import QTest

    from microrans.gui.widgets import SciEdit
    old = QLocale()
    QLocale.setDefault(QLocale(QLocale.French, QLocale.France))
    try:
        for typed, value in (("0.5", 0.5), ("0,5", 0.5), ("-0.25", -0.25),
                             ("2.5E-4", 2.5e-4), ("1,5e-3", 1.5e-3), ("1 000", 1000.0)):
            e = SciEdit(1.0)
            e.clear()
            QTest.keyClicks(e, typed)
            assert e.value() == pytest.approx(value), (typed, e.text())
    finally:
        QLocale.setDefault(old)


def test_probe_points_roundtrip(win):
    """Audit C2 : les sondes [[x, y], ...] d'un cas devenaient le texte « [[3.0, ... »
    (4 exemples plantaient au lancement dans l'interface)."""
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 2, "y0": 0, "y1": 1, "nx": 8, "ny": 4},
           "physics": {"nu": 0.01},
           "boundary": {"left": {"type": "wall"}, "right": {"type": "wall"},
                        "bottom": {"type": "wall"}, "top": {"type": "wall"}},
           "output": {"probes": [[3.0, 0.0], [5.0, 1.5]]}}
    win.load_cfg(cfg)
    assert win.probes_edit.text() == "3 0 ; 5 1.5"
    win._store_forms()
    assert win.cfg["output"]["probes"] == [[3.0, 0.0], [5.0, 1.5]]
    win.probes_edit.setText("0.5 0.25; 1 0.75")
    win._store_forms()
    assert win.cfg["output"]["probes"] == [[0.5, 0.25], [1.0, 0.75]]
    win.probes_edit.setText("")
    win._store_forms()
    assert "probes" not in win.cfg["output"]


def test_switch_to_triangles_from_rectangle_case_meshes(win):
    """Audit C3 : passer un cas sans [domain] en « Triangles » écrivait [domain] sans
    type → KeyError 'type' à la génération (signalé par un utilisateur)."""
    from microrans.mesh2d.builder import build_mesh
    win.new_case()
    combo = win.mesh_type
    combo.setCurrentIndex(combo.findData("unstructured"))      # comme un clic
    win._store_forms()
    cfg = win.cfg
    assert "type" not in cfg.get("domain", {}) or cfg["domain"]["type"] == "rectangle"
    cfg["mesh"].update(h_max=4.0, h_surface=0.4)
    mesh = build_mesh(cfg)
    names = {p.name for p in mesh.patches}
    assert {"inlet", "outlet", "cylinder"} <= names


def test_error_dialog_text_is_readable():
    """Audit M3 : l'interface affichait l'exception brute (« KeyError : 'type' »)."""
    from microrans.gui.app import user_message
    assert user_message("ValueError : Le maillage en O demande exactement un corps.\n\nTB") \
        .startswith("Le maillage en O demande exactement un corps.")
    assert "paramètre manquant" in user_message("KeyError : 'path'\n\nTB").lower()
    m = user_message("TypeError : unsupported operand\n\nTB")
    assert "erreur interne" in m.lower() and "Journal" in m
    assert "mémoire insuffisante" in user_message("MemoryError : \n\nTB").lower()


def test_mesh_failure_dialog(win, monkeypatch):
    """Audit M3 : un échec du maillage s'affichait comme une exception brute, sous le
    titre « Erreur »."""
    import time

    from PySide6.QtWidgets import QApplication, QMessageBox
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[1:3]))
    monkeypatch.setattr(win, "quiet", False)
    win.new_case()
    win.mesh_type.setCurrentIndex(win.mesh_type.findData("unstructured"))
    win.cfg["bodies"][0]["radius"] = -0.5
    win._load_bodies()
    win.generate_mesh()
    t0 = time.time()
    while win.thread is not None and time.time() - t0 < 30:
        QApplication.processEvents()
        time.sleep(0.01)
    assert len(shown) == 1
    title, text = shown[0]
    assert title == "Le maillage n'a pas pu aboutir"
    assert text.startswith("Géométrie « cylinder » : radius doit être > 0") and "Journal" in text


def test_toml_typo_reported(win, monkeypatch):
    """Audit C4 : une clé mal orthographiée dans l'onglet TOML était ignorée en silence."""
    from PySide6.QtWidgets import QMessageBox

    from microrans.tomlio import dumps
    shown = []
    monkeypatch.setattr(QMessageBox, "information", lambda *a: shown.append(a[1:3]))
    monkeypatch.setattr(win, "quiet", False)
    win.new_case()
    cfg = dict(win.cfg)
    cfg["solver"] = {**cfg.get("solver", {}), "max_iters": 50}
    win.toml_edit.setPlainText(dumps(cfg))
    win._apply_toml()
    assert len(shown) == 1 and shown[0][0] == "À vérifier"
    assert "vouliez-vous dire « max_iter » ?" in shown[0][1]
    # cas de l'interface, sans faute : aucune fenêtre
    shown.clear()
    win.new_case()
    win.toml_edit.setPlainText(dumps(win.cfg))
    win._apply_toml()
    assert shown == []


def test_mesh_only_example_run_explained(win, monkeypatch):
    """Audit M4 : « Lancer » sur un exemple de maillage seul maillait ~1 min puis échouait
    (conditions aux limites manquantes) sans dire pourquoi."""
    from PySide6.QtWidgets import QMessageBox

    from microrans.cli import examples_dir
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[1:3]))
    monkeypatch.setattr(win, "quiet", False)
    win.open_case(examples_dir() / "mesh_naca_multi.toml")
    win.run_2d()
    assert win.thread is None                       # rien n'est lancé
    assert len(shown) == 1 and shown[0][0] == "Maillage à générer d'abord"
    assert "ne décrit qu'un maillage" in shown[0][1]


def test_toml_incomplete_case_still_loads(win, monkeypatch):
    """Un cas en cours d'écriture (viscosité pas encore donnée) se charge ; le manque est
    signalé (le contrôle complet a lieu au lancement)."""
    from PySide6.QtWidgets import QMessageBox

    from microrans.tomlio import dumps
    shown = []
    monkeypatch.setattr(QMessageBox, "information", lambda *a: shown.append(a[1:3]))
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[1:3]))
    monkeypatch.setattr(win, "quiet", False)
    win.new_case()
    cfg = dict(win.cfg)
    cfg["physics"] = {k: v for k, v in cfg["physics"].items() if k not in ("nu", "reynolds")}
    cfg["solver"] = {**cfg.get("solver", {}), "max_iter": 7}
    win.toml_edit.setPlainText(dumps(cfg))
    win._apply_toml()
    assert win.cfg["solver"]["max_iter"] == 7               # chargé
    assert len(shown) == 1 and shown[0][0] == "À vérifier"
    assert "À corriger avant de lancer : [physics] : donner la viscosité" in shown[0][1]
