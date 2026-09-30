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
