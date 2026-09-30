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


def test_transition_example_roundtrip(win):
    """Cas T3A : modèle sst_gamma et convection de la turbulence conservés par les formulaires."""
    from microrans.cli import examples_dir
    from microrans.mesh2d.builder import load_config
    win.load_cfg(load_config(examples_dir() / "plaque_plane_transition_t3a.toml"))
    assert win.model_combo.currentData() == "sst_gamma"
    win._store_forms()
    assert win.cfg["physics"]["model"] == "sst_gamma"
    assert win.cfg["solver"]["convection_turb"] == "linearUpwindLimited"
