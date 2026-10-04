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
    # audit U4 : seuls les paramètres de la loi sont affichés, sans « défaut » trompeur
    shown = {k for k, w in win.visc_fields.items() if win.visc_form.isRowVisible(w)}
    assert shown == {"nu0", "nu_inf", "lambda", "n", "nu_min", "nu_max"}
    assert win.visc_fields["lambda"].placeholderText() == "obligatoire"
    assert win.visc_fields["nu_max"].placeholderText() == "auto"
    assert win.visc_formula.text().startswith("ν = ν∞ + (ν₀ − ν∞)")
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
    assert not any(win.visc_form.isRowVisible(w) for w in win.visc_fields.values())
    assert not win.visc_form.isRowVisible(win.visc_formula)
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


def test_pages_fit_a_narrow_settings_panel(win):
    """Audit U2 : pages de 492 à 852 px de large minimum, tableau et formulaires coupés."""
    for i in range(win.pages.count()):
        page = win.pages.widget(i).widget()
        assert page.minimumSizeHint().width() <= 480, win.nav.item(i).text()


def test_bc_table_shows_only_useful_cells(win):
    """Audit U1 : dix colonnes toujours affichées, valeurs sans effet modifiables en silence
    (p d'une paroi, U d'une sortie…), colonne Type trop étroite pour son libellé."""
    import copy

    from PySide6.QtCore import Qt
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 2, "y0": 0, "y1": 1, "nx": 8, "ny": 4,
                    "names": {"left": "inlet", "right": "outlet", "bottom": "wall",
                              "top": "wall"}},
           "physics": {"nu": 0.01},
           "boundary": {"inlet": {"type": "inlet", "U": [1.0, 0.0]},
                        "outlet": {"type": "outlet"}, "wall": {"type": "wall"}}}
    win.load_cfg(copy.deepcopy(cfg))
    t = win.bc_table

    def shown():
        return [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())
                if not t.isColumnHidden(c)]
    rows = {t.item(r, 0).text(): r for r in range(t.rowCount())}

    def editable(name, col):
        return bool(t.item(rows[name], col).flags() & Qt.ItemIsEditable)
    assert shown() == ["Frontière", "Type", "Ux", "Uy", "p", "débit Q"]
    assert not editable("wall", 4) and editable("outlet", 4) and not editable("outlet", 2)
    # type changé : cellules réévaluées, colonne Type élargie au nouveau libellé
    cb = t.cellWidget(rows["inlet"], 1)
    cb.setCurrentIndex(cb.findData("pressure_inlet"))      # comme un clic (signal émis)
    assert editable("inlet", 4) and not editable("inlet", 2)
    assert win.cfg["boundary"]["inlet"] == {"type": "pressure_inlet", "p0": 0.0}
    assert "débit Q" not in shown()                 # plus aucune entrée en vitesse
    fm = t.fontMetrics()
    assert t.columnWidth(1) >= fm.horizontalAdvance(cb.currentText()) + 30
    # virgule décimale acceptée dans les cellules
    t.item(rows["outlet"], 4).setText("1,5")
    assert win.cfg["boundary"]["outlet"]["p"] == 1.5
    # thermique activée : colonnes T et q affichées
    win.load_cfg({**cfg, "energy": {"Pr": 0.7}})
    assert shown() == ["Frontière", "Type", "Ux", "Uy", "p", "T", "flux q", "débit Q"]


def test_derived_viscosity_shown(win):
    """Audit U3 : en mode « Nombre de Reynolds », le champ ν grisé affichait 0.01 (valeur
    par défaut) pour Re = 20 ; changer de mode doit garder la même viscosité."""
    from microrans.tomlio import loads as toml_loads

    from microrans.cli import examples_dir
    cfg = toml_loads((examples_dir() / "cylindre_re20.toml").read_text(encoding="utf-8"))
    ph = cfg["physics"]
    assert "reynolds" in ph and "nu" not in ph
    win.load_cfg(cfg)
    UL = ph.get("reference_velocity", 1.0) * ph.get("reference_length", 1.0)
    assert win.nu_edit.value() == pytest.approx(UL / ph["reynolds"], rel=1e-5)
    assert not win.nu_edit.isEnabled()
    win.re_edit.setText("40")                       # Re modifié : ν suit
    assert win.nu_edit.value() == pytest.approx(UL / 40, rel=1e-5)
    win.nu_mode.setCurrentIndex(win.nu_mode.findData("nu"))
    assert win.cfg["physics"]["nu"] == pytest.approx(UL / 40, rel=1e-5)
    assert "reynolds" not in win.cfg["physics"]
    assert win.re_edit.value() == pytest.approx(40) and not win.re_edit.isEnabled()


def test_non_newtonian_case_keeps_its_reference_viscosity(win):
    """Audit U11 : sans ν dans le cas (loi non newtonienne : ν de référence tiré de la loi),
    l'interface affichait 0.01 et l'écrivait dans le cas à la première modification."""
    from microrans.tomlio import loads as toml_loads

    from microrans.cli import examples_dir
    cfg = toml_loads((examples_dir() / "sang_carreau_artere.toml").read_text(
        encoding="utf-8"))
    assert "nu" not in cfg["physics"]
    win.load_cfg(cfg)
    assert win.nu_edit.text() == "" and win.nu_edit.placeholderText().startswith("auto")
    win.visc_fields["lambda"].setText("3.3")
    assert "nu" not in win.cfg["physics"]
    assert win.cfg["physics"]["viscosity"]["lambda"] == 3.3


def test_form_values_not_rounded(win):
    """Audit U12 : les champs affichaient 6 chiffres (ν = 1/550 → 0.00181818) et cette valeur
    arrondie était réécrite dans le cas à la première modification."""
    from microrans.gui.widgets import SciEdit
    for v in (1 / 550, 0.000632455532, 3000, 0.1, 1e-05, -9.81):
        e = SciEdit(v)
        assert e.value() == v
    assert SciEdit(3000.0).text() == "3000" and SciEdit(1e-05).text() == "1e-05"


def test_added_body_does_not_overlap(win):
    """Audit U5 : « Ajouter » posait un cercle de rayon 0.5 en (0, 0), exactement sur le
    cylindre ; en maillage en O, il créait un second corps refusé au maillage."""
    from microrans.tomlio import loads as toml_loads

    from microrans.cli import examples_dir
    from microrans.fv2d.validate import check_case
    from microrans.gui.widgets import set_combo
    cfg = toml_loads((examples_dir() / "cylindre_re20.toml").read_text(encoding="utf-8"))
    win.load_cfg(cfg)
    assert win.cfg["mesh"]["type"] == "ogrid" and not win.body_buttons[0].isEnabled()
    set_combo(win.mesh_type, "hybrid")
    win._mesh_type_changed()
    assert win.body_buttons[0].isEnabled()
    win._body_add()
    win._body_add()
    new = win.cfg["bodies"][1:]
    assert [b["name"] for b in new] == ["body2", "body3"]
    assert new[0]["center"][0] - new[0]["radius"] > 0.5      # à droite du cylindre
    win._store_forms()
    check_case(win.cfg)                                     # ni recouvrement ni collision
    # type changé : même place, même taille (avant : valeurs par défaut, sur le cylindre)
    win.body_list.setCurrentRow(1)
    for kind in ("rectangle", "ellipse", "naca", "circle"):
        win.b_type.setCurrentIndex(win.b_type.findData(kind))
        assert win.cfg["bodies"][1]["type"] == kind
        check_case(win.cfg)
    set_combo(win.mesh_type, "ogrid")
    win._mesh_type_changed()
    win._store_forms()
    assert not win.body_buttons[0].isEnabled()
    with pytest.raises(ValueError, match="le maillage en O entoure exactement un corps "
                                         r"\(3 donné\(s\)\)"):
        check_case(win.cfg)


def test_home_examples_grouped(win):
    """Audit U6 : accueil classé, description et bouton d'ouverture (le double-clic seul
    n'était pas découvrable)."""
    from PySide6.QtCore import Qt
    tree = win.examples
    first = tree.topLevelItem(0)
    assert first.text(0) == "Commencer ici" and not first.flags() & Qt.ItemIsSelectable
    assert tree.topLevelItem(tree.topLevelItemCount() - 1).text(0).startswith("Maillage seul")
    assert not win.example_open.isEnabled()
    tree.setCurrentItem(first.child(0))
    assert win.example_open.isEnabled()
    assert "cavite_re100.toml" in win.example_info.text()
    assert "microrans run2d" not in win.example_info.text()


def test_turbulence_settings_follow_model(win):
    """Audit U9 : lois de paroi proposées avec k-ε et transition (refus au lancement) ;
    réglages de turbulence actifs en laminaire."""
    from microrans.tomlio import loads as toml_loads

    from microrans.cli import examples_dir
    cfg = toml_loads((examples_dir() / "cavite_re100.toml").read_text(encoding="utf-8"))
    win.load_cfg(cfg)
    wf = win.wall_treat_combo.model().item(win.wall_treat_combo.findData("wall_function"))
    assert not win.wall_treat_combo.isEnabled() and not win.turb_fields[0].isEnabled()
    win.model_combo.setCurrentIndex(win.model_combo.findData("sst"))
    assert win.wall_treat_combo.isEnabled() and wf.isEnabled()
    assert win.cfg["physics"]["model"] == "sst"
    win.model_combo.setCurrentIndex(win.model_combo.findData("ke"))
    assert not wf.isEnabled() and "y⁺ ≈ 1" in wf.toolTip()


def test_sweep_section_loaded(win):
    """Audit U13 : [sweep] du cas ignorée par l'interface (exemple de polaire : « Lancer »
    ne calculait qu'un point, balayage avec les valeurs par défaut -4:12:2)."""
    from microrans.tomlio import loads as toml_loads

    from microrans.cli import examples_dir
    cfg = toml_loads((examples_dir() / "naca0012_polaire.toml").read_text(encoding="utf-8"))
    win.load_cfg(cfg)
    assert win.sweep_param.currentData() == "physics.angle_of_attack"
    assert win.sweep_values.text() == "-4:14:2" and win.sweep_cont.isChecked()
    assert "Lancer le balayage" in win.run_info.text()
    win.load_cfg(toml_loads((examples_dir() / "cavite_re100.toml").read_text(
        encoding="utf-8")))
    assert win.run_info.text() == ""


def test_mesh_only_example_after_meshing_explained(win, monkeypatch):
    """Exemple de maillage seul, maillage généré puis « Lancer » : la physique manque (avant
    U11, l'interface inventait ν = 0.01 et lançait un calcul) ; le refus dit quoi faire."""
    import copy

    from PySide6.QtWidgets import QMessageBox

    from microrans.cli import examples_dir
    shown = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a: shown.append(a[1:3]))
    monkeypatch.setattr(win, "quiet", False)
    win.open_case(examples_dir() / "mesh_cylindre_hybride.toml")
    cfg = copy.deepcopy(win.cfg)
    cfg["boundary"] = {"cylinder": {"type": "wall"}}            # frontières après maillage
    assert not win._case_ok(cfg)
    assert shown[-1][0] == "Réglages à corriger"
    assert "donner la viscosité nu" in shown[-1][1]
    assert "ne décrivait qu'un maillage" in shown[-1][1]


def test_mesh_quality_alert_shown(win):
    from test_report import _sheared_mesh
    win._mesh_done(_sheared_mesh())
    assert "ATTENTION : non-orthogonalité max 78.7° &gt; 70°" in win.mesh_info.text()


def test_help_menu_links_the_guides(win):
    # audit D1-D3 : guides accessibles depuis l'interface (fichiers présents dans docs/)
    from pathlib import Path
    menu = [m for m in win.menuBar().actions() if m.text() == "&Aide"][0].menu()
    texts = [a.text() for a in menu.actions()]
    assert texts[:4] == ["Premier calcul (tutoriel)",
                         "Dépannage : divergence, erreurs, résultats douteux",
                         "Glossaire des termes", "Toutes les clés du fichier de cas"]
    docs = Path(__file__).parents[1] / "docs"
    for name in ("tutoriel.md", "depannage.md", "glossaire.md", "reference_cas.md"):
        assert (docs / name).is_file(), name


def test_out_of_list_values_kept(win):
    # avant : convection_U = "linearUpwindLimited" (absent de la liste) devenait linearUpwind
    # en silence ; une valeur hors liste est désormais gardée (et vérifiée avant le calcul)
    import copy
    from microrans.tomlio import loads as toml_loads

    from microrans.cli import examples_dir
    cfg = toml_loads((examples_dir() / "cavite_re100.toml").read_text(encoding="utf-8"))
    cfg["solver"].update(convection_U="linearUpwindLimited", solver_p="mon_solveur")
    win.load_cfg(copy.deepcopy(cfg))
    win._store_forms()
    assert win.cfg["solver"]["convection_U"] == "linearUpwindLimited"
    assert win.cfg["solver"]["solver_p"] == "mon_solveur"
    combo = next(w for p, w, _, _ in win.binder.items if p == ("solver", "solver_p"))
    assert combo.currentText() == "mon_solveur (valeur du fichier)"
    n = combo.count()
    cfg["solver"]["solver_p"] = "auto"
    win.load_cfg(copy.deepcopy(cfg))
    win._store_forms()
    assert win.cfg["solver"]["solver_p"] == "auto" and combo.count() == n - 1


def _wait(win, timeout=120):
    import time

    from PySide6.QtWidgets import QApplication
    t0 = time.time()
    while win.thread is not None and time.time() - t0 < timeout:
        QApplication.processEvents()
        time.sleep(0.02)


def test_3d_box_case_form_roundtrip(win):
    """Lot E : pavé 3D dans l'interface — 3 composantes (raffinement, vitesses), colonne Uz
    des conditions limites, cas relu valide ; le raffinement absent vaut 1 en z (pas 0)."""
    import copy

    from microrans.cli import examples_dir
    from microrans.fv2d.validate import check_case
    win.open_case(examples_dir() / "conduite_carree_3d.toml")
    assert win._dim() == 3 and not win.bc_table.isColumnHidden(10)
    assert win.box_z.isVisibleTo(win) and not win.box_extrude.isVisibleTo(win)
    before = copy.deepcopy(win.cfg["mesh"])
    win._store_forms()
    m = win.cfg["mesh"]
    assert m["grading"] == [1.0, 1.0, 1.0]
    assert {k: m[k] for k in before} == before                   # rien de perdu
    assert check_case(copy.deepcopy(win.cfg)) == []


def test_extrude_checkbox_switches_case_to_3d(win):
    """Lot E : un cas 2D extrudé depuis l'interface (case à cocher) devient 3D : vitesses à 3
    composantes, bloc [mesh.extrude] écrit ; décoché, le cas redevient 2D."""
    import copy

    from microrans.cli import examples_dir
    from microrans.fv2d.validate import check_case
    win.open_case(examples_dir() / "cavite_re100.toml")
    assert win._dim() == 2 and win.bc_table.isColumnHidden(10)
    win.extrude_on.setChecked(True)
    win.ext_z0.set_value(0.0)
    win.ext_z1.set_value(0.1)
    win.ext_nz.setValue(2)
    win._store_forms()
    ext = win.cfg["mesh"]["extrude"]
    assert (ext["z0"], ext["z1"], ext["nz"]) == (0.0, 0.1, 2)
    assert ext["periodic"] == [["back", "front"]]
    lid = next(b for b in win.cfg["boundary"].values() if b.get("U"))
    assert len(lid["U"]) == 3 and not win.bc_table.isColumnHidden(10)
    assert check_case(copy.deepcopy(win.cfg)) == []
    win.extrude_on.setChecked(False)
    win._store_forms()
    assert "extrude" not in win.cfg["mesh"]
    lid = next(b for b in win.cfg["boundary"].values() if b.get("U"))
    assert len(lid["U"]) == 2 and win.bc_table.isColumnHidden(10)


def test_renamed_extrusion_faces_keep_periodicity(win):
    """C13 (audit 2) : canal_turbulent_3d nomme ses faces z « z0 » / « z1 » ; l'interface
    comparait aux noms back / front littéraux, effaçait la périodicité au premier
    enregistrement du formulaire et calculait un canal fermé (C_d du fond 0.00791 au lieu
    de 0.008889)."""
    import copy

    from microrans.cli import examples_dir
    from microrans.tomlio import loads
    path = examples_dir() / "canal_turbulent_3d.toml"
    ref = loads(path.read_text(encoding="utf-8"))["mesh"]["extrude"]
    win.open_case(path)
    assert win.ext_ends.currentData() == "periodic"
    win._store_forms()
    ext = win.cfg["mesh"]["extrude"]
    assert ext["periodic"] == ref["periodic"] == [["z0", "z1"]]
    assert "patch_types" not in ext or not set(ext["patch_types"]) & {"z0", "z1"}
    assert copy.deepcopy(ext) == ref


def test_run_after_mesh_settings_change_remeshes(win, monkeypatch, tmp_path):
    """U14 (audit 2) : maillage généré, puis réglages du maillage modifiés, puis « Lancer » :
    le calcul tournait sur l'ancien maillage sans le dire."""
    from microrans.cli import examples_dir
    monkeypatch.setenv("MICRORANS_RESULTS", str(tmp_path))
    win.open_case(examples_dir() / "cavite_re100.toml")
    win.cfg["mesh"].update(nx=6, ny=6)
    win.cfg["solver"].update(max_iter=2)
    win.load_cfg(win.cfg)
    win.quiet = True
    win.errors.clear()
    win.generate_mesh()
    _wait(win)
    assert win.mesh.n_cells == 36
    win.cfg["mesh"].update(nx=8, ny=8)               # comme un changement du formulaire
    win.load_cfg(win.cfg)
    win.log_view.clear()
    win.run_2d()
    _wait(win)
    assert not win.errors and win.summary["n_cells"] == 64
    assert "le calcul remaille" in win.log_view.toPlainText()
    win.log_view.clear()
    win.run_2d()                                     # réglages inchangés : même maillage
    _wait(win)
    assert win.summary["n_cells"] == 64
    assert "remaille" not in win.log_view.toPlainText()


def test_3d_case_meshed_run_and_plotted(win, monkeypatch, tmp_path):
    """Lot E : cas 3D de bout en bout dans l'interface (avant : refusé) — maillage, calcul,
    champs dans un plan z (Uz compris), plan hors du domaine expliqué au lieu d'une figure
    vide."""
    from microrans.cli import examples_dir
    monkeypatch.setenv("MICRORANS_RESULTS", str(tmp_path))
    win.open_case(examples_dir() / "cavite_cubique_re100_3d.toml")
    win.cfg["mesh"].update(nx=4, ny=4, nz=4)
    win.cfg["solver"].update(max_iter=3)
    win.cfg["output"]["directory"] = str(tmp_path / "cav3d")
    win.load_cfg(win.cfg)
    win.quiet = True
    win.errors.clear()
    win.generate_mesh()
    _wait(win)
    assert win.mesh is not None and win.mesh.dim == 3 and win.mesh.n_cells == 64
    win.run_2d()
    _wait(win)
    assert not win.errors and win.summary["iterations"] == 3
    keys = [win.field_combo.itemData(i) for i in range(win.field_combo.count())]
    assert "Uz" in keys and win.slice_value.isEnabled() and win.slice_axis.isEnabled()
    for i in range(len(keys)):
        win.field_combo.setCurrentIndex(i)
        win.plot_field()
    assert "plan z = 0.5" in win.canvas.fig.axes[0].get_title()
    win.line_start.set_value((0.5, 0.0, 0.5))
    win.line_end.set_value((0.5, 1.0, 0.5))
    win.plot_line()
    win.slice_value.set_value(5.0)
    win.plot_field()
    texts = [t.get_text() for t in win.canvas.fig.axes[0].texts]
    assert any("hors du domaine" in t for t in texts)
    win.draw_mesh()
    win.slice_value.set_value(None)
    # lot E4 : coupes x = cte et y = cte (vue y, z et x, z), vorticité normale au plan
    for axis, labels in (("x", ("y", "z")), ("y", ("x", "z"))):
        win.slice_axis.setCurrentIndex(win.slice_axis.findData(axis))
        win.slice_value.set_value(0.3)
        win.vec_check.setChecked(True)
        for key in ("U_mag", "vorticity"):
            win.field_combo.setCurrentIndex(win.field_combo.findData(key))
            win.plot_field()
            ax = win.canvas.fig.axes[0]
            assert f"plan {axis} = 0.3" in ax.get_title()
            assert (ax.get_xlabel(), ax.get_ylabel()) == labels
        assert f"ω_{axis}" in ax.get_title()
        win.draw_mesh()
        assert any(f"Coupe {axis} = 0.3" in a.get_title() for a in win.canvas.fig.axes)
    win.slice_value.set_value(-1.0)
    win.plot_field()
    texts = [t.get_text() for t in win.canvas.fig.axes[0].texts]
    assert any("plan y = -1 hors du domaine" in t for t in texts)
    win.vec_check.setChecked(False)
    win.slice_axis.setCurrentIndex(0)
    win.slice_value.set_value(None)
    assert not win.errors
