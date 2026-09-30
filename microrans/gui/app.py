"""Fenêtre principale de microrans (PySide6).

Parcours « à la Fluent / SimFlow » : Accueil → Maillage → Physique → Conditions limites →
Numérique → Calcul → Résultats, plus un module Canal 1D. Tout est stocké dans un cas au
format TOML (le même que la ligne de commande) : l'onglet « Fichier de cas » montre et
permet d'éditer ce TOML (fonctions avancées : maillage multi-blocs, raffinements...).
"""
from __future__ import annotations

import copy
import json
import os
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QThread, QUrl
from PySide6.QtGui import QAction, QDesktopServices, QFont, QKeySequence
from PySide6.QtWidgets import (QApplication, QCheckBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QListWidget, QListWidgetItem, QMainWindow,
                               QMessageBox, QPlainTextEdit, QProgressBar, QPushButton,
                               QScrollArea, QSplitter, QStackedWidget, QTableWidget,
                               QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from .. import __version__
from ..cli import examples_dir
from ..fv2d.solver import BC_TYPES, TIME_SCHEMES
from ..linalg import SOLVERS
from ..models import MODELS
from ..solver import TIME_SCHEMES_1D
from ..tomlio import dumps, loads
from .widgets import (Binder, PlotCanvas, SciEdit, Vec2, Worker, apply_plot_style, combo,
                      debounce, set_combo)

APP_NAME = "microrans"
PAGES = ["Accueil", "Canal 1D", "Maillage", "Physique", "Conditions limites", "Numérique",
         "Calcul", "Résultats"]
MODEL_LABELS = [("laminar", "Laminaire"), ("sa", "Spalart-Allmaras"), ("ke", "k-ε (Launder-Sharma)"),
                ("kw", "k-ω (Wilcox 2006)"), ("sst", "k-ω SST (Menter)"),
                ("sst_gamma", "k-ω SST + transition γ (Menter 2015)")]
BC_LABELS = {"wall": "Paroi", "inlet": "Entrée (vitesse)", "outlet": "Sortie (pression)",
             "symmetry": "Symétrie", "farfield": "Champ lointain",
             "axis": "Axe (axisymétrique)", "pressure_inlet": "Entrée (pression totale)"}
MESH_TYPES = [("rectangle", "Rectangle structuré"), ("ogrid", "Structuré en O autour d'un corps"),
              ("unstructured", "Triangles (non structuré)"),
              ("hybrid", "Hybride : couches de quadrilatères + triangles"),
              ("file", "Importer un fichier (.msh Gmsh, .su2)"),
              ("blocks", "Multi-blocs (édition dans l'onglet TOML)")]
BODY_TYPES = [("circle", "Cercle"), ("rectangle", "Rectangle"), ("ellipse", "Ellipse"),
              ("naca", "Profil NACA 4 chiffres"), ("file", "Contour importé (.dat .csv .svg .dxf)")]
FIELD_LABELS = {"T": "température T", "U_mag": "|U|", "Ux": "U_x", "Uy": "U_y", "p": "pression p", "vorticity":
                "vorticité ω_z", "nut_over_nu": "ν_t / ν", "k": "k", "omega": "ω", "eps": "ε",
                "nu_tilde": "ν̃", "wall_distance": "distance à la paroi",
                "viscosity": "viscosité ν (non newtonien)", "shear_rate": "taux de cisaillement γ̇"}
VISCOSITY_MODELS = [("newtonian", "Newtonien (ν constante)"), ("power_law", "Loi puissance"),
                    ("carreau", "Carreau"), ("cross", "Cross"),
                    ("herschel_bulkley", "Herschel-Bulkley (seuil)"), ("bingham", "Bingham (seuil)"),
                    ("casson", "Casson (sang…)")]
# paramètres (clé, libellé) de chaque loi ; nu_min / nu_max communs
VISCOSITY_PARAMS = {"power_law": ("K", "n"), "carreau": ("nu0", "nu_inf", "lambda", "n"),
                    "cross": ("nu0", "nu_inf", "m", "n"), "herschel_bulkley": ("tau_y", "K", "n"),
                    "bingham": ("tau_y", "K"), "casson": ("tau_y", "nu_inf")}
VISCOSITY_LABELS = {"K": "K (consistance, m²/sⁿ⁻¹…)", "n": "n (indice)", "nu0": "ν₀ (γ̇ → 0)",
                    "nu_inf": "ν∞ (γ̇ → ∞)", "lambda": "λ (temps, s)", "m": "m (temps, s)",
                    "tau_y": "τ_y / ρ (seuil)", "nu_max": "ν max (bouchon / borne)",
                    "nu_min": "ν min (borne)"}

DEFAULT_CASE = {
    "mesh": {"type": "rectangle", "x0": 0.0, "x1": 1.0, "y0": 0.0, "y1": 1.0, "nx": 48,
             "ny": 48, "names": {"left": "walls", "right": "walls", "bottom": "walls",
                                  "top": "lid"}},
    "physics": {"nu": 0.01, "model": "laminar", "reference_velocity": 1.0,
                "reference_length": 1.0},
    "boundary": {"lid": {"type": "wall", "U": [1.0, 0.0]}, "walls": {"type": "wall"}},
    "solver": {"mode": "steady", "max_iter": 2000, "tol": 1e-6},
    "output": {"directory": "results/nouveau_cas"},
}


def results_root() -> Path:
    return Path(os.environ.get("MICRORANS_RESULTS", Path.home() / "microrans_resultats"))


def _form(parent=None):
    box = QGroupBox(parent) if isinstance(parent, str) else QGroupBox()
    if isinstance(parent, str):
        box.setTitle(parent)
    lay = QFormLayout(box)
    lay.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
    return box, lay


def _note(text):
    lab = QLabel(text)
    lab.setWordWrap(True)
    lab.setStyleSheet("color:#52514e;")
    return lab


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} {__version__} — RANS/URANS 1D-2D")
        self.resize(1400, 860)
        apply_plot_style()
        self.cfg = copy.deepcopy(DEFAULT_CASE)
        self.case_path: Path | None = None
        self.mesh = None
        self.solver = None
        self.summary = None
        self.history: list = []
        self.thread = None
        self.worker = None
        self._syncing = False
        self.errors: list[str] = []
        self.quiet = False                      # pas de boîtes de dialogue (auto-test)
        self._build_ui()
        self._build_menu()
        self.load_cfg(self.cfg)
        self.nav.setCurrentRow(0)

    # ================================================================== interface
    def _build_ui(self):
        self.nav = QListWidget()
        self.nav.setMaximumWidth(190)
        f = QFont()
        f.setPointSize(f.pointSize() + 1)
        self.nav.setFont(f)
        for i, p in enumerate(PAGES):
            self.nav.addItem(QListWidgetItem(f"{i}. {p}" if i else p))
        self.pages = QStackedWidget()
        self.binder = Binder(self)
        self.binder.changed.connect(self._form_changed)
        builders = [self._page_home, self._page_1d, self._page_mesh, self._page_physics,
                    self._page_bc, self._page_numerics, self._page_run, self._page_results]
        for b in builders:
            w = b()
            sc = QScrollArea()
            sc.setWidgetResizable(True)
            sc.setWidget(w)
            self.pages.addWidget(sc)
        self.nav.currentRowChanged.connect(self.pages.setCurrentIndex)
        self.pages.setMinimumWidth(430)

        self.canvas = PlotCanvas()
        self.toml_edit = QPlainTextEdit()
        self.toml_edit.setFont(QFont("Monospace"))
        tb = QWidget()
        tl = QVBoxLayout(tb)
        tl.addWidget(_note("Fichier de cas complet (même format que « microrans run2d »). "
                           "Modifiez-le puis cliquez « Appliquer » : toutes les options "
                           "avancées sont accessibles ici (maillage multi-blocs, "
                           "raffinements, expressions de vitesse en x, y...)."))
        tl.addWidget(self.toml_edit, 1)
        row = QHBoxLayout()
        b = QPushButton("Appliquer le TOML")
        b.clicked.connect(self._apply_toml)
        row.addWidget(b)
        row.addStretch(1)
        tl.addLayout(row)
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setFont(QFont("Monospace"))
        self.log_view.setMaximumBlockCount(20000)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.canvas, "Visualisation")
        self.tabs.addTab(tb, "Fichier de cas (TOML)")
        self.tabs.addTab(self.log_view, "Journal")

        split = QSplitter(Qt.Horizontal)
        split.addWidget(self.nav)
        split.addWidget(self.pages)
        split.addWidget(self.tabs)
        split.setStretchFactor(2, 1)
        split.setSizes([170, 460, 800])
        self.setCentralWidget(split)
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(260)
        self.progress.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress)
        self._toml_timer = debounce(self, 250, self._refresh_toml)
        self.canvas.message("Bienvenue !\nOuvrez un exemple ou créez un cas.")

    def _build_menu(self):
        m = self.menuBar().addMenu("&Fichier")
        for text, key, fn in [("Nouveau cas 2D", QKeySequence.New, self.new_case),
                              ("Ouvrir un cas…", QKeySequence.Open, self.open_case),
                              ("Enregistrer", QKeySequence.Save, self.save_case),
                              ("Enregistrer sous…", QKeySequence.SaveAs, self.save_case_as)]:
            a = QAction(text, self)
            a.setShortcut(key)
            a.triggered.connect(fn)
            m.addAction(a)
        ex = m.addMenu("Exemples")
        for f in sorted(examples_dir().glob("*.toml")):
            a = QAction(f.stem, self)
            a.triggered.connect(lambda _=False, p=f: self.open_case(p))
            ex.addAction(a)
        m.addSeparator()
        a = QAction("Quitter", self)
        a.setShortcut(QKeySequence.Quit)
        a.triggered.connect(self.close)
        m.addAction(a)
        h = self.menuBar().addMenu("&Aide")
        a = QAction("Documentation (GitHub)", self)
        a.triggered.connect(lambda: QDesktopServices.openUrl(
            QUrl("https://github.com/Kayoouu/Microrans#readme")))
        h.addAction(a)
        a = QAction("À propos", self)
        a.triggered.connect(lambda: QMessageBox.about(
            self, "À propos", f"<b>{APP_NAME} {__version__}</b><br>Solveur RANS/URANS 1D et 2D "
            "(volumes finis) avec mailleur intégré.<br>Modèles : Spalart-Allmaras, k-ε, k-ω, "
            "k-ω SST.<br>Code libre, résultats vérifiés : voir le README."))
        h.addAction(a)

    # ------------------------------------------------------------------ pages
    def _page_home(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        title = QLabel(f"<h2>{APP_NAME}</h2><p>Simulation d'écoulements turbulents "
                       "RANS / URANS en 1D (canal) et 2D (volumes finis, maillages "
                       "structurés, non structurés ou hybrides).</p>")
        title.setWordWrap(True)
        lay.addWidget(title)
        lay.addWidget(_note("Pour débuter : double-cliquez sur un exemple, puis suivez les "
                            "étapes 2 à 7 dans la colonne de gauche. « 6. Calcul » lance la "
                            "simulation ; « 7. Résultats » affiche les champs et les efforts."))
        self.examples = QListWidget()
        for f in sorted(examples_dir().glob("*.toml")):
            first = f.read_text(encoding="utf-8").splitlines()[0].lstrip("# ").strip()
            it = QListWidgetItem(f"{f.stem}\n    {first}")
            it.setData(Qt.UserRole, str(f))
            self.examples.addItem(it)
        self.examples.itemDoubleClicked.connect(
            lambda it: self.open_case(Path(it.data(Qt.UserRole))))
        lay.addWidget(QLabel("<b>Exemples</b> (double-clic pour ouvrir)"))
        lay.addWidget(self.examples, 1)
        row = QHBoxLayout()
        for text, fn in [("Nouveau cas 2D", self.new_case), ("Ouvrir un cas…", self.open_case),
                         ("Canal 1D", lambda: self.nav.setCurrentRow(1))]:
            b = QPushButton(text)
            b.clicked.connect(fn)
            row.addWidget(b)
        lay.addLayout(row)
        return w

    def _page_1d(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.addWidget(_note("Canal plan turbulent pleinement développé (1D) : RANS stationnaire "
                            "ou URANS à gradient de pression pulsé. Très rapide : idéal pour "
                            "comparer les modèles de turbulence."))
        box, f = _form("Modèles")
        self.c1d_models = {}
        row = QHBoxLayout()
        for key, label in MODEL_LABELS[1:]:
            cb = self._check(label, key == "sa")
            self.c1d_models[key] = cb
            row.addWidget(cb)
        f.addRow(row)
        lay.addWidget(box)
        box, f = _form("Écoulement et maillage")
        self.c1d_retau = SciEdit(395.0)
        self.c1d_n = SciEdit(192)
        self.c1d_y1 = SciEdit(0.2)
        f.addRow("Re_τ = u_τ h / ν", self.c1d_retau)
        f.addRow("Nombre de mailles", self.c1d_n)
        f.addRow("y⁺ de la 1re maille", self.c1d_y1)
        self.c1d_mode = combo([("rans", "RANS stationnaire"), ("urans", "URANS pulsé")])
        f.addRow("Calcul", self.c1d_mode)
        lay.addWidget(box)
        box, f = _form("URANS : forçage f = 1 + A sin(ωt)")
        self.c1d_omega = SciEdit(0.01)
        self.c1d_amp = SciEdit(10.0)
        self.c1d_steps = SciEdit(64)
        self.c1d_scheme = combo([(k, f"{k} — {v}") for k, v in TIME_SCHEMES_1D.items()])
        set_combo(self.c1d_scheme, "sdirk2")
        f.addRow("ω⁺ (pulsation, unités de paroi)", self.c1d_omega)
        f.addRow("Amplitude A", self.c1d_amp)
        f.addRow("Pas par période", self.c1d_steps)
        f.addRow("Schéma en temps", self.c1d_scheme)
        lay.addWidget(box)
        b = QPushButton("Lancer le calcul 1D")
        b.clicked.connect(self.run_1d)
        lay.addWidget(b)
        self.c1d_status = _note("")
        lay.addWidget(self.c1d_status)
        lay.addStretch(1)
        return w

    def _page_mesh(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        B = self.binder
        box, f = _form("Type de maillage")
        self.mesh_type = B.combo(("mesh", "type"), MESH_TYPES, "rectangle")
        self.mesh_type.currentIndexChanged.connect(self._mesh_type_changed)
        f.addRow(self.mesh_type)
        lay.addWidget(box)
        # rectangle
        self.box_rect, f = _form("Rectangle")
        f.addRow("x min / x max", self._pair(B.sci(("mesh", "x0"), 0.0), B.sci(("mesh", "x1"), 1.0)))
        f.addRow("y min / y max", self._pair(B.sci(("mesh", "y0"), 0.0), B.sci(("mesh", "y1"), 1.0)))
        f.addRow("Cellules nx / ny", self._pair(B.int(("mesh", "nx"), 48, 1), B.int(("mesh", "ny"), 48, 1)))
        f.addRow("Progression (dernière/1re) x, y", B.vec(("mesh", "grading"), (1.0, 1.0)))
        for side, lab in [("left", "gauche"), ("right", "droite"), ("bottom", "bas"), ("top", "haut")]:
            f.addRow(f"Nom du bord {lab}", B.text(("mesh", "names", side), side))
        lay.addWidget(self.box_rect)
        # O-grid
        self.box_ogrid, f = _form("Maillage en O (un seul corps)")
        f.addRow("Cellules autour du corps", B.int(("mesh", "n_around"), 128, 8))
        f.addRow("Cellules radiales", B.int(("mesh", "n_radial"), 64, 4))
        f.addRow("Rayon du champ lointain", B.sci(("mesh", "farfield_radius"), 20.0))
        f.addRow("Épaisseur 1re maille", B.sci(("mesh", "first_height"), 1e-3))
        lay.addWidget(self.box_ogrid)
        # non structuré / hybride
        self.box_unst, f = _form("Tailles (triangles)")
        f.addRow("Taille max (loin)", B.sci(("mesh", "h_max"), 1.0))
        f.addRow("Taille sur les corps", B.sci(("mesh", "h_surface"), 0.05))
        f.addRow("Croissance avec la distance", B.sci(("mesh", "growth"), 0.2))
        lay.addWidget(self.box_unst)
        self.box_layers, f = _form("Couches limites (hybride)")
        f.addRow("Nombre de couches", B.int(("mesh", "layers", "n"), 10, 1))
        f.addRow("Épaisseur 1re couche", B.sci(("mesh", "layers", "first_height"), 1e-3))
        f.addRow("Raison géométrique", B.sci(("mesh", "layers", "ratio"), 1.2))
        lay.addWidget(self.box_layers)
        self.box_domain, f = _form("Domaine de calcul (rectangle)")
        f.addRow("x min / x max", self._pair(B.sci(("domain", "x0"), -10.0), B.sci(("domain", "x1"), 30.0)))
        f.addRow("y min / y max", self._pair(B.sci(("domain", "y0"), -10.0), B.sci(("domain", "y1"), 10.0)))
        for side, lab, d in [("left", "gauche", "inlet"), ("right", "droite", "outlet"),
                             ("bottom", "bas", "bottom"), ("top", "haut", "top")]:
            f.addRow(f"Nom du bord {lab}", B.text(("domain", "names", side), d))
        lay.addWidget(self.box_domain)
        self.box_file, f = _form("Fichier de maillage")
        self.mesh_path = B.text(("mesh", "path"), "")
        bb = QPushButton("Parcourir…")
        bb.clicked.connect(self._browse_mesh)
        f.addRow(self._pair(self.mesh_path, bb))
        lay.addWidget(self.box_file)
        # corps
        self.box_bodies = QGroupBox("Objets (corps solides)")
        bl = QVBoxLayout(self.box_bodies)
        self.body_list = QListWidget()
        self.body_list.setMaximumHeight(110)
        self.body_list.currentRowChanged.connect(self._body_selected)
        bl.addWidget(self.body_list)
        row = QHBoxLayout()
        for text, fn in [("Ajouter", self._body_add), ("Supprimer", self._body_remove)]:
            b = QPushButton(text)
            b.clicked.connect(fn)
            row.addWidget(b)
        bl.addLayout(row)
        bf = QFormLayout()
        self.b_type = combo(BODY_TYPES)
        self.b_name = QLineEdit_("body")
        self.b_center = Vec2((0.0, 0.0))
        self.b_radius = SciEdit(0.5)
        self.b_ab = Vec2((1.0, 0.5))
        self.b_corners = Vec2((-0.5, -0.5))
        self.b_corners2 = Vec2((0.5, 0.5))
        self.b_code = QLineEdit_("0012")
        self.b_chord = SciEdit(1.0)
        self.b_path = QLineEdit_("")
        self.b_inc = SciEdit(0.0)
        self.b_translate = Vec2((0.0, 0.0))
        self.b_rows = {}
        for key, label, wdg in [("type", "Type", self.b_type), ("name", "Nom du patch", self.b_name),
                                ("center", "Centre", self.b_center), ("radius", "Rayon", self.b_radius),
                                ("ab", "Demi-axes a, b", self.b_ab), ("c0", "Coin (x0, y0)", self.b_corners),
                                ("c1", "Coin (x1, y1)", self.b_corners2), ("code", "Code NACA", self.b_code),
                                ("chord", "Corde", self.b_chord), ("path", "Fichier du contour", self.b_path),
                                ("incidence", "Incidence (°)", self.b_inc),
                                ("translate", "Translation", self.b_translate)]:
            lab = QLabel(label)
            bf.addRow(lab, wdg)
            self.b_rows[key] = (lab, wdg)
        bl.addLayout(bf)
        for wdg in (self.b_center, self.b_radius, self.b_ab, self.b_corners, self.b_corners2,
                    self.b_chord, self.b_inc, self.b_translate):
            wdg.valueChanged.connect(self._body_changed)
        for wdg in (self.b_name, self.b_code, self.b_path):
            wdg.textChanged.connect(self._body_changed)
        self.b_type.currentIndexChanged.connect(self._body_changed)
        lay.addWidget(self.box_bodies)
        lay.addWidget(self.binder.check(
            ("mesh", "cut_axis"), False,
            "Couper à l'axe y = 0 (axisymétrique : ne garder que y > 0, frontière « axis »)"))
        row = QHBoxLayout()
        b = QPushButton("Générer le maillage")
        b.clicked.connect(self.generate_mesh)
        row.addWidget(b)
        b = QPushButton("Exporter…")
        b.clicked.connect(self.export_mesh)
        row.addWidget(b)
        lay.addLayout(row)
        self.mesh_info = _note("")
        lay.addWidget(self.mesh_info)
        lay.addStretch(1)
        return w

    def _page_physics(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        B = self.binder
        box, f = _form("Fluide (incompressible, grandeurs cinématiques ρ = 1)")
        self.nu_mode = combo([("nu", "Viscosité ν"), ("reynolds", "Nombre de Reynolds")])
        self.nu_mode.currentIndexChanged.connect(self._nu_mode_changed)
        f.addRow("Définir par", self.nu_mode)
        self.nu_edit = B.sci(("physics", "nu"), 0.01)
        self.re_edit = B.sci(("physics", "reynolds"), 100.0)
        f.addRow("ν (m²/s)", self.nu_edit)
        f.addRow("Re = U L / ν", self.re_edit)
        f.addRow("Vitesse de référence U", B.sci(("physics", "reference_velocity"), 1.0))
        f.addRow("Longueur de référence L", B.sci(("physics", "reference_length"), 1.0))
        f.addRow("Incidence α (°)", B.sci(("physics", "angle_of_attack"), None, True, "0"))
        f.addRow(B.check(("physics", "axisymmetric"), False,
                         "Axisymétrique : x = axe de révolution, y = rayon (tuyau, jet, sphère…)"))
        f.addRow(_note("Incidence : l'écoulement amont (entrées, champ lointain, vitesse "
                       "initiale) est tourné de α ; Cd et Cl sont donnés dans les axes de "
                       "l'écoulement."))
        lay.addWidget(box)
        box, f = _form("Viscosité (fluide non newtonien)")
        self.visc_combo = B.combo(("physics", "viscosity", "model"), VISCOSITY_MODELS,
                                  "newtonian")
        self.visc_combo.currentIndexChanged.connect(self._visc_changed)
        f.addRow("Loi", self.visc_combo)
        self.visc_fields = {}
        for key in ("K", "n", "nu0", "nu_inf", "lambda", "m", "tau_y", "nu_max", "nu_min"):
            wdg = B.sci(("physics", "viscosity", key), None, True, "défaut")
            self.visc_fields[key] = wdg
            f.addRow(VISCOSITY_LABELS[key], wdg)
        f.addRow(_note("ν = ν(γ̇), grandeurs cinématiques (divisées par ρ) ; laminaire "
                       "uniquement. Sans ν ci-dessus, ν de référence = ν(U/L). Écoulement "
                       "entraîné par une force (sans entrée) : relaxation U = 1 conseillée."))
        lay.addWidget(box)
        box, f = _form("Turbulence")
        self.model_combo = B.combo(("physics", "model"), MODEL_LABELS, "laminar")
        f.addRow("Modèle", self.model_combo)
        f.addRow("Intensité turbulente amont", B.sci(("turbulence", "intensity"), 0.001))
        f.addRow("Rapport ν_t/ν amont", B.sci(("turbulence", "viscosity_ratio"), 0.1))
        f.addRow("Traitement pariétal", B.combo(("solver", "wall_treatment"), [
            ("resolved", "Résolu jusqu'à la paroi (y⁺ ≈ 1)"),
            ("wall_function", "Lois de paroi, Spalding (y⁺ ≈ 30 à 300)")], "resolved"))
        f.addRow(_note("SA : ν̃ = 3ν en amont (recommandation NASA TMR). Lois de paroi : SA, "
                       "k-ω, SST (pas le k-ε bas-Reynolds) ; maillages 3 à 10× plus légers "
                       "près des parois, précision ~2-5 % (voir README)."))
        lay.addWidget(box)
        self.box_energy, f = _form("Thermique (équation de l'énergie, Boussinesq)")
        self.energy_on = self._check("Résoudre la température", False)
        self.energy_on.toggled.connect(self._energy_toggled)
        f.addRow(self.energy_on)
        self.energy_fields = [
            ("Prandtl Pr = ν/α", B.sci(("energy", "Pr"), 0.71)),
            ("Prandtl turbulent Pr_t", B.sci(("energy", "Pr_t"), 0.85)),
            ("β (dilatation ; 0 = sans flottabilité)", B.sci(("energy", "beta"), 0.0)),
            ("Gravité (gx, gy)", B.vec(("energy", "gravity"), (0.0, -9.81))),
            ("Température de référence", B.sci(("energy", "T_ref"), 0.0)),
            ("ΔT de référence (Nusselt)", B.sci(("energy", "delta_T"), 1.0))]
        for lab, wdg in self.energy_fields:
            f.addRow(lab, wdg)
        f.addRow(_note("Conditions de paroi : colonne T (température imposée) ou q (flux "
                       "entrant) de la page Conditions limites ; sinon adiabatique."))
        lay.addWidget(self.box_energy)
        box, f = _form("Scalaires transportés (concentration, polluant, âge du fluide…)")
        self.scalar_table = QTableWidget(0, 5)
        self.scalar_table.setHorizontalHeaderLabels(["Nom", "Diffusivité D (m²/s)", "Sc_t",
                                                     "Source S", "Valeur initiale"])
        self.scalar_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.scalar_table.setMinimumHeight(110)
        self.scalar_table.itemChanged.connect(lambda *_: self._scalars_changed())
        f.addRow(self.scalar_table)
        row = QHBoxLayout()
        for text, fn in (("Ajouter", self._scalar_add), ("Supprimer", self._scalar_remove)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            row.addWidget(b)
        row.addStretch(1)
        f.addRow(row)
        f.addRow(_note("∂c/∂t + ∇·(Uc) = ∇·((D + ν_t/Sc_t)∇c) + S. Valeurs aux frontières : "
                       "colonne « scalaires » de la page Conditions limites (ex. c=1) ; "
                       "défaut 0 en entrée, flux nul ailleurs. Source S = 1 : âge moyen du "
                       "fluide (temps de séjour). Bilan par frontière dans summary.json."))
        lay.addWidget(box)
        box, f = _form("Conditions initiales et forces")
        f.addRow("Vitesse initiale (Ux, Uy)", B.vec(("initial", "U"), (0.0, 0.0)))
        f.addRow("Perturbation du sillage", B.sci(("initial", "perturbation"), None, True, "0"))
        f.addRow("Force volumique (fx, fy)", B.vec(("physics", "body_force"), (0.0, 0.0)))
        self.restart_path = B.text(("initial", "restart"), "")
        self.restart_path.setPlaceholderText("vide : démarrer de l'état initial ci-dessus")
        bb = QPushButton("Parcourir…")
        bb.clicked.connect(self._browse_restart)
        f.addRow("Repartir d'un calcul (checkpoint.npz)", self._pair(self.restart_path, bb))
        f.addRow(_note("Même maillage : reprise exacte. Maillage différent : les champs sont "
                       "interpolés (ex. démarrer un maillage fin depuis un calcul grossier)."))
        lay.addWidget(box)
        lay.addStretch(1)
        return w

    def _page_bc(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.addWidget(_note("Une ligne par frontière (patch) du maillage. Paroi : adhérence "
                            "(U = paroi mobile) ; Entrée : U imposée, ou débit Q (vitesse "
                            "normale uniforme ; U ignorée) ; Entrée (pression totale) : colonne "
                            "p = p0 ; Sortie : p imposée ; Champ lointain : U∞ en entrée, p∞ en "
                            "sortie selon le signe de U∞·n. Les composantes de U acceptent des "
                            "expressions en x, y (ex. 6*y*(1-y))."))
        self.bc_table = QTableWidget(0, 9)
        self.bc_table.setHorizontalHeaderLabels(["Patch", "Type", "Ux", "Uy", "p", "T", "flux q",
                                                 "débit Q", "scalaires"])
        self.bc_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.bc_table.itemChanged.connect(lambda *_: self._bc_changed())
        lay.addWidget(self.bc_table, 1)
        return w

    def _page_numerics(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        B = self.binder
        box, f = _form("Type de calcul")
        self.mode_combo = B.combo(("solver", "mode"), [("steady", "Stationnaire (RANS)"),
                                                       ("transient", "Instationnaire (URANS)")],
                                  "steady")
        self.mode_combo.currentIndexChanged.connect(self._mode_changed)
        f.addRow(self.mode_combo)
        lay.addWidget(box)
        self.box_steady, f = _form("Stationnaire : SIMPLE / SIMPLEC")
        f.addRow("Algorithme", B.combo(("solver", "algorithm"), [("SIMPLEC", "SIMPLEC (recommandé)"),
                                                                   ("SIMPLE", "SIMPLE")], "SIMPLEC"))
        f.addRow("Sous-relaxation U (vide = auto)", B.sci(("solver", "relax_U"), None, True, "auto"))
        f.addRow("Sous-relaxation turbulence", B.sci(("solver", "relax_turb"), 0.8))
        f.addRow("Itérations max", B.int(("solver", "max_iter"), 3000, 1))
        f.addRow("Tolérance des résidus", B.sci(("solver", "tol"), 1e-5))
        f.addRow("Arrêt sur efforts stabilisés (vide = non)",
                 B.sci(("solver", "monitor_tol"), None, True, "ex. 1e-5"))
        f.addRow("Pseudo-transitoire local, Courant (vide = non)",
                 B.sci(("solver", "pseudo_cfl"), None, True, "ex. 200 (canaux très fins)"))
        f.addRow("Démarrage multigrille (niveaux grossiers)",
                 B.int(("solver", "fmg_levels"), 0, 0, 4))
        f.addRow(_note("Multigrille : le cas est d'abord résolu sur des maillages 2, 4… fois "
                       "plus grossiers. Mesuré : cavité 128² 2.5× plus rapide, cylindre 1.4× ; "
                       "sans gain sur la plaque plane turbulente. Maillages générés seulement."))
        lay.addWidget(self.box_steady)
        self.box_transient, f = _form("Instationnaire")
        opts = [("auto", "auto — choix automatique (recommandé)")] + [
            (k, f"{k} — {v['label']}") for k, v in TIME_SCHEMES.items()]
        f.addRow("Schéma en temps", B.combo(("solver", "time_scheme"), opts, "auto"))
        f.addRow("Pas de temps Δt", B.sci(("solver", "dt"), 0.01))
        f.addRow("Temps final", B.sci(("solver", "t_end"), 10.0))
        f.addRow(B.check(("solver", "adjust_dt"), False, "Pas de temps adaptatif (Courant)"))
        f.addRow("Courant visé", B.sci(("solver", "max_co"), 1.0))
        f.addRow("PIMPLE : boucles externes / corrections p",
                 self._pair(B.int(("solver", "n_outer"), 2, 1), B.int(("solver", "n_corr"), 2, 1)))
        f.addRow(_note("Implicite (euler, backward, crankNicolson) : robuste, grands Δt. "
                       "Explicite (rk2-rk4, ab2) : moins cher par pas mais Δt limité par le "
                       "Courant et la diffusion près des parois. Voir docs/schemas_temps.png."))
        lay.addWidget(self.box_transient)
        box, f = _form("Schémas spatiaux et solveurs linéaires")
        f.addRow("Convection U", B.combo(("solver", "convection_U"),
                                         [("linearUpwind", "linearUpwind (ordre 2)"),
                                          ("upwind", "upwind (ordre 1, robuste)")], "linearUpwind"))
        f.addRow("Convection turbulence", B.combo(
            ("solver", "convection_turb"),
            [("upwind", "upwind (ordre 1, robuste)"),
             ("linearUpwindLimited", "linearUpwindLimited (ordre 2 limité, conseillé en transition)"),
             ("linearUpwind", "linearUpwind (ordre 2)")], "upwind"))
        sol = [(s, s) for s in SOLVERS]
        f.addRow("Solveur pression", B.combo(("solver", "solver_p"), sol, "auto"))
        f.addRow("Solveur vitesse / turbulence", B.combo(("solver", "solver_U"), sol, "auto"))
        f.addRow("Matériel de calcul", B.combo(("solver", "backend"),
                                                [("cpu", "CPU (NumPy / SciPy)"),
                                                 ("gpu", "GPU NVIDIA (CuPy, expérimental)")],
                                                "cpu"))
        f.addRow(_note("GPU : nécessite une carte NVIDIA, CUDA et CuPy (pip install "
                       "cupy-cuda12x) ; utile seulement au-delà de ~10⁵ cellules. Chemin "
                       "testé par émulation, pas sur une vraie carte : voir le README."))
        lay.addWidget(box)
        box, f = _form("Sorties")
        f.addRow("Dossier de résultats", B.text(("output", "directory"), "results/cas"))
        f.addRow(_note(f"Chemin relatif : placé dans {results_root()}"))
        self.probes_edit = B.text(("output", "probes"), "")
        self.probes_edit.setPlaceholderText("x y ; x y …  (vide : aucune)")
        f.addRow("Sondes (Ux, Uy, p à chaque itération)", self.probes_edit)
        f.addRow("Moyennes temporelles à partir de t =",
                 B.sci(("output", "average_from"), None, True, "vide : non (instationnaire)"))
        f.addRow("Animation GIF (instationnaire)", B.combo(
            ("output", "animate"), [(None, "aucune"), ("vorticity", "vorticité"),
                                    ("U_mag", "|U|"), ("p", "pression"),
                                    ("T", "température"), ("viscosity", "viscosité")], None))
        f.addRow(_note("Animation d'un scalaire : animate = \"nom\" dans l'onglet TOML."))
        lay.addWidget(box)
        lay.addStretch(1)
        return w

    def _page_run(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.addWidget(_note("Le calcul tourne en arrière-plan : l'interface reste utilisable. "
                            "Les résidus (stationnaire) ou les efforts (instationnaire) "
                            "s'affichent en direct."))
        row = QHBoxLayout()
        self.run_btn = QPushButton("Lancer le calcul")
        self.run_btn.clicked.connect(self.run_2d)
        self.stop_btn = QPushButton("Arrêter")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop)
        self.continue_btn = QPushButton("Continuer le calcul précédent")
        self.continue_btn.setToolTip(
            "Repart de la sauvegarde checkpoint.npz du dossier de résultats. Stationnaire : "
            "« Itérations max » itérations de plus. Instationnaire : jusqu'au nouveau temps final.")
        self.continue_btn.clicked.connect(self.continue_2d)
        row.addWidget(self.run_btn)
        row.addWidget(self.continue_btn)
        row.addWidget(self.stop_btn)
        lay.addLayout(row)
        lay.addWidget(_note("Le calcul est sauvegardé (checkpoint.npz) toutes les 5 minutes, "
                            "à la fin et à l'arrêt : un calcul arrêté ou interrompu peut être "
                            "poursuivi avec « Continuer »."))
        box, f = _form("Polaire / balayage d'un paramètre")
        self.sweep_param = combo([("physics.angle_of_attack", "Incidence α (°) : polaire"),
                                  ("physics.reynolds", "Nombre de Reynolds"),
                                  ("physics.nu", "Viscosité ν")])
        self.sweep_values = QLineEdit_("-4:12:2")
        self.sweep_values.setToolTip("début:fin:pas (fin incluse) ou liste a, b, c")
        self.sweep_cont = QCheckBox("Continuation : chaque point part du précédent")
        self.sweep_cont.setChecked(True)
        self.sweep_btn = QPushButton("Lancer le balayage")
        self.sweep_btn.clicked.connect(self.run_sweep_gui)
        f.addRow("Paramètre", self.sweep_param)
        f.addRow("Valeurs", self.sweep_values)
        f.addRow(self.sweep_cont)
        f.addRow(self.sweep_btn)
        f.addRow(_note("Un calcul complet par valeur, sur le même maillage. Résultats : "
                       "balayage.csv, polaire.png et un sous-dossier par point. Points non "
                       "convergés (ex. après le décrochage) : marqueurs creux."))
        lay.addWidget(box)
        self.run_info = _note("")
        lay.addWidget(self.run_info)
        lay.addStretch(1)
        return w

    def _page_results(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        box, f = _form("Champ")
        self.field_combo = combo([("U_mag", "|U|")])
        self.cmap_combo = combo([(c, c) for c in ("viridis", "RdBu_r", "magma", "coolwarm",
                                                  "turbo", "cividis")])
        self.zoom_check = self._check("Zoom sur les corps", True)
        self.mesh_check = self._check("Superposer le maillage", False)
        self.vec_check = self._check("Vecteurs vitesse", False)
        f.addRow("Grandeur", self.field_combo)
        f.addRow("Palette", self.cmap_combo)
        f.addRow(self.zoom_check)
        f.addRow(self.mesh_check)
        f.addRow(self.vec_check)
        b = QPushButton("Tracer le champ")
        b.clicked.connect(self.plot_field)
        f.addRow(b)
        lay.addWidget(box)
        box, f = _form("Distributions pariétales")
        self.wall_combo = combo([])
        self.wall_q = combo([("Cp", "Coefficient de pression Cp"), ("Cf", "Frottement Cf"),
                             ("yplus", "y⁺"), ("q", "Flux de chaleur q"),
                             ("Tw", "Température de paroi")])
        f.addRow("Paroi", self.wall_combo)
        f.addRow("Grandeur", self.wall_q)
        b = QPushButton("Tracer")
        b.clicked.connect(self.plot_wall)
        f.addRow(b)
        lay.addWidget(box)
        box, f = _form("Profil le long d'une ligne (grandeur choisie ci-dessus)")
        self.line_start, self.line_end = Vec2((0.0, 0.0)), Vec2((1.0, 0.0))
        f.addRow("Début (x, y)", self.line_start)
        f.addRow("Fin (x, y)", self.line_end)
        b = QPushButton("Tracer le profil")
        b.clicked.connect(self.plot_line)
        f.addRow(b)
        lay.addWidget(box)
        b = QPushButton("Historique (résidus / efforts)")
        b.clicked.connect(lambda: self._plot_history(final=True))
        lay.addWidget(b)
        self.summary_view = QPlainTextEdit()
        self.summary_view.setReadOnly(True)
        self.summary_view.setMinimumHeight(180)
        lay.addWidget(QLabel("<b>Résumé</b>"))
        lay.addWidget(self.summary_view)
        b = QPushButton("Ouvrir l'animation")
        b.clicked.connect(self.open_animation)
        lay.addWidget(b)
        b = QPushButton("Ouvrir le dossier de résultats")
        b.clicked.connect(lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.out_dir()))))
        lay.addWidget(b)
        lay.addStretch(1)
        return w

    # ------------------------------------------------------------------ petits outils
    def _check(self, text, value):
        from PySide6.QtWidgets import QCheckBox
        c = QCheckBox(text)
        c.setChecked(value)
        return c

    @staticmethod
    def _pair(a, b):
        w = QWidget()
        lay = QHBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(a)
        lay.addWidget(b)
        return w

    def log(self, text):
        self.log_view.appendPlainText(text)

    def out_dir(self) -> Path:
        d = Path(self.cfg.get("output", {}).get("directory") or "results/cas")
        return d if d.is_absolute() else results_root() / d

    def base_dir(self) -> Path:
        return self.case_path.parent if self.case_path else Path.cwd()

    # ================================================================== cas
    def load_cfg(self, cfg):
        self._syncing = True
        self.cfg = cfg
        ph = cfg.setdefault("physics", {})
        set_combo(self.nu_mode, "reynolds" if "reynolds" in ph and "nu" not in ph else "nu")
        self.binder.load(cfg)
        self._load_bodies()
        self._mesh_type_changed()
        self._nu_mode_changed()
        self._mode_changed()
        self.energy_on.setChecked("energy" in cfg)
        for _, wdg in self.energy_fields:
            wdg.setEnabled("energy" in cfg)
        self._visc_changed()
        self._fill_scalar_table()
        self._fill_bc_table()
        self._syncing = False
        self._refresh_toml()

    def _form_changed(self):
        if self._syncing:
            return
        self._store_forms()
        self._toml_timer.start()

    def _store_forms(self):
        self.binder.store(self.cfg)
        if not self.energy_on.isChecked():
            self.cfg.pop("energy", None)
        ph = self.cfg.setdefault("physics", {})
        v = ph.get("viscosity") or {}
        model = v.get("model", "newtonian")
        if model == "newtonian":
            ph.pop("viscosity", None)
        else:
            keep = {"model", "nu_min", "nu_max", "relax", "a", "gamma_ref",
                    *VISCOSITY_PARAMS.get(model, ())}
            ph["viscosity"] = {k: val for k, val in v.items() if k in keep}
        if self.nu_mode.currentData() == "nu":
            ph.pop("reynolds", None)
        else:
            ph.pop("nu", None)
        self._prune_mesh_keys()

    def _prune_mesh_keys(self):
        m = self.cfg.setdefault("mesh", {})
        kind = m.get("type", "rectangle")
        keep = {"rectangle": {"x0", "x1", "y0", "y1", "nx", "ny", "grading", "names", "periodic",
                              "patch_types"},
                "ogrid": {"n_around", "n_radial", "farfield_radius", "first_height", "center"},
                "unstructured": {"h_max", "h_surface", "growth", "refinements", "max_iter"},
                "hybrid": {"h_max", "h_surface", "growth", "refinements", "max_iter", "layers"},
                "file": {"path", "patch_types"}}.get(kind)
        if keep is None:
            return
        managed = {"x0", "x1", "y0", "y1", "nx", "ny", "grading", "names", "n_around",
                   "n_radial", "farfield_radius", "first_height", "h_max", "h_surface", "growth",
                   "layers", "path"}
        for k in list(m):
            if k in managed and k not in keep:
                m.pop(k)
        if kind not in ("unstructured", "hybrid"):
            self.cfg.pop("domain", None)
        if kind in ("rectangle", "file"):
            self.cfg.pop("bodies", None)

    def _refresh_toml(self):
        self.toml_edit.setPlainText(dumps(self.cfg))

    def _apply_toml(self):
        try:
            cfg = loads(self.toml_edit.toPlainText())
        except Exception as exc:                     # noqa: BLE001
            QMessageBox.warning(self, "TOML invalide", str(exc))
            return
        self.load_cfg(cfg)
        self.statusBar().showMessage("Cas mis à jour depuis le TOML", 4000)

    def new_case(self):
        self.case_path = None
        self.mesh = None
        self.load_cfg(copy.deepcopy(DEFAULT_CASE))
        self.nav.setCurrentRow(2)

    def open_case(self, path=None):
        if not path:
            path, _ = QFileDialog.getOpenFileName(self, "Ouvrir un cas", str(examples_dir()),
                                                  "Cas (*.toml *.json)")
            if not path:
                return
        path = Path(path)
        try:
            text = path.read_text(encoding="utf-8")
            cfg = json.loads(text) if path.suffix == ".json" else loads(text)
        except Exception as exc:                     # noqa: BLE001
            QMessageBox.warning(self, "Lecture impossible", str(exc))
            return
        self.case_path = path
        self.mesh = None
        self.solver = None
        self.load_cfg(cfg)
        self.setWindowTitle(f"{APP_NAME} — {path.name}")
        self.log(f"Cas ouvert : {path}")
        self.nav.setCurrentRow(2)
        if cfg.get("mesh", {}).get("type") not in ("blocks",):
            self.canvas.message(f"{path.stem}\nCliquez « Générer le maillage ».")

    def save_case(self):
        if self.case_path is None or examples_dir() in self.case_path.parents:
            return self.save_case_as()
        self._store_forms()
        self.case_path.write_text(dumps(self.cfg), encoding="utf-8")
        self.statusBar().showMessage(f"Enregistré : {self.case_path}", 4000)

    def save_case_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "Enregistrer le cas",
                                              str(results_root() / "cas.toml"), "Cas (*.toml)")
        if path:
            self.case_path = Path(path)
            self.save_case()

    # ------------------------------------------------------------------ maillage
    def _mesh_type_changed(self):
        kind = self.mesh_type.currentData()
        self.box_rect.setVisible(kind == "rectangle")
        self.box_ogrid.setVisible(kind == "ogrid")
        self.box_unst.setVisible(kind in ("unstructured", "hybrid"))
        self.box_layers.setVisible(kind == "hybrid")
        self.box_domain.setVisible(kind in ("unstructured", "hybrid"))
        self.box_file.setVisible(kind == "file")
        self.box_bodies.setVisible(kind in ("ogrid", "unstructured", "hybrid"))
        for box, on in [(self.box_rect, kind == "rectangle"), (self.box_ogrid, kind == "ogrid"),
                        (self.box_unst, kind in ("unstructured", "hybrid")),
                        (self.box_layers, kind == "hybrid"),
                        (self.box_domain, kind in ("unstructured", "hybrid")),
                        (self.box_file, kind == "file")]:
            box.setEnabled(on)
        if kind in ("ogrid", "unstructured", "hybrid") and not self.cfg.get("bodies"):
            self.cfg["bodies"] = [{"type": "circle", "center": [0.0, 0.0], "radius": 0.5,
                                   "name": "cylinder"}]
            self._load_bodies()

    def _browse_mesh(self):
        path, _ = QFileDialog.getOpenFileName(self, "Maillage", "", "Maillages (*.msh *.su2)")
        if path:
            self.mesh_path.setText(path)

    def _load_bodies(self):
        self.body_list.blockSignals(True)
        self.body_list.clear()
        for b in self.cfg.get("bodies", []):
            self.body_list.addItem(f"{b.get('name', b['type'])} ({b['type']})")
        self.body_list.blockSignals(False)
        if self.cfg.get("bodies"):
            self.body_list.setCurrentRow(0)
            self._body_selected(0)

    def _body_selected(self, row):
        bodies = self.cfg.get("bodies", [])
        if not 0 <= row < len(bodies):
            return
        b = bodies[row]
        self._body_sync = True
        set_combo(self.b_type, b.get("type", "circle"))
        self.b_name.setText(b.get("name", ""))
        self.b_center.set_value(b.get("center", (0.0, 0.0)))
        self.b_radius.set_value(b.get("radius", 0.5))
        self.b_ab.set_value((b.get("a", 1.0), b.get("b", 0.5)))
        self.b_corners.set_value((b.get("x0", -0.5), b.get("y0", -0.5)))
        self.b_corners2.set_value((b.get("x1", 0.5), b.get("y1", 0.5)))
        self.b_code.setText(str(b.get("code", "0012")))
        self.b_chord.set_value(b.get("chord", 1.0))
        self.b_path.setText(b.get("path", ""))
        self.b_inc.set_value(b.get("incidence", 0.0))
        self.b_translate.set_value(b.get("translate", (0.0, 0.0)))
        self._body_sync = False
        self._body_fields()

    def _body_fields(self):
        t = self.b_type.currentData()
        show = {"circle": {"center", "radius"}, "rectangle": {"c0", "c1"},
                "ellipse": {"center", "ab"}, "naca": {"code", "chord"}, "file": {"path"}}[t]
        for key, (lab, wdg) in self.b_rows.items():
            vis = key in ("type", "name", "incidence", "translate") or key in show
            lab.setVisible(vis)
            wdg.setVisible(vis)

    def _body_changed(self):
        self._body_fields()
        if getattr(self, "_body_sync", False):
            return
        row = self.body_list.currentRow()
        bodies = self.cfg.setdefault("bodies", [])
        if not 0 <= row < len(bodies):
            return
        t = self.b_type.currentData()
        b = {"type": t, "name": self.b_name.text().strip() or t}
        if t == "circle":
            b.update(center=self.b_center.value(), radius=self.b_radius.value() or 0.5)
        elif t == "rectangle":
            (x0, y0), (x1, y1) = self.b_corners.value(), self.b_corners2.value()
            b.update(x0=x0, y0=y0, x1=x1, y1=y1)
        elif t == "ellipse":
            a, bb = self.b_ab.value()
            b.update(center=self.b_center.value(), a=a, b=bb)
        elif t == "naca":
            b.update(code=self.b_code.text().strip() or "0012", chord=self.b_chord.value() or 1.0)
        else:
            b.update(path=self.b_path.text().strip())
        if self.b_inc.value():
            b["incidence"] = self.b_inc.value()
        if any(self.b_translate.value()):
            b["translate"] = self.b_translate.value()
        bodies[row] = b
        self.body_list.item(row).setText(f"{b['name']} ({t})")
        self._toml_timer.start()

    def _body_add(self):
        self.cfg.setdefault("bodies", []).append(
            {"type": "circle", "center": [0.0, 0.0], "radius": 0.5,
             "name": f"body{len(self.cfg['bodies']) + 1}"})
        self._load_bodies()
        self.body_list.setCurrentRow(len(self.cfg["bodies"]) - 1)
        self._toml_timer.start()

    def _body_remove(self):
        row = self.body_list.currentRow()
        if 0 <= row < len(self.cfg.get("bodies", [])):
            self.cfg["bodies"].pop(row)
            self._load_bodies()
            self._toml_timer.start()

    def generate_mesh(self):
        self._store_forms()
        cfg = copy.deepcopy(self.cfg)
        base = self.base_dir()
        if cfg.get("mesh", {}).get("type") == "file":
            p = Path(cfg["mesh"].get("path", ""))
            if not p.is_absolute() and self.case_path is None:
                cfg["mesh"]["path"] = str(p.resolve())

        def job(worker):
            from ..mesh2d.builder import build_mesh
            print("Génération du maillage…")
            return build_mesh(cfg, base_dir=base, verbose=True)
        self._start(job, self._mesh_done, "Maillage en cours…")

    def _mesh_done(self, mesh):
        self.mesh = mesh
        self.solver = None
        q = mesh.quality()
        types = ", ".join(f"{v} {k}" for k, v in q["cell_types"].items())
        self.mesh_info.setText(
            f"<b>{q['n_cells']} cellules</b> ({types})<br>non-orthogonalité max "
            f"{q['non_orthogonality_max_deg']:.1f}°, moyenne {q['non_orthogonality_mean_deg']:.1f}°"
            f"<br>asymétrie max {q['skewness_max']:.2f}, allongement max "
            f"{q['aspect_ratio_max']:.0f}<br>patches : "
            + ", ".join(f"{n} ({t})" for n, t, _ in mesh.all_boundary_patches()))
        self._sync_bc_with_mesh()
        self.draw_mesh()
        self.log(f"Maillage : {q['n_cells']} cellules.")

    def draw_mesh(self):
        from ..mesh2d.plot import plot_mesh
        if self.mesh is None:
            return
        ax = self.canvas.axes()
        plot_mesh(self.mesh, ax=ax, zoom=self._zoom() if self.zoom_check.isChecked() else None,
                  linewidth=0.25 if self.mesh.n_cells < 30000 else 0.1)
        self.canvas.draw()
        self.tabs.setCurrentIndex(0)

    def export_mesh(self):
        if self.mesh is None:
            QMessageBox.information(self, "Exporter", "Générez d'abord le maillage.")
            return
        path, filt = QFileDialog.getSaveFileName(
            self, "Exporter le maillage", str(results_root() / "maillage.msh"),
            "Gmsh (*.msh);;SU2 (*.su2);;VTK (*.vtk);;OpenFOAM (dossier polyMesh) (*)")
        if not path:
            return
        from ..mesh2d.io import write_mesh, write_openfoam
        try:
            if "OpenFOAM" in filt:
                write_openfoam(self.mesh, Path(path))
            else:
                write_mesh(self.mesh, Path(path))
            self.statusBar().showMessage(f"Maillage exporté : {path}", 5000)
        except Exception as exc:                     # noqa: BLE001
            QMessageBox.warning(self, "Export impossible", str(exc))

    def _zoom(self):
        from ..fv2d.post import _zoom
        if self.solver is not None:
            steady = self.cfg.get("solver", {}).get("mode", "steady") == "steady"
            return _zoom(self.solver, 4.0 if steady else 14.0)
        if self.mesh is None:
            return None
        walls = [n for n, t, _ in self.mesh.all_boundary_patches()
                 if self.cfg.get("boundary", {}).get(n, {}).get("type") == "wall"
                 or t == "wall"]
        pts = [self.mesh.points[e].reshape(-1, 2) for n, t, e in self.mesh.all_boundary_patches()
               if n in walls]
        if not pts:
            return None
        pts = np.vstack(pts)
        L = float(max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1])))
        bb = self.mesh.bbox()
        if L > 0.5 * max(bb[2] - bb[0], bb[3] - bb[1]):
            return None
        (x0, y0), (x1, y1) = pts.min(axis=0), pts.max(axis=0)
        return (x0 - L, x1 + 4 * L, y0 - 1.5 * L, y1 + 1.5 * L)

    # ------------------------------------------------------------------ physique / CL
    def _nu_mode_changed(self):
        by_nu = self.nu_mode.currentData() == "nu"
        self.nu_edit.setEnabled(by_nu)
        self.re_edit.setEnabled(not by_nu)
        if not self._syncing:
            self._form_changed()

    def _energy_toggled(self, on=None):
        on = self.energy_on.isChecked()
        for _, wdg in self.energy_fields:
            wdg.setEnabled(on)
        if self._syncing:
            return
        if on:
            self.cfg.setdefault("energy", {})
            self.binder.store(self.cfg)
        else:
            self.cfg.pop("energy", None)
        self._toml_timer.start()

    def _visc_changed(self, *_):
        model = self.visc_combo.currentData() or "newtonian"
        used = set(VISCOSITY_PARAMS.get(model, ()))
        if model != "newtonian":
            used |= {"nu_max", "nu_min"}
        for key, wdg in self.visc_fields.items():
            wdg.setEnabled(key in used)
        if not self._syncing:
            self._form_changed()

    def _fill_scalar_table(self):
        self._scalar_sync = True
        sc = self.cfg.get("scalars", {}) or {}
        self.scalar_table.setRowCount(0)
        self.scalar_table.setRowCount(len(sc))
        for r, (name, sp) in enumerate(sc.items()):
            D = sp.get("diffusivity", "")
            if D == "" and "schmidt" in sp:
                D = f"Sc={sp['schmidt']}"
            for c, val in enumerate((name, D, sp.get("Sc_t", ""), sp.get("source", ""),
                                     sp.get("initial", ""))):
                self.scalar_table.setItem(r, c, QTableWidgetItem(str(val)))
        self._scalar_sync = False

    def _scalar_add(self):
        sc = self.cfg.setdefault("scalars", {})
        i = 1
        while f"c{i}" in sc:
            i += 1
        sc[f"c{i}"] = {"diffusivity": 1e-3}
        self._fill_scalar_table()
        self._toml_timer.start()

    def _scalar_remove(self):
        r = self.scalar_table.currentRow()
        sc = self.cfg.get("scalars", {})
        if r < 0 or not sc:
            return
        sc.pop(list(sc)[r], None)
        if not sc:
            self.cfg.pop("scalars", None)
        self._fill_scalar_table()
        self._toml_timer.start()

    def _scalars_changed(self):
        if getattr(self, "_scalar_sync", False):
            return
        old = self.cfg.get("scalars", {}) or {}
        new = {}

        def cell(r, c):
            it = self.scalar_table.item(r, c)
            txt = it.text().strip() if it else ""
            if not txt:
                return None
            try:
                return float(txt)
            except ValueError:
                return txt                          # expression en x, y
        for r, prev in enumerate(old.values()):
            name = cell(r, 0)
            if not isinstance(name, str):
                name = f"c{r + 1}"
            sp = {k: v for k, v in prev.items() if k in ("scheme",)}
            D = cell(r, 1)
            if isinstance(D, str) and D.lower().startswith("sc="):
                sp["schmidt"] = float(D[3:])
            elif D is not None:
                sp["diffusivity"] = D
            for key, c in (("Sc_t", 2), ("source", 3), ("initial", 4)):
                val = cell(r, c)
                if val is not None:
                    sp[key] = val
            new[name] = sp
        self.cfg["scalars"] = new
        self._toml_timer.start()

    def _mode_changed(self):
        steady = self.mode_combo.currentData() == "steady"
        self.box_steady.setVisible(steady)
        self.box_transient.setVisible(not steady)
        self.box_steady.setEnabled(steady)
        self.box_transient.setEnabled(not steady)

    def _patch_names(self):
        if self.mesh is not None:
            return [n for n, t, _ in self.mesh.all_boundary_patches()
                    if n not in self._periodic_names()]
        return list(self.cfg.get("boundary", {}))

    def _periodic_names(self):
        if self.mesh is None:
            return set()
        return {p for pair in self.mesh.periodic_pairs for p in pair[:2]}

    @staticmethod
    def _guess_bc(name, ptype):
        n = name.lower()
        if "farfield" in n or "far" in n:
            return {"type": "farfield", "U": [1.0, 0.0]}
        if "inlet" in n or "entree" in n or "inflow" in n:
            return {"type": "inlet", "U": [1.0, 0.0]}
        if "outlet" in n or "sortie" in n or "outflow" in n:
            return {"type": "outlet", "p": 0.0}
        if n in ("axis", "axe"):
            return {"type": "axis"}
        if "sym" in n or ptype == "symmetry":
            return {"type": "symmetry"}
        if n in ("top", "bottom"):
            return {"type": "symmetry"}
        return {"type": "wall"}

    def _sync_bc_with_mesh(self):
        bnd = self.cfg.setdefault("boundary", {})
        names = self._patch_names()
        types = {n: t for n, t, _ in self.mesh.all_boundary_patches()}
        for n in names:
            if n not in bnd:
                bnd[n] = self._guess_bc(n, types.get(n, "patch"))
                self.log(f"Condition par défaut pour '{n}' : {bnd[n]['type']} (à vérifier)")
        for n in list(bnd):
            if n not in names:
                bnd.pop(n)
        self._fill_bc_table()
        self._refresh_toml()

    def _fill_bc_table(self):
        self._bc_sync = True
        bnd = self.cfg.get("boundary", {})
        names = self._patch_names() or list(bnd)
        self.bc_table.setRowCount(0)             # supprime aussi les anciennes listes
        self.bc_table.setRowCount(len(names))
        for r, n in enumerate(names):
            spec = bnd.get(n, {"type": "wall"})
            it = QTableWidgetItem(n)
            it.setFlags(it.flags() & ~Qt.ItemIsEditable)
            self.bc_table.setItem(r, 0, it)
            cb = combo([(t, BC_LABELS[t]) for t in BC_TYPES])
            set_combo(cb, spec.get("type", "wall"))
            cb.currentIndexChanged.connect(lambda *_: self._bc_changed())
            self.bc_table.setCellWidget(r, 1, cb)
            U = spec.get("U", ["", ""])
            self.bc_table.setItem(r, 2, QTableWidgetItem(str(U[0])))
            self.bc_table.setItem(r, 3, QTableWidgetItem(str(U[1])))
            self.bc_table.setItem(r, 4, QTableWidgetItem(str(spec.get("p", spec.get("p0", "")))))
            self.bc_table.setItem(r, 5, QTableWidgetItem(str(spec.get("T", ""))))
            self.bc_table.setItem(r, 6, QTableWidgetItem(str(spec.get("q", ""))))
            self.bc_table.setItem(r, 7, QTableWidgetItem(str(spec.get("flow_rate", ""))))
            self.bc_table.setItem(r, 8, QTableWidgetItem("; ".join(
                f"{k}={v}" for k, v in (spec.get("scalars") or {}).items())))
        self._bc_sync = False

    def _bc_changed(self):
        if getattr(self, "_bc_sync", False):
            return
        bnd = {}

        def num(txt):
            txt = txt.strip()
            if not txt:
                return None
            try:
                return float(txt)
            except ValueError:
                return txt                          # expression en x, y
        old = self.cfg.get("boundary", {})
        for r in range(self.bc_table.rowCount()):
            name = self.bc_table.item(r, 0).text()
            t = self.bc_table.cellWidget(r, 1).currentData()
            # clés sans colonne (profil de débit, flux de scalaires…) conservées
            spec = {k: v for k, v in old.get(name, {}).items()
                    if k in ("profile", "scalar_flux")}
            spec["type"] = t
            ux = num(self.bc_table.item(r, 2).text() if self.bc_table.item(r, 2) else "")
            uy = num(self.bc_table.item(r, 3).text() if self.bc_table.item(r, 3) else "")
            Q = num(self.bc_table.item(r, 7).text() if self.bc_table.item(r, 7) else "")
            if t == "inlet" and Q is not None:
                spec["flow_rate"] = Q
            elif t in ("inlet", "farfield") or (t == "wall" and (ux or uy)):
                spec["U"] = [ux if ux is not None else 0.0, uy if uy is not None else 0.0]
            p = num(self.bc_table.item(r, 4).text() if self.bc_table.item(r, 4) else "")
            if t in ("outlet", "farfield") and p is not None:
                spec["p"] = p
            elif t == "pressure_inlet":
                spec["p0"] = p if p is not None else 0.0
            T = num(self.bc_table.item(r, 5).text() if self.bc_table.item(r, 5) else "")
            q = num(self.bc_table.item(r, 6).text() if self.bc_table.item(r, 6) else "")
            if T is not None and t in ("wall", "inlet", "farfield", "pressure_inlet"):
                spec["T"] = T
            elif q is not None and t == "wall":
                spec["q"] = q
            sc = self.bc_table.item(r, 8).text() if self.bc_table.item(r, 8) else ""
            vals = {}
            for part in sc.replace(",", ";").split(";"):
                if "=" in part:
                    k, v = (x.strip() for x in part.split("=", 1))
                    vals[k] = num(v)
            if vals:
                spec["scalars"] = vals
            bnd[name] = spec
        self.cfg["boundary"] = bnd
        self._toml_timer.start()

    # ================================================================== exécution
    def _start(self, job, on_done, message, on_progress=None):
        if self.thread is not None:
            QMessageBox.information(self, "Occupé", "Un calcul est déjà en cours.")
            return
        self.thread = QThread(self)
        self.worker = Worker(job)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        # méthodes liées (pas de lambda) : Qt les exécute dans le fil de l'interface
        self._on_done = on_done
        self.worker.log.connect(self.log)
        if on_progress:
            self.worker.progress.connect(on_progress)
        self.worker.done.connect(self._finish)
        self.worker.failed.connect(self._failed)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self.statusBar().showMessage(message)
        self.run_btn.setEnabled(False)
        self.continue_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.thread.start()

    def _cleanup(self):
        if self.thread is not None:
            self.thread.quit()
            self.thread.wait()
        self.thread = self.worker = None
        self.progress.setVisible(False)
        self.run_btn.setEnabled(True)
        self.continue_btn.setEnabled(True)
        self.sweep_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _finish(self, res):
        on_done = self._on_done
        self._cleanup()
        self.statusBar().showMessage("Terminé", 5000)
        on_done(res)

    def _failed(self, msg):
        self.errors.append(msg)
        self._cleanup()
        self.log(msg)
        self.statusBar().showMessage("Échec", 5000)
        if not self.quiet:
            QMessageBox.warning(self, "Erreur", msg.split("\n\n")[0])

    def stop(self):
        if self.worker is not None:
            self.worker.stop_requested = True
            self.statusBar().showMessage("Arrêt demandé…")

    def run_sweep_gui(self):
        from ..fv2d.sweep import parse_values
        key = self.sweep_param.currentData()
        try:
            values = parse_values(self.sweep_values.text())
        except ValueError as exc:
            QMessageBox.warning(self, "Valeurs", f"Valeurs invalides : {exc}")
            return
        if not values:
            return
        self._store_forms()
        cfg = copy.deepcopy(self.cfg)
        cfg.setdefault("output", {})["directory"] = str(self.out_dir())
        if cfg.get("solver", {}).get("mode", "steady") != "steady":
            QMessageBox.information(self, "Balayage", "Le balayage fonctionne aussi en "
                                    "instationnaire, mais chaque point est alors long. "
                                    "Il est conseillé de passer en stationnaire.")
        mesh, base, cont = self.mesh, self.base_dir(), self.sweep_cont.isChecked()
        self.sweep_rows, self._sweep_key, self._sweep_n = [], key, len(values)

        def job(worker):
            from ..fv2d.sweep import run_sweep
            rows = run_sweep(cfg, key, values, base_dir=base, out_dir=cfg["output"]["directory"],
                             continuation=cont, verbose=True, plot=True,
                             callback=lambda s, n: worker.stop_requested, mesh=mesh,
                             on_point=lambda row: worker.report(row, force=True))
            return rows
        self.nav.setCurrentRow(6)
        self.tabs.setCurrentIndex(0)
        self.canvas.message(f"Balayage : {len(values)} calculs…")
        self.sweep_btn.setEnabled(False)
        self._start(job, self._sweep_done, "Balayage en cours…", self._sweep_progress)

    def _sweep_progress(self, row):
        from ..fv2d.sweep import plot_sweep
        self.sweep_rows.append(row)
        n = len(self.sweep_rows)
        self.progress.setRange(0, self._sweep_n)
        self.progress.setValue(n)
        self.run_info.setText(f"Point {n}/{self._sweep_n} : {self._sweep_key} = "
                              f"{row[self._sweep_key]:g} "
                              f"({'convergé' if row['converged'] else 'NON convergé'})")
        self.canvas.fig.clear()
        plot_sweep(self.sweep_rows, self._sweep_key, None, fig=self.canvas.fig)
        self.canvas.draw()

    def _sweep_done(self, rows):
        from ..fv2d.sweep import plot_sweep
        self.sweep_btn.setEnabled(True)
        self.sweep_rows = list(rows)             # tableau complet (référence)
        if rows:
            self.canvas.fig.clear()
            plot_sweep(rows, self._sweep_key, None, fig=self.canvas.fig)
            self.canvas.draw()
        cols = [k for k in (rows[0] if rows else {}) if not k.startswith(("Cd_p", "Cd_v"))]
        lines = ["\t".join(cols)] + ["\t".join(f"{r.get(c, ''):.5g}" if isinstance(
            r.get(c), float) else str(r.get(c, "")) for c in cols) for r in rows]
        self.summary_view.setPlainText("\n".join(lines))
        self.run_info.setText(f"Balayage terminé : {len(rows)} points. Tableau : "
                              f"{self.out_dir() / 'balayage.csv'}")
        self.log(f"Balayage écrit dans {self.out_dir()}")

    def _browse_restart(self):
        path, _ = QFileDialog.getOpenFileName(self, "Fichier de reprise", str(self.out_dir()),
                                              "Reprise microrans (*.npz)")
        if path:
            self.restart_path.setText(path)
            self._store_forms()

    def continue_2d(self):
        ck = self.out_dir() / "checkpoint.npz"
        if not ck.is_file():
            QMessageBox.information(self, "Rien à continuer",
                                    f"Aucune sauvegarde dans {self.out_dir()} : lancez d'abord "
                                    "un calcul (ou choisissez un fichier dans Conditions "
                                    "initiales).")
            return
        self.run_2d(restart=ck)

    def run_2d(self, restart=None):
        self._store_forms()
        cfg = copy.deepcopy(self.cfg)
        cfg.setdefault("output", {})["directory"] = str(self.out_dir())
        if restart:
            cfg.setdefault("initial", {})["restart"] = str(restart)
        self._it0 = None
        mesh = self.mesh
        base = self.base_dir()
        self.history = []
        steady = cfg.get("solver", {}).get("mode", "steady") == "steady"
        self._run_meta = (steady, cfg.get("solver", {}).get("max_iter", 3000),
                          cfg.get("solver", {}).get("t_end", 1.0))

        def job(worker):
            from ..fv2d.case import run_case

            def cb(s, n):
                # enregistrement du pas courant (efforts en axes écoulement, calculés par
                # run_case)
                src = s.history if steady else s.series
                rec = {k: v for k, v in src[-1].items() if k.startswith(
                    ("iteration", "time", "dt", "Cd_", "Cl_")) or steady} if src else {}
                worker.report(rec)
                return worker.stop_requested
            return run_case(cfg, base_dir=base, out_dir=cfg["output"]["directory"], verbose=True,
                            plot=True, callback=cb, mesh=mesh, return_solver=True)
        self.nav.setCurrentRow(6)
        self.tabs.setCurrentIndex(0)
        self.canvas.message("Calcul en cours…")
        self._start(job, self._run_done, "Calcul en cours…", self._on_progress)

    def _on_progress(self, rec):
        self.history.append(rec)
        steady, max_iter, t_end = self._run_meta
        self.progress.setRange(0, 1000)
        if steady and "iteration" in rec:
            if self._it0 is None:
                self._it0 = rec["iteration"] - 1          # suite d'un calcul repris
            self.progress.setValue(int(1000 * (rec["iteration"] - self._it0) / max(max_iter, 1)))
            self.run_info.setText(f"itération {rec['iteration']} — " + ", ".join(
                f"{k} {v:.1e}" for k, v in rec.items() if k != "iteration"))
        elif "time" in rec:
            self.progress.setValue(int(1000 * rec["time"] / max(t_end, 1e-300)))
            self.run_info.setText(f"t = {rec['time']:.4g} (Δt = {rec['dt']:.3g})")
        self._plot_history()

    def _plot_history(self, final=False):
        hist = self.history
        if final and self.solver is not None and self.solver.history:
            hist = self.solver.history
        if not hist:
            return
        ax = self.canvas.axes()
        colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
        if "iteration" in hist[-1]:
            it = [h["iteration"] for h in hist]
            keys = [k for k in hist[-1] if k != "iteration"]
            for i, k in enumerate(keys):
                ax.semilogy(it, [max(h.get(k, np.nan), 1e-300) for h in hist],
                            color=colors[i % len(colors)], label=k)
            ax.set(xlabel="itération", ylabel="résidu normalisé", title="Convergence")
        else:
            t = [h["time"] for h in hist]
            keys = [k for k in hist[-1] if k.startswith(("Cd_", "Cl_"))]
            for i, k in enumerate(keys):
                ax.plot(t, [h.get(k, np.nan) for h in hist], color=colors[i % len(colors)],
                        label=k)
            ax.set(xlabel="t", ylabel="coefficient", title="Efforts")
        ax.grid(True)
        if len(keys) > 1:
            ax.legend(loc="best")
        self.canvas.draw()
        self.tabs.setCurrentIndex(0)

    def _run_done(self, res):
        summary, solver = res
        self.summary, self.solver = summary, solver
        self.mesh = solver.mesh
        text = json.dumps(summary, indent=2, ensure_ascii=False)
        self.summary_view.setPlainText(text)
        fields = solver.fields()
        self.field_combo.clear()
        for k in ["U_mag", "Ux", "Uy", "p", "vorticity"] + [k for k in fields if k not in (
                "U", "U_mag", "p", "vorticity")]:
            self.field_combo.addItem(FIELD_LABELS.get(k, k), k)
        self.wall_combo.clear()
        for p in self.mesh.patches:
            if p.type == "wall":
                self.wall_combo.addItem(p.name, p.name)
        conv = summary.get("converged")
        self.run_info.setText(("Calcul terminé" if conv in (True, None) else "Arrêté / non convergé")
                              + f" en {summary.get('wall_time_s', 0):.1f} s. Résultats : "
                              f"{self.out_dir()}")
        self.log(f"Résultats écrits dans {self.out_dir()}")
        self.nav.setCurrentRow(7)
        self.plot_field()

    def plot_field(self):
        from ..fv2d.post import _body_size
        from ..mesh2d.plot import plot_field, plot_mesh
        s = self.solver
        if s is None:
            self.canvas.message("Aucun résultat : lancez d'abord un calcul.")
            return
        key = self.field_combo.currentData() or "U_mag"
        f = s.fields()
        if key == "Ux":
            val = s.U[:, 0]
        elif key == "Uy":
            val = s.U[:, 1]
        elif key == "vorticity":
            g = s.grad_U(s.U)
            val = g[:, 1, 0] - g[:, 0, 1]
        else:
            val = f[key]
        zoom = self._zoom() if self.zoom_check.isChecked() else None
        vmin = vmax = None
        if key == "vorticity":
            _, L = _body_size(s)
            far = s.mesh.wall_distance > 0.2 * L if L else np.ones(len(val), dtype=bool)
            lim = float(np.percentile(np.abs(val[far] if np.any(far) else val), 99))
            vmin, vmax = -lim, lim
        ax = self.canvas.axes()
        cmap = self.cmap_combo.currentData()
        if key in ("vorticity", "p") and cmap == "viridis":
            cmap = "RdBu_r"
        mirror = None
        if s.axisymmetric:                           # image miroir par rapport à l'axe
            mirror = -1 if key in ("Uy", "vorticity") else 1
        plot_field(s.mesh, val, ax=ax, cmap=cmap, zoom=zoom, vmin=vmin, vmax=vmax,
                   title=FIELD_LABELS.get(key, key),
                   vectors=s.U if self.vec_check.isChecked() else None, mirror=mirror)
        if self.mesh_check.isChecked():
            plot_mesh(s.mesh, ax=ax, zoom=zoom, linewidth=0.15, show_patches=False)
            ax.collections[-1].set_facecolor("none")
            ax.set_title(FIELD_LABELS.get(key, key))
        self.canvas.draw()
        self.tabs.setCurrentIndex(0)

    def open_animation(self):
        path = (self.summary or {}).get("animation")
        if not path or not Path(path).is_file():
            QMessageBox.information(self, "Animation", "Pas d'animation : choisir une grandeur "
                                    "dans Numérique → Sorties → Animation, puis lancer un calcul "
                                    "instationnaire.")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def plot_line(self):
        from ..fv2d.sampling import Sampler, line_points
        s = self.solver
        if s is None:
            self.canvas.message("Aucun résultat : lancez d'abord un calcul.")
            return
        key = self.field_combo.currentData() or "U_mag"
        dist, pts = line_points(self.line_start.value(), self.line_end.value(), 300)
        smp = Sampler(s, pts)
        vals = smp.sample([key]).get(key)
        if vals is None:                            # grandeur sans reconstruction : cellule
            f = s.fields()
            if key == "vorticity":
                g = s.grad_U(s.U)
                f[key] = g[:, 1, 0] - g[:, 0, 1]
            cell = np.asarray(f[key])
            vals = np.where(smp.ok, cell[smp._c], np.nan)
        ax = self.canvas.axes()
        ax.plot(dist, vals, color="#2a78d6", lw=1.6)
        ax.set(xlabel="abscisse le long de la ligne", ylabel=FIELD_LABELS.get(key, key),
               title=f"Profil de {FIELD_LABELS.get(key, key)}")
        ax.grid(True, alpha=0.3)
        if not smp.ok.all():
            ax.text(0.01, 0.01, "points hors du domaine non tracés", transform=ax.transAxes,
                    fontsize=8, color="#666666")
        self.canvas.draw()
        self.tabs.setCurrentIndex(0)

    def plot_wall(self):
        s = self.solver
        name = self.wall_combo.currentData()
        if s is None or not name:
            return
        xf, tau, yp = s.wall_shear(name)
        q = 0.5 * s.U_ref ** 2
        pb = s.boundary_p(s.p)[s.patch_slices[name]]
        what = self.wall_q.currentData()
        if what in ("q", "Tw"):
            if s.energy is None:
                self.canvas.message("Pas de thermique dans ce calcul.")
                return
            Tw, qw = s.wall_heat_flux(name)
            y = qw if what == "q" else Tw
        else:
            y = {"Cp": pb / q, "Cf": tau / q, "yplus": yp}[what]
        o = np.argsort(xf[:, 0])
        ax = self.canvas.axes()
        ax.plot(xf[o, 0], y[o], "o", ms=3, color="#2a78d6")
        ax.set(xlabel="x", ylabel=what, title=f"{what} sur '{name}'")
        if what == "Cp":
            ax.invert_yaxis()
        ax.grid(True)
        self.canvas.draw()

    # ------------------------------------------------------------------ 1D
    def run_1d(self):
        models = [k for k, cb in self.c1d_models.items() if cb.isChecked()]
        if not models:
            QMessageBox.information(self, "Canal 1D", "Cochez au moins un modèle.")
            return
        mode = self.c1d_mode.currentData()
        out = results_root() / f"canal1d_{mode}"
        argv = [mode, "-m", *models, "--re-tau", str(self.c1d_retau.value() or 395),
                "--n-cells", str(int(self.c1d_n.value() or 192)),
                "--y1plus", str(self.c1d_y1.value() or 0.2), "-o", str(out)]
        if mode == "urans":
            argv += ["--omega-plus", str(self.c1d_omega.value() or 0.01),
                     "--amplitude", str(self.c1d_amp.value() or 10),
                     "--steps-per-period", str(int(self.c1d_steps.value() or 64)),
                     "--scheme", self.c1d_scheme.currentData()]

        def job(worker):
            from ..cli import main
            print("microrans " + " ".join(argv))
            code = main(argv)
            return out, models, mode, code
        self.c1d_status.setText("Calcul en cours…")
        self._start(job, self._run_1d_done, "Calcul 1D en cours…")

    def _run_1d_done(self, res):
        out, models, mode, code = res
        self.c1d_status.setText(f"Terminé (code {code}). Résultats : {out}")
        imgs = sorted(out.glob("comparaison*.png")) or sorted(out.glob(f"{models[0]}/*.png"))
        if not imgs:
            self.canvas.message("Pas de figure produite (voir le journal).")
            return
        import matplotlib.image as mpimg
        ax = self.canvas.axes()
        ax.imshow(mpimg.imread(str(imgs[0])))
        ax.axis("off")
        self.canvas.draw()
        self.tabs.setCurrentIndex(0)


def QLineEdit_(text=""):
    from PySide6.QtWidgets import QLineEdit
    return QLineEdit(text)


# ============================================================================ lancement
def _selftest(win: MainWindow, app, shot: str | None) -> int:
    """Test de fumée (CI, exécutable) : cavité 16×16, quelques itérations, tracés."""
    import time
    win.open_case(examples_dir() / "cavite_re100.toml")
    win.cfg["mesh"].update(nx=16, ny=16)
    win.cfg["solver"].update(max_iter=60)
    win.cfg["output"]["directory"] = str(results_root() / "selftest")
    win.load_cfg(win.cfg)
    win.quiet = True
    errors = win.errors
    its = []
    for action in (win.generate_mesh, win.run_2d, win.continue_2d):
        action()
        t0 = time.time()
        while win.thread is not None and time.time() - t0 < 300:
            app.processEvents()
            time.sleep(0.02)
        if errors:
            print("SELFTEST ÉCHEC :", errors[0])
            return 1
        if win.summary is not None:
            its.append(win.summary["iterations"])
    # « Continuer » : reprise exacte, 60 itérations de plus
    ok = (win.solver is not None and len(its) == 2 and its[0] > 10 and its[1] > its[0]
          and win.summary.get("restart", {}).get("mode") == "exact")
    # balayage (2 viscosités) : tableau et figure
    set_combo(win.sweep_param, "physics.nu")
    win.sweep_values.setText("0.01, 0.02")
    win.run_sweep_gui()
    t0 = time.time()
    while win.thread is not None and time.time() - t0 < 300:
        app.processEvents()
        time.sleep(0.02)
    ok = ok and not errors and len(win.sweep_rows) == 2 and (
        win.out_dir() / "balayage.csv").is_file()
    if errors:
        print("SELFTEST ÉCHEC :", errors[0])
    for i in range(win.field_combo.count()):
        win.field_combo.setCurrentIndex(i)
        win.plot_field()
        win.line_start.set_value((0.5, 0.0))
        win.line_end.set_value((0.5, 1.0))
        win.plot_line()
    app.processEvents()
    if shot:
        win.grab().save(shot)
    print("SELFTEST", "OK" if ok else "ÉCHEC", its, len(win.sweep_rows))
    return 0 if ok else 1


def run(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName(APP_NAME)
    app.setStyle("Fusion")
    icon = Path(__file__).with_name("icon.png")
    if icon.exists():
        from PySide6.QtGui import QIcon
        app.setWindowIcon(QIcon(str(icon)))
    win = MainWindow()
    if "--selftest" in argv:
        shot = argv[argv.index("--selftest") + 1] if len(argv) > argv.index("--selftest") + 1 \
            else None
        return _selftest(win, app, shot)
    win.show()
    for a in argv:
        if a.endswith((".toml", ".json")) and Path(a).exists():
            win.open_case(Path(a))
    return app.exec()
