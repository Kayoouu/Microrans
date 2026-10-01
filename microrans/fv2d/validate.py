"""Vérification d'un fichier de cas 2D avant calcul.

Clé ou section inconnue (faute de frappe) : avertissement, avec la clé connue la plus
proche (« max_iters » → « max_iter ? », « [solveur] » → « [solver] ? »). Clé lue par l'autre
solveur seulement (incompressible / compressible) ou sans effet pour le type de maillage
choisi : avertissement. Structure impossible (liste attendue, table attendue) : erreur.

SCHEMA est la liste unique des clés connues ; un test vérifie qu'elle couvre tous les
réglages des solveurs (Settings, CompressibleSettings) et tous les exemples fournis.
"""
from __future__ import annotations

import difflib
import warnings
from dataclasses import dataclass, field

INC, COMP = "incompressible", "compressible"


class CaseWarning(UserWarning):
    """Clé du fichier de cas inconnue ou sans effet."""


@dataclass(frozen=True)
class Key:
    doc: str
    only: str | None = None                 # INC | COMP : lue par un seul solveur
    types: tuple = ()                       # [mesh] : types de maillage qui l'utilisent


@dataclass(frozen=True)
class Table:
    keys: dict = field(default_factory=dict)   # nom -> Key | Table ; "*" : noms libres
    doc: str = ""
    only: str | None = None
    types: tuple = ()
    many: bool = False                      # liste de tables : [[bodies]], [[porous]]…
    free: bool = False                      # contenu libre (non vérifié ici)


K = Key
RECT, BLOCKS, OGRID, TRI, HYB, FILE = ("rectangle", "blocks", "ogrid", "unstructured",
                                       "hybrid", "file")
_NAMES = Table({"left": K("nom du côté x = x0"), "right": K("nom du côté x = x1"),
                "bottom": K("nom du côté y = y0"), "top": K("nom du côté y = y1")},
               "noms des frontières du rectangle")
# géométries : [domain], [[bodies]], [[mesh.refinements]] shape (geometry.shape_from_dict)
_SHAPE = {
    "type": K("circle | rectangle | ellipse | polygon | naca | spline | file"),
    "name": K("nom du corps = nom de la frontière (condition aux limites)"),
    "patch_type": K("type de frontière (défaut wall pour un corps)"),
    "center": K("centre [x, y] (circle, ellipse ; défaut [0, 0])"),
    "radius": K("rayon (circle)"),
    "x0": K("rectangle : x min"), "y0": K("rectangle : y min"),
    "x1": K("rectangle : x max"), "y1": K("rectangle : y max"),
    "names": _NAMES,
    "a": K("ellipse : demi-axe selon x"), "b": K("ellipse : demi-axe selon y"),
    "points": K("polygon, spline : [[x, y], …]"),
    "code": K("naca : 4 chiffres (défaut \"0012\")"),
    "chord": K("naca : corde (défaut 1)"),
    "n": K("naca, spline : nombre de points du contour"),
    "closed_te": K("naca : bord de fuite fermé (défaut true)"),
    "path": K("file : fichier de contour (.dat Selig/Lednicer, .csv, .dxf)"),
    "sharp_angle": K("file : angle (°) au-delà duquel un sommet est un coin (défaut 60)"),
    "angle": K("rotation (°, sens trigonométrique)"),
    "incidence": K("incidence (°) : rotation de −incidence (profil cabré)"),
    "rotation_center": K("centre de rotation [x, y] (défaut [0, 0])"),
    "scale": K("facteur d'échelle (défaut 1)"),
    "translate": K("translation [dx, dy]"),
}
_TRI_HYB = (TRI, HYB)

SCHEMA = Table({
    "mesh": Table({
        "type": K("rectangle | blocks | ogrid | unstructured | hybrid | file "
                  "(défaut unstructured)"),
        "preset": K("maillage prédéfini (cavity, channel, cylinder-ogrid…) : remplace type"),
        "cut_axis": K("garde la moitié y > 0 (axisymétrique autour d'un corps)"),
        "x0": K("x min", types=(RECT,)), "x1": K("x max", types=(RECT,)),
        "y0": K("y min", types=(RECT,)), "y1": K("y max", types=(RECT,)),
        "nx": K("nombre de mailles selon x", types=(RECT,)),
        "ny": K("nombre de mailles selon y", types=(RECT,)),
        "grading": K("resserrement [gx, gy] (rapport dernière / première maille)",
                     types=(RECT,)),
        "names": Table(_NAMES.keys, _NAMES.doc, types=(RECT,)),
        "patch_types": K("types des frontières {nom = \"wall\" | \"patch\" | …}",
                         types=(RECT, FILE)),
        "periodic": K("paires de frontières périodiques [[\"a\", \"b\"], …]",
                      types=(RECT, BLOCKS)),
        "vertices": K("sommets [[x, y], …]", types=(BLOCKS,)),
        "blocks": Table({"vertices": K("4 indices de sommets (sens trigonométrique)"),
                         "cells": K("[nx, ny]"), "grading": K("[gx, gy]")},
                        "blocs", types=(BLOCKS,), many=True),
        "edges": Table({"type": K("arc | polyline | spline"),
                        "vertices": K("[i, j] : arête courbe entre deux sommets"),
                        "point": K("arc : point intermédiaire [x, y]"),
                        "points": K("polyline, spline : points intermédiaires")},
                       "arêtes courbes", types=(BLOCKS,), many=True),
        "patches": Table({"*": Table({"type": K("wall | patch | symmetry | empty"),
                                      "faces": K("[[i, j], …] : paires de sommets")})},
                         "frontières nommées", types=(BLOCKS,)),
        "n_around": K("mailles autour du corps (défaut 128)", types=(OGRID,)),
        "n_radial": K("mailles dans la direction radiale (défaut 64)", types=(OGRID,)),
        "farfield_radius": K("rayon du champ lointain (défaut 20)", types=(OGRID,)),
        "first_height": K("hauteur de la 1re maille à la paroi (défaut 1e-3)",
                          types=(OGRID,)),
        "center": K("centre du maillage en O (défaut : centre du corps)", types=(OGRID,)),
        "h_max": K("taille maximale des triangles (défaut 1)", types=_TRI_HYB),
        "h_surface": K("taille des mailles sur les corps (défaut 0.05)", types=_TRI_HYB),
        "growth": K("croissance de la taille avec la distance au corps (défaut 0.2)",
                    types=_TRI_HYB),
        "max_iter": K("itérations du mailleur DistMesh (défaut 300)", types=_TRI_HYB),
        "refinements": Table({"shape": Table(_SHAPE, "zone (même syntaxe que les corps)"),
                              "h": K("taille des mailles dans la zone"),
                              "growth": K("croissance hors de la zone (défaut 0.2)")},
                             "zones de raffinement", types=_TRI_HYB, many=True),
        "layers": Table({"n": K("nombre de couches (défaut 10)"),
                         "first_height": K("hauteur de la 1re couche (défaut 1e-3)"),
                         "ratio": K("rapport de croissance (défaut 1.2)")},
                        "couches de quadrilatères aux parois", types=(HYB,)),
        "path": K("fichier .msh (Gmsh) ou .su2", types=(FILE,)),
    }, "maillage"),
    "domain": Table(_SHAPE, "domaine extérieur (défaut : rectangle [-10, 30] × [-10, 10])",
                    types=_TRI_HYB),
    "bodies": Table(_SHAPE, "corps (obstacles)", types=(OGRID, TRI, HYB), many=True),
    "physics": Table({
        "compressible": K("true : solveur compressible (grandeurs SI, section [flow])"),
        "nu": K("viscosité cinématique ν (m²/s)", INC),
        "reynolds": K("nombre de Reynolds (ν = U_ref L_ref / Re) — au lieu de nu", INC),
        "reference_velocity": K("vitesse de référence U_ref (défaut 1)", INC),
        "reference_length": K("longueur de référence L_ref (défaut 1)"),
        "reference_area": K("axisymétrique : aire de référence de Cd (défaut π L²/4)", INC),
        "model": K("laminar | sa | ke | kw | sst | sst_gamma (défaut laminar)", INC),
        "model_options": Table({}, "options du modèle (sa : ft2 ; sst_gamma : "
                               "kato_launder)", INC, free=True),
        "body_force": K("force volumique [fx, fy] (conduite périodique)", INC),
        "angle_of_attack": K("incidence (°) de l'écoulement amont"),
        "axisymmetric": K("true : axisymétrique (x = axe, y = rayon)"),
        "swirl": K("true : rotation propre (axisymétrique)", INC),
        "viscosity": Table({
            "model": K("newtonian | power_law | carreau | cross | herschel_bulkley | "
                       "bingham | casson"),
            "K": K("indice de consistance"), "n": K("indice de comportement"),
            "nu0": K("viscosité à cisaillement nul"), "nu_inf": K("viscosité à "
                                                                  "cisaillement infini"),
            "lambda": K("carreau : temps caractéristique"), "a": K("carreau : exposant "
                                                                   "(défaut 2)"),
            "m": K("cross : temps caractéristique"), "tau_y": K("contrainte seuil"),
            "nu_min": K("borne inférieure de ν"), "nu_max": K("borne supérieure de ν"),
            "relax": K("sous-relaxation de ν (défaut 0.5)"),
            "gamma_ref": K("taux de cisaillement de référence (défaut U_ref / L_ref)")},
            "fluide non newtonien", INC),
    }, "physique"),
    "flow": Table({
        "mach": K("nombre de Mach amont"), "velocity": K("vitesse [u, v] (m/s)"),
        "pressure": K("pression (Pa)"), "temperature": K("température (K)"),
        "density": K("masse volumique (kg/m³)"), "gamma": K("γ (défaut 1.4)"),
        "gas_constant": K("R (J/kg/K, défaut 287.058)"), "prandtl": K("Pr (défaut 0.72)"),
        "viscosity": K("inviscid | constant | sutherland"), "mu": K("viscosité (Pa·s)"),
        "reynolds": K("Re = ρ U L_ref / μ (au lieu de mu)"),
        "T_ref": K("Sutherland : température de référence (défaut 273.15 K)"),
        "sutherland_S": K("Sutherland : constante S (défaut 110.4 K)"),
    }, "écoulement amont (compressible)", COMP),
    "initial": Table({
        "U": K("vitesse initiale [ux, uy] (valeurs ou formules en x, y)"),
        "perturbation": K("amplitude d'un tourbillon initial (déclenche le lâcher)", INC),
        "perturbation_center": K("centre du tourbillon initial (défaut [1.5, 0])", INC),
        "restart": K("fichier checkpoint.npz de reprise"),
        "restart_mode": K("exact (défaut) | fields"),
        "restart_shift_U": K("décalage de vitesse à la reprise", INC),
        "rho": K("masse volumique (formule possible)", COMP),
        "density": K("= rho", COMP), "p": K("pression (formule possible)", COMP),
        "pressure": K("= p", COMP), "T": K("température (formule possible)", COMP),
        "temperature": K("= T", COMP), "mach": K("Mach initial (direction amont)", COMP),
    }, "état initial"),
    "turbulence": Table({"intensity": K("intensité turbulente amont (défaut 0.001)"),
                         "viscosity_ratio": K("ν_t/ν amont (défaut 0.1)")},
                        "turbulence amont", INC),
    "energy": Table({
        "enabled": K("false : thermique désactivée"), "Pr": K("nombre de Prandtl"),
        "Pr_t": K("Prandtl turbulent (défaut 0.85)"),
        "beta": K("dilatation thermique (Boussinesq)"), "T_ref": K("température de "
                                                                   "référence"),
        "gravity": K("gravité [gx, gy]"), "T0": K("température initiale"),
        "delta_T": K("écart de température de référence (nombre de Nusselt)"),
        "source": K("source de chaleur (valeur ou formule)")}, "thermique", INC),
    "scalars": Table({"*": Table({
        "diffusivity": K("diffusivité D"), "schmidt": K("Sc = ν / D (au lieu de "
                                                         "diffusivity)"),
        "Sc_t": K("Schmidt turbulent (défaut 0.7)"), "source": K("source (formule)"),
        "initial": K("valeur initiale"), "scheme": K("upwind | linearUpwind")})},
        "scalaires passifs [scalars.<nom>]", INC),
    "porous": Table({
        "name": K("nom de la zone"), "region": K("rectangle | circle | expression"),
        "x0": K("rectangle : x min"), "x1": K("rectangle : x max"),
        "y0": K("rectangle : y min"), "y1": K("rectangle : y max"),
        "center": K("circle : centre"), "radius": K("circle : rayon"),
        "expression": K("condition en x, y"), "darcy": K("coefficient de Darcy (1/m²)"),
        "permeability": K("perméabilité K = 1/darcy (m²)"),
        "forchheimer": K("coefficient de Forchheimer (1/m)"),
        "angle": K("orientation des axes principaux (°)")}, "zones poreuses", INC,
        many=True),
    "actuator_disk": Table({
        "name": K("nom"), "x0": K("début (x)"), "x1": K("fin (x)"), "radius": K("rayon"),
        "hub_radius": K("rayon du moyeu"), "thrust": K("poussée"),
        "thrust_coefficient": K("coefficient de poussée"), "mode": K("turbine | propeller"),
        "torque": K("couple"), "distribution": K("uniform | optimal"),
        "center": K("centre")}, "disques actuateurs", INC, many=True),
    "boundary": Table({"*": Table({
        "type": K("type de condition (incompressible : wall, inlet, outlet, symmetry, "
                  "farfield, axis, pressure_inlet ; compressible : farfield, inlet, outlet, "
                  "supersonic_inlet, supersonic_outlet, slip_wall, symmetry, wall)"),
        "U": K("vitesse [ux, uy] (valeurs ou formules en x, y)"),
        "p": K("pression"), "p0": K("pression totale"), "T": K("température"),
        "q": K("flux de chaleur pariétal", INC),
        "flow_rate": K("débit (au lieu de U)", INC),
        "profile": K("profil du débit : uniform | parabolic", INC),
        "scalars": Table({}, "valeurs des scalaires", INC, free=True),
        "scalar_flux": Table({}, "flux des scalaires", INC, free=True),
        "U_theta": K("vitesse de rotation (swirl)", INC),
        "omega": K("paroi tournante (rad/s, swirl)", INC),
        "mach": K("Mach", COMP), "pressure": K("= p", COMP), "temperature": K("= T", COMP),
        "density": K("masse volumique", COMP), "velocity": K("vitesse [u, v]", COMP),
        "angle": K("direction de l'écoulement (°)", COMP),
        "direction": K("direction [dx, dy]", COMP),
        "total_pressure": K("= p0", COMP), "T0": K("température totale", COMP),
        "total_temperature": K("= T0", COMP),
        "vortex": K("correction de tourbillon : point [x, y]", COMP)})},
        "conditions aux limites [boundary.<frontière>]"),
    "solver": Table({}, "réglages numériques"),         # complété plus bas (Settings)
    "output": Table({
        "directory": K("dossier des résultats"),
        "probes": K("sondes [[x, y], …]"),
        "lines": Table({"name": K("nom"), "start": K("[x, y]"), "end": K("[x, y]"),
                        "n": K("nombre de points")}, "profils sur des segments",
                       many=True),
        "forces": K("frontières où calculer les efforts (défaut : parois)"),
        "moment_center": K("centre des moments (défaut [0, 0])"),
        "checkpoint": K("écrit checkpoint.npz (défaut true)"),
        "checkpoint_minutes": K("sauvegarde périodique (min, défaut 5)"),
        "vtk": K("écrit fields.vtk (défaut true)"), "plots": K("figures (défaut true)"),
        "average_from": K("instationnaire : moyennes à partir de t", INC),
        "animate": K("instationnaire : grandeur animée (vorticity, U_mag, p…)", INC),
        "animate_every": K("une image tous les N pas", INC),
        "animate_fps": K("images par seconde (défaut 15)", INC),
        "vtk_every": K("instationnaire : fields_N.vtk tous les N pas", INC),
        "nusselt": K("bulk : Nusselt local par température de mélange", INC)},
        "sorties"),
    "sweep": Table({"parameter": K("clé balayée (ex. physics.reynolds)"),
                    "values": K("liste de valeurs"), "range": K("[début, fin, pas]"),
                    "continuation": K("part du point précédent (défaut true)"),
                    "jobs": K("calculs en parallèle (défaut 1)")}, "balayage"),
})

_SOLVER_CASE = {
    "mode": K("steady (défaut) | transient"),
    "dt": K("pas de temps (instationnaire)"), "t_end": K("temps final (instationnaire)"),
    "log_every": K("affichage tous les N itérations / pas (défaut 100)"),
    "fmg_levels": K("démarrage multigrille : niveaux grossiers (défaut 0)", INC),
    "fmg_tol": K("démarrage multigrille : tolérance des niveaux grossiers", INC),
}
_SOLVER_DOC = {
    # incompressible (fv2d/solver.py, Settings)
    "algorithm": "SIMPLE | SIMPLEC (défaut) | coupled",
    "relax_U": "sous-relaxation de la vitesse (défaut auto)",
    "relax_p": "sous-relaxation de la pression (défaut 1)",
    "relax_turb": "sous-relaxation de la turbulence (défaut 0.8)",
    "relax_T": "sous-relaxation de la température (défaut 0.9)",
    "relax_scalar": "sous-relaxation des scalaires (défaut 1)",
    "pseudo_dt": "pas de pseudo-temps (stationnaire pseudo-transitoire)",
    "pseudo_cfl": "CFL de pseudo-temps (stationnaire pseudo-transitoire)",
    "convection_U": "upwind | linearUpwind (défaut)",
    "convection_turb": "upwind (défaut) | linearUpwind",
    "convection_T": "upwind | linearUpwind (défaut)",
    "solver_p": "auto | direct | amg | bicgstab | pyamg",
    "solver_U": "auto | direct | amg | bicgstab | pyamg",
    "solver_turb": "auto | direct | amg | bicgstab | pyamg",
    "time_scheme": "auto | euler | backward | crankNicolson | rk1 | rk2 | rk3 | rk4 | ab2",
    "adjust_dt": "pas de temps adaptatif (défaut false)",
    "max_co": "Courant visé si adjust_dt (défaut 1)",
    "max_dt": "pas de temps maximal si adjust_dt",
    "cn_theta": "Crank-Nicolson : θ (défaut 0.5)",
    "ddt_phi_coeff": "correction ddtCorr (défaut 1)",
    "wall_treatment": "resolved (défaut) | wall_function",
    "backend": "cpu (défaut) | cuda | rocm | intel",
    "numba": "noyaux Numba (défaut false)", "threads": "fils de calcul Numba (défaut 1)",
    "n_outer": "boucles externes PIMPLE (défaut 2)",
    "n_corr": "corrections de pression (défaut 2)",
    "turbulence_every_outer": "turbulence à chaque boucle externe (défaut false)",
    "n_nonorth": "corrections non orthogonales (défaut 1)",
    "nonorth_limit": "limiteur de la correction non orthogonale (défaut 0.5)",
    # communs
    "max_iter": "itérations maximales (stationnaire)",
    "tol": "critère de convergence sur les résidus",
    "monitor_tol": "arrêt quand les efforts varient de moins de monitor_tol",
    "monitor_window": "fenêtre (itérations) de monitor_tol",
    # compressible (fv2d/compressible.py, CompressibleSettings)
    "flux": "roe (défaut) | hllc", "order": "1 | 2 (défaut)",
    "limiter": "venkatakrishnan (défaut) | barth_jespersen | none",
    "venkat_k": "seuil du limiteur de Venkatakrishnan (défaut 0.05)",
    "limiter_freeze": "limiteur gelé après N itérations (défaut 0 : jamais)",
    "entropy_fix": "correction d'entropie de Harten (défaut 0.1)",
    "cfl": "nombre CFL", "steady_scheme": "rk3 (défaut) | rk5 | implicit",
    "cfl_max": "implicite : CFL maximal", "cfl_growth": "implicite : croissance du CFL",
    "cfl_adapt": "implicite : réduction du CFL si le résidu croît",
    "cfl_cuts": "implicite : divisions du plafond de CFL (défaut 3)",
    "linear_sweeps": "implicite : balayages SGS", "implicit_jacobian": "roe | rusanov",
    "linear_solver": "gmres (défaut) | sgs", "linear_iter": "gmres : itérations",
    "linear_tol": "gmres : tolérance", "first_order_iter": "itérations d'ordre 1 au départ",
    "viscous_factor": "coefficient du pas de temps visqueux",
}


def _solver_table() -> Table:
    from dataclasses import fields as dc_fields

    from .compressible_case import _SETTINGS
    from .solver import Settings
    inc = {f.name for f in dc_fields(Settings)} | {"log_every"}
    keys = dict(_SOLVER_CASE)
    for name in sorted(inc | _SETTINGS):
        if name in keys:
            continue
        only = None if (name in inc and name in _SETTINGS) else (INC if name in inc else COMP)
        keys[name] = K(_SOLVER_DOC.get(name, ""), only)
    return Table(keys, "réglages numériques")


# ------------------------------------------------------------------ vérification
# noms de sections en français (fréquents chez les débutants)
_FRENCH = {"maillage": "mesh", "domaine": "domain", "corps": "bodies", "physique": "physics",
           "ecoulement": "flow", "écoulement": "flow", "initiale": "initial",
           "conditions": "boundary", "frontieres": "boundary", "frontières": "boundary",
           "limites": "boundary", "solveur": "solver", "sortie": "output",
           "sorties": "output", "energie": "energy", "énergie": "energy",
           "thermique": "energy", "scalaires": "scalars", "poreux": "porous",
           "balayage": "sweep", "turbulence_amont": "turbulence",
           # clés
           "vitesse": "U", "pression": "p", "viscosite": "nu", "viscosité": "nu",
           "debit": "flow_rate", "débit": "flow_rate", "rayon": "radius",
           "centre": "center", "nom": "name", "modele": "model", "modèle": "model",
           "iterations": "max_iter", "itérations": "max_iter", "tolerance": "tol",
           "tolérance": "tol", "pas_de_temps": "dt", "temps_final": "t_end",
           "sondes": "probes", "dossier": "directory"}
# conseils pour une clé lue par l'autre solveur
_HINTS = {("physics", "nu"): "viscosité du compressible : [flow] mu ou reynolds",
          ("physics", "reynolds"): "viscosité du compressible : [flow] mu ou reynolds",
          ("physics", "model"): "le solveur compressible est laminaire (ou Euler) : pas de "
                                "modèle de turbulence",
          ("boundary", "pressure"): "incompressible : p",
          ("boundary", "temperature"): "incompressible : T",
          ("boundary", "total_pressure"): "incompressible : p0"}


def _suggest(key: str, known, section=False) -> str:
    known = [k for k in known if k != "*"]
    low = {k.lower(): k for k in known}
    if key.lower() in low:
        best = low[key.lower()]
    elif key.lower() in _FRENCH and _FRENCH[key.lower()] in known:
        best = _FRENCH[key.lower()]
    else:
        m = difflib.get_close_matches(key, known, n=1, cutoff=0.6)
        best = m[0] if m else None
    if not best:
        return ""
    return f" — vouliez-vous dire [{best}] ?" if section else f" — vouliez-vous dire « {best} » ?"


def _where(parts) -> str:
    s, after = "", ""
    for p in parts:
        if isinstance(p, int):
            s, after = f"[[{s.strip('[]')}]] n° {p + 1}", ""
        elif s.startswith("[["):
            after += f", {p}"
        else:
            s = f"{s}.{p}" if s else p
    return (s if s.startswith("[[") else f"[{s}]") + after


class _Check:
    def __init__(self, kind: str, mesh_type: str):
        self.kind, self.mesh_type = kind, mesh_type
        self.warnings: list[str] = []
        self.errors: list[str] = []

    def table(self, value, spec: Table, parts):
        if spec.free:
            return
        if spec.many:
            if isinstance(value, dict):
                name = ".".join(str(p) for p in parts)
                self.errors.append(f"[{name}] : liste de tables attendue — écrire "
                                   f"[[{name}]] (double crochet) pour chaque élément.")
                return
            if not isinstance(value, list):
                self.errors.append(f"{_where(parts)} : liste de tables attendue.")
                return
            for i, item in enumerate(value):
                if not isinstance(item, dict):
                    self.errors.append(f"{_where(parts + [i])} : table attendue "
                                       f"(reçu {item!r}).")
                else:
                    self.keys(item, spec, parts + [i])
            return
        if not isinstance(value, dict):
            self.errors.append(f"{_where(parts)} : table attendue (reçu {value!r}).")
            return
        self.keys(value, spec, parts)

    def keys(self, value: dict, spec: Table, parts):
        for k, v in value.items():
            sub = spec.keys.get(k, spec.keys.get("*"))
            here = parts + [k]
            if sub is None:
                if not parts and not isinstance(v, (dict, list)):
                    home = [s for s, t in spec.keys.items()
                            if isinstance(t, Table) and k in t.keys]
                    self.warnings.append(f"{k} = {v!r} : clé hors de toute section, ignorée"
                                         + (f" — à placer sous [{home[0]}] ?" if home
                                            else _suggest(k, spec.keys)))
                elif not parts:
                    self.warnings.append(f"Section [{k}] inconnue, ignorée"
                                         f"{_suggest(k, spec.keys, section=True)}")
                else:
                    self.warnings.append(f"{_where(parts)} {k} : clé inconnue, ignorée"
                                         f"{_suggest(k, spec.keys)}")
                continue
            what = f"Section [{k}]" if not parts else f"{_where(parts)} {k}"
            if sub.only and sub.only != self.kind:
                hint = _HINTS.get((str(parts[0]) if parts else "", k))
                self.warnings.append(f"{what} : sans effet avec le solveur {self.kind} "
                                     f"(lue seulement par le solveur {sub.only})"
                                     + (f" ; {hint}." if hint else "."))
                continue
            if (sub.types and parts == ["mesh"] and self.mesh_type
                    and self.mesh_type not in sub.types):
                self.warnings.append(f"{what} : sans effet pour un maillage "
                                     f"« {self.mesh_type} » (utilisée par : "
                                     f"{', '.join(sub.types)}).")
                continue
            if isinstance(sub, Table):
                self.table(v, sub, here)


def _mesh_type(m) -> str:
    """Type de maillage, ou "" si inconnu / préréglage (pas de vérification par type)."""
    from ..mesh2d.builder import MESH_TYPES
    if not isinstance(m, dict) or "preset" in m:
        return ""
    t = str(m.get("type", "unstructured")).lower()
    return t if t in MESH_TYPES else ""


def check_case(cfg: dict, mesh_only: bool = False) -> list[str]:
    """Avertissements (liste de textes) ; ValueError si la structure est impossible.
    mesh_only : fichier de maillage seul (microrans mesh) : [mesh], [domain], [[bodies]]."""
    if not isinstance(cfg, dict):
        raise ValueError("Le fichier de cas doit contenir des sections [mesh], [physics]…")
    kind = COMP if (isinstance(cfg.get("physics"), dict)
                    and cfg["physics"].get("compressible")) else INC
    m = cfg.get("mesh", {})
    mesh_type = _mesh_type(m)
    schema = Table({**SCHEMA.keys, "solver": _solver_table()})
    if mesh_only:                                # maillage seul : autres sections ignorées
        schema = Table({k: v for k, v in schema.keys.items()
                        if k in ("mesh", "domain", "bodies")})
        if "mesh" in cfg:
            cfg = {k: v for k, v in cfg.items() if k in schema.keys}
        else:                                    # fichier de maillage : clés à la racine
            cfg = {"mesh": {k: v for k, v in cfg.items()
                            if k not in ("domain", "bodies")},
                   **{k: v for k, v in cfg.items() if k in ("domain", "bodies")}}
            mesh_type = _mesh_type(cfg["mesh"])
    c = _Check(kind, mesh_type)
    for sec in ("domain", "bodies"):              # sans effet pour ce type de maillage
        spec = schema.keys.get(sec)
        if sec in cfg and spec is not None and mesh_type and mesh_type not in spec.types:
            c.warnings.append(f"Section [{sec}] : sans effet pour un maillage "
                              f"« {mesh_type} » (utilisée par : {', '.join(spec.types)}).")
    c.keys({k: v for k, v in cfg.items()
            if not (k in ("domain", "bodies") and mesh_type
                    and mesh_type not in schema.keys[k].types)}, schema, [])
    if c.errors:
        raise ValueError("\n".join(c.errors))
    return c.warnings


def warn_case(cfg: dict, mesh_only: bool = False) -> list[str]:
    """check_case, avertissements émis (warnings.warn, catégorie CaseWarning)."""
    out = check_case(cfg, mesh_only)
    for w in out:
        warnings.warn(w, CaseWarning, stacklevel=2)
    return out
