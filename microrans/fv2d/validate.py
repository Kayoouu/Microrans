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
RECT, BLOCKS, OGRID, TRI, HYB, FILE, BOX = ("rectangle", "blocks", "ogrid", "unstructured",
                                            "hybrid", "file", "box")
_NAMES = Table({"left": K("nom du côté x = x0"), "right": K("nom du côté x = x1"),
                "bottom": K("nom du côté y = y0"), "top": K("nom du côté y = y1")},
               "noms des frontières du rectangle")
_NAMES3 = Table({**_NAMES.keys, "back": K("box : nom du côté z = z0"),
                 "front": K("box : nom du côté z = z1")},
                "noms des frontières du rectangle (du pavé : + back, front)", types=(RECT, BOX))
_RB = (RECT, BOX)
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
    "trailing_edge": K("naca : closed (défaut, −0.1036, workshops High-Order CFD) | open "
                       "(équation d'origine, bord de fuite épais) | sharp (prolongée jusqu'à "
                       "épaisseur nulle, Vassberg & Jameson 2010) ; prioritaire sur closed_te"),
    "path": K("file : fichier de contour (.dat Selig/Lednicer, .csv « x,y » ou « x;y » à "
              "virgule décimale, .dxf)"),
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
        "type": K("rectangle | blocks | ogrid | unstructured | hybrid | file | box "
                  "(défaut unstructured ; box : pavé 3D)"),
        "preset": K("maillage prédéfini (cavity, channel, cylinder-ogrid…) : remplace type"),
        "cut_axis": K("garde la moitié y > 0 (axisymétrique autour d'un corps)"),
        "x0": K("x min", types=_RB), "x1": K("x max", types=_RB),
        "y0": K("y min", types=_RB), "y1": K("y max", types=_RB),
        "z0": K("z min", types=(BOX,)), "z1": K("z max", types=(BOX,)),
        "nx": K("nombre de mailles selon x", types=_RB),
        "ny": K("nombre de mailles selon y", types=_RB),
        "nz": K("nombre de mailles selon z", types=(BOX,)),
        "grading": K("resserrement [gx, gy] (box : [gx, gy, gz] ; rapport dernière / "
                     "première maille)", types=_RB),
        "names": _NAMES3,
        "patch_types": K("types des frontières {nom = \"wall\" | \"patch\" | …}",
                         types=(RECT, FILE, BOX)),
        "periodic": K("paires de frontières périodiques [[\"a\", \"b\"], …]",
                      types=(RECT, BLOCKS, BOX)),
        "extrude": Table({
            "z0": K("z de départ (défaut 0)"), "z1": K("z d'arrivée (défaut 1)"),
            "nz": K("nombre de couches (défaut 1)"),
            "grading": K("resserrement selon z (rapport dernière / première couche)"),
            "names": Table({"back": K("nom de la face z = z0 (défaut back)"),
                            "front": K("nom de la face z = z1 (défaut front)")},
                           "noms des faces d'extrémité"),
            "patch_types": K("types des frontières {nom = \"wall\" | \"symmetry\" | …} "
                             "(défaut : type 2D ; patch pour back et front)"),
            "periodic": K("paires périodiques supplémentaires [[\"back\", \"front\"]]")},
            "3D : extrusion du maillage 2D selon z (quadrilatères → hexaèdres, triangles "
            "→ prismes)"),
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
        "reference_velocity": K("vitesse de référence U_ref : ν = U L / Re, coefficients "
                                "(défaut : vitesse d'entrée, sinon de paroi mobile, sinon 1)",
                                INC),
        "reference_length": K("longueur de référence L_ref (défaut 1)"),
        "reference_area": K("axisymétrique : aire de référence de Cd (défaut π L²/4)", INC),
        "model": K("laminar | sa | ke | kw | sst | sst_gamma (défaut laminar)", INC),
        "model_options": Table({}, "options du modèle (sa : ft2 ; sst_gamma : "
                               "kato_launder)", INC, free=True),
        "body_force": K("force volumique [fx, fy] (3D : [fx, fy, fz] ; conduite périodique)", INC),
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
        "U": K("vitesse initiale [ux, uy] (3D : [ux, uy, uz] ; valeurs ou formules en x, y, z)"),
        "perturbation": K("amplitude d'un tourbillon initial (déclenche le lâcher)", INC),
        "perturbation_center": K("centre du tourbillon initial (défaut [1.5, 0] ; 3D : [x, y] "
                                 "= tube selon z, ou [x, y, z])", INC),
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
        "gravity": K("gravité [gx, gy] (3D : [gx, gy, gz])"), "T0": K("température initiale"),
        "delta_T": K("écart de température de référence (nombre de Nusselt)"),
        "source": K("source de chaleur (valeur ou formule)")}, "thermique", INC),
    "scalars": Table({"*": Table({
        "diffusivity": K("diffusivité D"), "schmidt": K("Sc = ν / D (au lieu de "
                                                         "diffusivity)"),
        "Sc_t": K("Schmidt turbulent (défaut 0.7)"), "source": K("source (formule)"),
        "initial": K("valeur initiale"),
        "scheme": K("upwind | linearUpwind | linearUpwindLimited (défaut)")})},
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
        "U": K("vitesse [ux, uy] (3D : [ux, uy, uz] ; valeurs ou formules en x, y, z)"),
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
        "probes": K("sondes [[x, y], …] (3D : [[x, y, z], …])"),
        "lines": Table({"name": K("nom"), "start": K("[x, y] (3D : [x, y, z])"),
                        "end": K("[x, y] (3D : [x, y, z])"),
                        "n": K("nombre de points")}, "profils sur des segments",
                       many=True),
        "forces": K("frontières où calculer les efforts (défaut : parois)"),
        "moment_center": K("centre des moments (défaut [0, 0] ; 3D : [x, y, z], moment autour de "
                           "l'axe z)"),
        "checkpoint": K("écrit checkpoint.npz (défaut true)"),
        "checkpoint_minutes": K("sauvegarde périodique (min, défaut 5)"),
        "vtk": K("écrit fields.vtk (défaut true)"),
        "vtk_format": K("binary (défaut : binaire, valeurs exactes, 30 à 50 fois plus rapide "
                        "à écrire, ~30 % plus petit) | ascii (texte, 10 chiffres)"),
        "plots": K("figures (défaut true)"),
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
    "convection_U": "upwind | linearUpwind (défaut) | linearUpwindLimited (gradient limité, "
                    "plus robuste sur maillage déformé)",
    "convection_turb": "upwind (défaut) | linearUpwind | linearUpwindLimited",
    "convection_T": "upwind | linearUpwind (défaut) | linearUpwindLimited",
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
    "venkat_k": "seuil du limiteur de Venkatakrishnan (défaut 0.05 ; profils transsoniques : "
                "0.3, le défaut limite aussi hors des chocs, voir docs/compressible.md § 3.5)",
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


def _closest(key: str, known):
    known = [k for k in known if k != "*"]
    low = {k.lower(): k for k in known}
    if key.lower() in low:
        return low[key.lower()]
    if key.lower() in _FRENCH and _FRENCH[key.lower()] in known:
        return _FRENCH[key.lower()]
    m = difflib.get_close_matches(key, known, n=1, cutoff=0.6)
    return m[0] if m else None


def _suggest(key: str, known, section=False) -> str:
    best = _closest(key, known)
    if not best:
        return ""
    return f" — vouliez-vous dire [{best}] ?" if section else f" — vouliez-vous dire « {best} » ?"


def missing_bc_message(missing, types, extra_hint: str = "") -> str:
    """Conditions aux limites manquantes (M17) : noms en clair, types possibles, et pour
    les faces d'extrusion (back, front) la condition d'un écoulement 2D."""
    ext = [n for n in missing if n in ("back", "front")]
    tip = (" ; faces d'extrusion " + " et ".join(ext) + " : type = \"symmetry\" (écoulement "
           "2D extrudé) ou [mesh.extrude] periodic = [[\"back\", \"front\"]]") if ext else ""
    return (f"Conditions aux limites manquantes pour : {', '.join(missing)} — ajouter "
            + " et ".join(f"[boundary.{n}]" for n in missing[:3])
            + ("…" if len(missing) > 3 else "")
            + f" avec type = un de : {', '.join(types)}{tip}{extra_hint}.")


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
    def __init__(self, kind: str, mesh_type: str, root: Table | None = None):
        self.kind, self.mesh_type, self.root = kind, mesh_type, root
        self.warnings: list[str] = []
        self.errors: list[str] = []

    def _home(self, key: str, here: str) -> list[str]:
        """Sections (premier niveau) qui connaissent `key` (M16 : clé mal placée)."""
        if self.root is None:
            return []
        return [s for s, t in self.root.keys.items()
                if s != here and isinstance(t, Table) and not t.free and key in t.keys]

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
                    home = self._home(k, str(parts[0]))
                    self.warnings.append(f"{_where(parts)} {k} : clé inconnue, ignorée"
                                         + (" — clé de " + " ou de ".join(
                                             f"[{h}]" for h in home) + " : la déplacer"
                                            if home and len(parts) == 1
                                            else _suggest(k, spec.keys)))
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


# ------------------------------------------------------------------ valeurs
def _isnum(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


class _Values:
    """Valeurs impossibles : messages groupés (tout est signalé en une fois)."""

    def __init__(self, errors: list):
        self.errors = errors

    def num(self, sec: str, d: dict, key: str, *, gt=None, ge=None, le=None, integer=False,
            required=False, why=""):
        if not isinstance(d, dict) or d.get(key) is None:
            if required:
                self.errors.append(f"[{sec}] {key} manquant{why}.")
            return None
        v = d[key]
        if integer and not (isinstance(v, int) and not isinstance(v, bool)):
            whole = _isnum(v) and float(v).is_integer()
            self.errors.append(f"[{sec}] {key} = {v!r} : nombre entier attendu"
                               + (f" (écrire {int(v)}, sans « .0 »)." if whole else "."))
            return None
        if not _isnum(v):
            self.errors.append(f"[{sec}] {key} = {v!r} : nombre attendu.")
            return None
        if gt is not None and not v > gt:
            self.errors.append(f"[{sec}] {key} = {v!r} : doit être > {gt:g}.")
        elif ge is not None and not v >= ge:
            self.errors.append(f"[{sec}] {key} = {v!r} : doit être ≥ {ge:g}.")
        elif le is not None and not v <= le:
            self.errors.append(f"[{sec}] {key} = {v!r} : doit être ≤ {le:g}.")
        return float(v)

    def boolean(self, sec: str, d: dict, key: str):
        if isinstance(d, dict) and d.get(key) is not None and not isinstance(d[key], bool):
            self.errors.append(f"[{sec}] {key} = {d[key]!r} : true ou false attendu (sans "
                               "guillemets).")

    def choice(self, sec: str, d: dict, key: str, choices, *, lower=True):
        if not isinstance(d, dict) or d.get(key) is None:
            return None
        v = d[key]
        norm = str(v).lower() if lower else str(v)
        ok = [c.lower() for c in choices] if lower else list(choices)
        if norm not in ok:
            best = _closest(str(v), list(choices))
            self.errors.append(f"[{sec}] {key} = {v!r} inconnu. Choix : {', '.join(choices)}"
                               + (f" — vouliez-vous dire « {best} » ?" if best else "."))
        return norm

    def vec2(self, sec: str, d: dict, key: str, dim: int = 2):
        """[a, b] (3D : [a, b, c]) : nombres ou formules en x, y (, z)."""
        if not isinstance(d, dict) or d.get(key) is None:
            return
        v = d[key]
        if not isinstance(v, (list, tuple)) or len(v) != dim or not all(
                _isnum(c) or isinstance(c, str) for c in v):
            self.errors.append(
                f"[{sec}] {key} = {v!r} : 2 composantes attendues [x, y] (nombres ou "
                "formules en x, y)." if dim == 2 else
                f"[{sec}] {key} = {v!r} : 3 composantes attendues [x, y, z] (maillage 3D ; "
                "nombres ou formules en x, y, z).")


def case_dim(cfg: dict) -> int:
    """Dimension du cas lue dans [mesh] : 3 pour type = "box" ou [mesh.extrude], 2 sinon."""
    m = cfg.get("mesh") if isinstance(cfg, dict) else None
    if not isinstance(m, dict):
        return 2
    return 3 if (str(m.get("type", "")).lower() == BOX or m.get("extrude")) else 2


# options 2D seulement : (section, clé, valeur qui l'active (None : présence), message)
_ONLY_2D = [
    ("physics", "axisymmetric", True, "axisymétrique"),
    ("physics", "swirl", True, "rotation propre (swirl)"),
    ("porous", None, None, "zones poreuses [[porous]]"),
    ("actuator_disk", None, None, "disques actuateurs [[actuator_disk]]"),
    ("solver", "algorithm", "coupled", "solveur couplé (algorithm = \"coupled\")"),
    ("output", "animate", None, "animation (animate)"),
    ("physics", "compressible", True, "solveur compressible"),
]


def _check_3d(cfg: dict, kind: str, errors: list):
    """Options disponibles en 2D seulement (refusées avant de mailler)."""
    used = []
    for sec, key, val, what in _ONLY_2D:
        d = cfg.get(sec)
        if d is None:
            continue
        if key is None:
            used.append(what)
        elif isinstance(d, dict) and d.get(key) is not None and (
                val is None or str(d[key]).lower() == str(val).lower()):
            used.append(what)
    if used:
        errors.append("Maillage 3D ([mesh] type = \"box\" ou [mesh.extrude]) : "
                      f"{', '.join(used)} disponible(s) en 2D seulement.")


_MODE = {"unsteady": "transient", "instationnaire": "transient", "transitoire": "transient",
         "urans": "transient", "stationnaire": "steady", "rans": "steady"}


def _check_values(cfg: dict, kind: str, mesh_type: str, mesh_only: bool, errors: list,
                  warns: list | None = None):
    V = _Values(errors)
    m = cfg.get("mesh")
    if m is None:
        errors.append("[mesh] manquante : décrire le maillage (type = \"rectangle\", …), voir "
                      "les exemples (microrans examples).")
    elif isinstance(m, dict):
        if mesh_type in (RECT, BOX):
            why = f" (maillage {mesh_type})"
            axes = "xy" if mesh_type == RECT else "xyz"
            for a in axes:
                lo, hi = (V.num("mesh", m, k, required=True, why=why)
                          for k in (f"{a}0", f"{a}1"))
                if lo is not None and hi is not None and hi <= lo:
                    errors.append(f"[mesh] {a}1 = {m[a + '1']} doit être > {a}0 = "
                                  f"{m[a + '0']}.")
            for a in axes:
                V.num("mesh", m, f"n{a}", ge=1, integer=True, required=True, why=why)
            g = m.get("grading")
            if isinstance(g, (list, tuple)) and len(g) == len(axes):
                for c in g:
                    if _isnum(c) and c <= 0:
                        errors.append(f"[mesh] grading = {g!r} : rapports > 0 attendus.")
            elif mesh_type == BOX and g is not None:
                errors.append(f"[mesh] grading = {g!r} : 3 rapports attendus [gx, gy, gz] "
                              "(maillage box).")
            names = m.get("names")
            if mesh_type == RECT and isinstance(names, dict):
                for k in ("back", "front"):
                    if k in names:
                        errors.append(f"[mesh.names] {k} : face en z, pour un maillage "
                                      "« box » seulement (rectangle : left, right, bottom, "
                                      "top).")
        elif mesh_type == OGRID:
            V.num("mesh", m, "n_around", ge=4, integer=True)
            V.num("mesh", m, "n_radial", ge=2, integer=True)
            V.num("mesh", m, "farfield_radius", gt=0)
            V.num("mesh", m, "first_height", gt=0)
        elif mesh_type in _TRI_HYB:
            for k in ("h_max", "h_surface", "growth"):
                V.num("mesh", m, k, gt=0)
            V.num("mesh", m, "max_iter", ge=1, integer=True)
            lay = m.get("layers")
            if mesh_type == HYB and isinstance(lay, dict):
                V.num("mesh.layers", lay, "n", ge=1, integer=True)
                V.num("mesh.layers", lay, "first_height", gt=0)
                V.num("mesh.layers", lay, "ratio", gt=0)
        elif mesh_type == FILE and not str(m.get("path") or "").strip():
            errors.append("[mesh] path manquant : fichier .msh (Gmsh) ou .su2 à importer.")
        elif mesh_type == BLOCKS:
            for k in ("vertices", "blocks"):
                if not m.get(k):
                    errors.append(f"[mesh] {k} manquant ou vide (maillage multi-blocs).")
        e = m.get("extrude")
        if e is not None:
            if not isinstance(e, dict):
                errors.append("[mesh.extrude] : table attendue (z0, z1, nz).")
            elif mesh_type == BOX:
                errors.append("[mesh.extrude] : sans objet pour un maillage « box » (déjà "
                              "en 3D).")
            else:
                V.num("mesh.extrude", e, "nz", ge=1, integer=True)
                V.num("mesh.extrude", e, "grading", gt=0)
                z0, z1 = e.get("z0", 0.0), e.get("z1", 1.0)
                if _isnum(z0) and _isnum(z1) and z1 <= z0:
                    errors.append(f"[mesh.extrude] z1 = {z1} doit être > z0 = {z0}.")
                for k in ("z0", "z1"):
                    V.num("mesh.extrude", e, k)
    if mesh_only:
        return
    dim = case_dim(cfg)
    if dim == 3:
        _check_3d(cfg, kind, errors)
    ph = cfg.get("physics") if isinstance(cfg.get("physics"), dict) else {}
    sc = cfg.get("solver") if isinstance(cfg.get("solver"), dict) else {}
    if kind == INC:
        V.num("physics", ph, "nu", gt=0)
        V.num("physics", ph, "reynolds", gt=0)
        visc = ph.get("viscosity")
        newtonian = not (isinstance(visc, dict)
                         and str(visc.get("model", "newtonian")).lower() != "newtonian")
        if "nu" in ph and "reynolds" in ph:
            errors.append(f"[physics] nu = {ph['nu']} et reynolds = {ph['reynolds']} donnés "
                          "ensemble : garder l'un des deux (nu serait utilisé, reynolds "
                          "ignoré).")
        elif "nu" not in ph and "reynolds" not in ph and newtonian:
            errors.append("[physics] : donner la viscosité nu (m²/s) ou le nombre de "
                          "Reynolds reynolds.")
        if not newtonian:                           # avant : refus au lancement seulement
            from .rheology import MODELS as LAWS, PARAMS
            law = V.choice("physics.viscosity", visc, "model", LAWS)
            need = PARAMS.get(law, ())
            missing = [k for k in need if visc.get(k) is None]
            if missing:
                errors.append(f"[physics.viscosity] loi {law} : paramètre(s) "
                              f"{', '.join(missing)} manquant(s), sans valeur par défaut "
                              f"(attendus : {', '.join(need)}).")
            for k in ("K", "n", "nu0", "lambda", "m"):
                if k in need:
                    V.num("physics.viscosity", visc, k, gt=0)
            for k in ("nu_inf", "tau_y", "nu_min", "nu_max", "a", "relax"):
                V.num("physics.viscosity", visc, k, ge=0)
        for k in ("reference_velocity", "reference_length", "reference_area"):
            V.num("physics", ph, k, gt=0)
        model = "laminar"
        if ph.get("model") is not None:
            from ..models import MODELS, canonical_name
            try:
                model = canonical_name(str(ph["model"]))
            except ValueError as exc:
                model = None
                errors.append(f"[physics] model : {exc}")
            else:
                import inspect
                opts = ph.get("model_options")
                if isinstance(opts, dict):
                    ok = [p for p in inspect.signature(MODELS[model].__init__).parameters
                          if p not in ("self", "grid", "nu")]
                    for k in opts:
                        if k not in ok:
                            errors.append(f"[physics.model_options] {k} : option inconnue "
                                          f"pour le modèle {model} (options : "
                                          f"{', '.join(ok) or 'aucune'}).")
        if not newtonian and model not in (None, "laminar"):   # avant : refus au lancement
            errors.append(f"[physics] model = \"{model}\" avec [physics.viscosity] (fluide non "
                          "newtonien) : non newtonien en laminaire uniquement (les modèles "
                          "de turbulence supposent un fluide newtonien).")
        if sc.get("wall_treatment") == "wall_function":   # avant : refus au lancement
            if model in ("ke", "sst_gamma"):
                errors.append(f"[solver] wall_treatment = \"wall_function\" : lois de paroi "
                              f"incompatibles avec le modèle {model} (" + (
                                  "k-ε bas-Reynolds" if model == "ke" else
                                  "la couche limite laminaire doit être résolue")
                              + ", y⁺ ≈ 1) : garder \"resolved\", ou SA, k-ω, SST.")
            elif model == "laminar" and warns is not None:
                warns.append("[solver] wall_treatment = \"wall_function\" : loi de paroi "
                             "turbulente, sans objet en laminaire (garder \"resolved\").")
        V.vec2("physics", ph, "body_force", dim)
    else:
        fl = cfg.get("flow") if isinstance(cfg.get("flow"), dict) else {}
        V.num("flow", fl, "mach", ge=0)
        for k in ("pressure", "temperature", "density", "gas_constant", "prandtl", "mu",
                  "reynolds"):
            V.num("flow", fl, k, gt=0)
        V.num("flow", fl, "gamma", gt=1)
        V.vec2("flow", fl, "velocity")
    V.num("physics", ph, "angle_of_attack")
    # [solver]
    mode = sc.get("mode", "steady")
    if str(mode).lower() not in ("steady", "transient"):
        alt = _MODE.get(str(mode).lower())
        errors.append(f"[solver] mode = {mode!r} inconnu : steady (stationnaire) ou "
                      f"transient (instationnaire)" + (f" — vouliez-vous dire « {alt} » ?"
                                                       if alt else "."))
    elif str(mode).lower() == "transient":
        V.num("solver", sc, "t_end", gt=0, required=True, why=" (mode transient : temps final)")
        V.num("solver", sc, "dt", gt=0, required=kind == INC,
              why=" (mode transient : pas de temps)")
    V.num("solver", sc, "max_iter", ge=1, integer=True)
    V.num("solver", sc, "tol", ge=0)
    V.num("solver", sc, "monitor_tol", gt=0)
    V.num("solver", sc, "monitor_window", ge=1, integer=True)
    oc = cfg.get("output") if isinstance(cfg.get("output"), dict) else {}
    vf = oc.get("vtk_format", "binary")
    if str(vf).lower() not in ("binary", "ascii"):
        errors.append(f"[output] vtk_format = {vf!r} inconnu : binary (défaut, valeurs exactes) "
                      "ou ascii (texte).")
    V.num("solver", sc, "log_every", ge=1, integer=True)
    if oc.get("probes") is not None:                # U16 : avant, refus au lancement
        from .sampling import parse_points
        try:
            parse_points(oc["probes"], dim)
        except (ValueError, TypeError) as exc:
            errors.append(f"[output] probes = {oc['probes']!r} : {exc}")
    if kind == INC:
        for k in ("relax_U", "relax_p", "relax_turb", "relax_T", "relax_scalar"):
            V.num("solver", sc, k, gt=0, le=1)
        for k in ("n_outer", "n_corr", "n_nonorth"):
            V.num("solver", sc, k, ge=1, integer=True)
        for k in ("max_co", "pseudo_cfl", "pseudo_dt", "max_dt"):
            V.num("solver", sc, k, gt=0)
        V.num("solver", sc, "fmg_levels", ge=0, integer=True)
        V.choice("solver", sc, "algorithm", ("SIMPLE", "SIMPLEC", "coupled"))
        V.choice("solver", sc, "wall_treatment", ("resolved", "wall_function"), lower=False)
        for k in ("convection_U", "convection_turb", "convection_T"):
            V.choice("solver", sc, k, ("upwind", "linearUpwind", "linearUpwindLimited"),
                     lower=False)
        from ..linalg import SOLVERS
        for k in ("solver_p", "solver_U", "solver_turb"):
            V.choice("solver", sc, k, tuple(SOLVERS), lower=False)
        from .solver import TIME_SCHEMES
        V.choice("solver", sc, "time_scheme", ("auto", *TIME_SCHEMES), lower=False)
        # clés « expert » (M18 : texte → erreur interne, valeurs hors bornes acceptées)
        th = sc.get("cn_theta")
        if _isnum(th) and not 0.5 <= th <= 1:
            errors.append(f"[solver] cn_theta = {th!r} : entre 0.5 (Crank-Nicolson, ordre 2) "
                          "et 1 (Euler implicite) ; θ < 0.5 est instable.")
        else:
            V.num("solver", sc, "cn_theta")
        for k in ("ddt_phi_coeff", "nonorth_limit"):
            V.num("solver", sc, k, ge=0, le=1)
        V.num("solver", sc, "threads", ge=0, integer=True)
        V.num("solver", sc, "fmg_tol", gt=0)
        for k in ("adjust_dt", "numba", "turbulence_every_outer"):
            V.boolean("solver", sc, k)
    else:
        for k in ("cfl", "cfl_max"):
            V.num("solver", sc, k, gt=0)
        V.num("solver", sc, "first_order_iter", ge=0, integer=True)
        # clés « expert » (M14, M18)
        from .compressible import FLUXES, LIMITERS, STEADY_SCHEMES
        V.choice("solver", sc, "flux", FLUXES)
        lim = str(sc.get("limiter", "")).lower().replace("-", "_")
        if lim not in ("venkat", "venkatakrishnan_wang", "barth", "bj"):
            V.choice("solver", sc, "limiter", LIMITERS)
        V.choice("solver", sc, "steady_scheme", STEADY_SCHEMES)
        V.choice("solver", sc, "implicit_jacobian", ("roe", "rusanov"))
        V.choice("solver", sc, "linear_solver", ("gmres", "sgs"))
        if sc.get("order") is not None and sc["order"] not in (1, 2):
            errors.append(f"[solver] order = {sc['order']!r} : 1 ou 2 attendu.")
        V.num("solver", sc, "venkat_k", gt=0)
        V.num("solver", sc, "entropy_fix", ge=0)
        V.num("solver", sc, "viscous_factor", gt=0)
        for k in ("cfl_growth", "cfl_adapt"):
            V.num("solver", sc, k, ge=1)
        V.num("solver", sc, "linear_tol", gt=0, le=1)
        for k in ("limiter_freeze", "cfl_cuts"):
            V.num("solver", sc, k, ge=0, integer=True)
        for k in ("linear_sweeps", "linear_iter"):
            V.num("solver", sc, k, ge=1, integer=True)
    # [boundary.*], [initial]
    bnd = cfg.get("boundary") if isinstance(cfg.get("boundary"), dict) else {}
    for name, spec in bnd.items():
        if not isinstance(spec, dict):
            continue
        sec = f"boundary.{name}"
        if not spec.get("type"):
            errors.append(f"[{sec}] type manquant (ex. type = \"wall\").")
            continue
        for k in ("U", "velocity", "direction"):
            V.vec2(sec, spec, k, dim)
        for k in (("q", "omega", "flow_rate") if kind == INC else ("angle",)):
            V.num(sec, spec, k)                     # pas de formule pour ces valeurs
        t = str(spec["type"]).lower()
        if kind == INC and t == "inlet" and "U" not in spec and "flow_rate" not in spec:
            errors.append(f"[{sec}] (inlet) : donner la vitesse U = [ux, uy] ou le débit "
                          "flow_rate.")
        if kind == INC and t == "farfield" and "U" not in spec:
            errors.append(f"[{sec}] (farfield) : donner la vitesse amont U = [ux, uy].")
    V.vec2("initial", cfg.get("initial") if isinstance(cfg.get("initial"), dict) else {}, "U",
           dim)
    en = cfg.get("energy")
    if kind == INC and isinstance(en, dict):
        V.num("energy", en, "Pr", gt=0)
        V.num("energy", en, "Pr_t", gt=0)
        V.vec2("energy", en, "gravity", dim)


# ordres de grandeur mesurés (cavité, SIMPLEC, solveur de pression AMG, un cœur) : mémoire
# ~250 Mo + ~1 Ko par cellule, ~15 µs par cellule et par itération
BIG_MESH = 500_000


def _size_warning(m: dict, mesh_type: str) -> str | None:
    """Avertissement pour un maillage structuré très gros (nombre de cellules exact avant
    maillage ; pas d'estimation fiable pour les triangles, raffinements compris)."""
    try:
        if mesh_type == RECT:
            n = int(m["nx"]) * int(m["ny"])
        elif mesh_type == OGRID:
            n = int(m.get("n_around", 128)) * int(m.get("n_radial", 64))
        elif mesh_type == BLOCKS:
            n = sum(int(b["cells"][0]) * int(b["cells"][1]) for b in m["blocks"])
        elif mesh_type == BOX:
            n = int(m["nx"]) * int(m["ny"]) * int(m["nz"])
        else:
            return None
        if isinstance(m.get("extrude"), dict):
            n *= int(m["extrude"].get("nz", 1))
    except (KeyError, TypeError, ValueError, IndexError):
        return None
    if n < BIG_MESH:
        return None
    count = f"{n:,}".replace(",", " ")
    if mesh_type == BOX or isinstance(m.get("extrude"), dict):
        # 3D mesuré (cavité cubique, 4 cœurs, machine virtuelle de charge variable) :
        # 10⁶ hexaèdres → pic 2.9 à 3.2 Go, 8.5 à 13.8 s par itération en laminaire (SST :
        # +40 %), préparation et écriture 78 s (maillage ~35 s + distance à la paroi exacte
        # 25 s, lot E2 ; avant : 3.4 à 5.3 min) ; exécutable et Python identiques à la
        # mesure près
        gb, lo, hi = 0.3 + n * 2.9e-6, n * 8.5e-6, n * 14e-6
        prep = f"{max(1, round(n * 6e-5 / 60))} à {max(2, round(n * 1.2e-4 / 60))}"
        return (f"[mesh] {count} cellules (3D) : prévoir ~{gb:.1f} Go de mémoire, ~{lo:.0f} à "
                f"{hi:.0f} s par itération en laminaire (turbulent : ~40 % de plus) et ~{prep} "
                "min de préparation (maillage et distance à la paroi) — ordre de grandeur "
                "mesuré, variable selon la machine ; plusieurs centaines d'itérations sont "
                "nécessaires. Régler d'abord le cas sur un maillage plus grossier.")
    gb = 0.25 + n * 1e-6
    it = n * 15e-6
    return (f"[mesh] {count} cellules : prévoir ~{gb:.1f} Go de mémoire et ~{it:.0f} s par "
            "itération (ordre de grandeur mesuré, variable selon la machine ; plusieurs "
            "centaines d'itérations sont nécessaires). Régler d'abord le cas sur un maillage "
            "plus grossier.")


def _check_bodies(cfg: dict, mesh_type: str, errors: list, warns: list):
    """Corps qui se recouvrent, trop proches ou au bord du domaine, détectés avant de mailler
    (avant : hybride refusé après 6 à 20 s « non manifold », ou maillage faux avec des
    couches dans un solide ; non structuré : corps disparu en silence)."""
    bodies = cfg.get("bodies")
    from ..mesh2d.geometry import NACA_TE
    for b in bodies if isinstance(bodies, list) else []:
        te = b.get("trailing_edge") if isinstance(b, dict) else None
        if te is not None and str(b.get("type", "")).lower() != "naca":
            warns.append(f"[[bodies]] « {b.get('name', b.get('type'))} » : trailing_edge "
                         f"sans effet (type = \"{b.get('type')}\" ; seulement pour type = "
                         "\"naca\").")                            # M19
        elif te is not None and str(te).lower() not in NACA_TE:
            errors.append(f"[[bodies]] « {b.get('name', b.get('type'))} » : trailing_edge = "
                          f"{te!r} inconnu (" + ", ".join(NACA_TE) + ").")
    if mesh_type == OGRID and (not isinstance(bodies, list) or len(bodies) != 1):
        errors.append(f"[[bodies]] : le maillage en O entoure exactement un corps "
                      f"({len(bodies) if isinstance(bodies, list) else 0} donné(s)). "
                      "Plusieurs corps : maillage « unstructured » ou « hybrid ».")
    if mesh_type not in (TRI, HYB) or not isinstance(bodies, list) or not bodies:
        return
    import numpy as np

    from ..mesh2d.builder import _outer
    from ..mesh2d.geometry import shape_from_dict
    shapes = []
    for i, b in enumerate(bodies):
        if not isinstance(b, dict) or b.get("type") == "file":
            continue                                # contour lu dans un fichier : non vérifié
        try:
            sh = shape_from_dict(dict(b))
            shapes.append((str(b.get("name", f"corps {i + 1}")), sh, sh.boundary_curve(n=400)))
        except Exception:                           # dimensions fausses : signalées ailleurs
            continue
    m = cfg.get("mesh", {})
    lay = m.get("layers", {}) if isinstance(m.get("layers"), dict) else {}
    try:
        n, h1, r = int(lay.get("n", 10)), float(lay.get("first_height", 1e-3)), \
            float(lay.get("ratio", 1.2))
        thick = h1 * n if abs(r - 1.0) < 1e-12 else h1 * (r ** n - 1.0) / (r - 1.0)
    except (TypeError, ValueError):
        thick = 0.0
    hyb = mesh_type == HYB
    for k, (na, a, pa) in enumerate(shapes):
        for nb, b, pb in shapes[k + 1:]:
            size = max(np.ptp(pa, axis=0).max(), np.ptp(pb, axis=0).max())
            gap = min(float(b.sdf(pa).min()), float(a.sdf(pb).min()))
            if gap < 1e-6 * size:
                txt = (f"[[bodies]] « {na} » et « {nb} » se recouvrent ou se touchent")
                if hyb:
                    errors.append(f"{txt} : impossible en maillage hybride (couches de paroi "
                                  "de l'un dans l'autre). Les écarter, ou maillage "
                                  "« unstructured » pour un obstacle composé.")
                else:
                    warns.append(f"{txt} : ils formeront un seul obstacle (frontières "
                                 "fusionnées ; une condition pour un corps entièrement "
                                 "recouvert restera sans frontière).")
            elif hyb and gap < 2.0 * thick:
                errors.append(f"[[bodies]] « {na} » et « {nb} » : écart {gap:.3g} < 2 × "
                              f"épaisseur des couches de paroi ({thick:.3g}) : couches en "
                              "collision. Écarter les corps ou réduire [mesh.layers] n, "
                              "first_height ou ratio.")
    try:
        outer = _outer(cfg)
    except Exception:
        return
    for na, a, pa in shapes:
        d = -outer.sdf(pa)                          # > 0 : dans le domaine
        if d.max() <= 0.0:
            txt = f"[[bodies]] « {na} » est entièrement hors du domaine de calcul [domain]"
            if hyb:
                errors.append(txt + ".")
            else:
                warns.append(txt + " : il sera ignoré.")
        elif hyb and d.min() < thick:
            errors.append(f"[[bodies]] « {na} » : à {max(d.min(), 0.0):.3g} du bord du "
                          f"domaine, moins que l'épaisseur des couches de paroi ({thick:.3g})"
                          " : impossible en maillage hybride. Agrandir [domain], déplacer le "
                          "corps, ou maillage « unstructured » (corps coupé par le bord "
                          "permis).")


def _mesh_type(m) -> str:
    """Type de maillage, ou "" si inconnu / préréglage (pas de vérification par type)."""
    from ..mesh2d.builder import MESH_TYPES
    if not isinstance(m, dict) or "preset" in m:
        return ""
    t = str(m.get("type", "unstructured")).lower()
    return t if t in MESH_TYPES else ""


def check_case(cfg: dict, mesh_only: bool = False, values: bool = True) -> list[str]:
    """Avertissements (liste de textes) ; ValueError si la structure ou une valeur est
    impossible. mesh_only : fichier de maillage seul (microrans mesh) : [mesh], [domain],
    [[bodies]] ; values = False : structure seulement (cas en cours d'écriture)."""
    if not isinstance(cfg, dict) or not cfg:
        raise ValueError("Fichier de cas vide : il faut au moins [mesh], [physics] et "
                         "[boundary.<frontière>] (partir d'un exemple : microrans examples).")
    if values and not mesh_only and "mesh" in cfg and not any(
            k in cfg for k in ("physics", "boundary", "flow")):
        raise ValueError("Ce fichier décrit seulement un maillage (ni [physics] ni "
                         "[boundary]) : le générer avec « microrans mesh <fichier> » ou le "
                         "bouton « Générer le maillage » de l'interface ; pour un calcul, "
                         "partir d'un exemple de calcul (microrans examples).")
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
    c = _Check(kind, mesh_type, schema)
    for sec in ("domain", "bodies"):              # sans effet pour ce type de maillage
        spec = schema.keys.get(sec)
        if sec in cfg and spec is not None and mesh_type and mesh_type not in spec.types:
            c.warnings.append(f"Section [{sec}] : sans effet pour un maillage "
                              f"« {mesh_type} » (utilisée par : {', '.join(spec.types)}).")
    c.keys({k: v for k, v in cfg.items()
            if not (k in ("domain", "bodies") and mesh_type
                    and mesh_type not in schema.keys[k].types)}, schema, [])
    if values and not c.errors:                   # structure lisible : valeurs
        _check_values(cfg, kind, mesh_type, mesh_only, c.errors, c.warnings)
        big = _size_warning(cfg.get("mesh"), mesh_type) if isinstance(cfg.get("mesh"),
                                                                      dict) else None
        if big:
            c.warnings.append(big)
        if not c.errors:
            _check_bodies(cfg, mesh_type, c.errors, c.warnings)
    if c.errors:                                  # une faute de frappe explique souvent l'erreur
        raise ValueError("\n".join(c.errors + [f"Remarque : {w}" for w in c.warnings]))
    return c.warnings


def warn_case(cfg: dict, mesh_only: bool = False) -> list[str]:
    """check_case, avertissements émis (warnings.warn, catégorie CaseWarning)."""
    out = check_case(cfg, mesh_only)
    for w in out:
        warnings.warn(w, CaseWarning, stacklevel=2)
    return out


# ------------------------------------------------------------------ référence (docs)
def _defaults() -> dict:
    """{clé de [solver] : (défaut incompressible, défaut compressible)} lus dans le code."""
    from dataclasses import MISSING, fields

    from .compressible import CompressibleSettings
    from .solver import Settings

    def get(cls):
        return {f.name: f.default for f in fields(cls) if f.default is not MISSING}
    inc, comp = get(Settings), get(CompressibleSettings)
    return {k: (inc.get(k, MISSING), comp.get(k, MISSING)) for k in set(inc) | set(comp)}


def _fmt_default(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "auto"
    if isinstance(v, float) and v == float("inf"):
        return "aucun"
    return f"`{v}`" if isinstance(v, str) else f"{v:g}" if isinstance(v, float) else str(v)


def _cell(text: str) -> str:
    return text.replace("|", "\\|")             # « a | b » dans un tableau Markdown


def _applies(item) -> str:
    out = []
    if item.only:
        out.append(f"{item.only} seulement")
    if item.types:
        out.append("maillage " + ", ".join(item.types))
    return " ; ".join(out)


def reference_markdown() -> str:
    """Référence de toutes les clés du fichier de cas, tirée de SCHEMA (docs/reference_cas.md ;
    un test vérifie que le fichier publié est à jour)."""
    from dataclasses import MISSING
    defaults = _defaults()
    lines = [
        "# Référence des clés du fichier de cas",
        "",
        "Document généré à partir de la liste des clés que le logiciel vérifie "
        "(`microrans/fv2d/validate.py`) : il contient exactement les clés reconnues. Ne pas "
        "le modifier à la main ; le régénérer avec",
        "`python -m microrans.fv2d.validate > docs/reference_cas.md`.",
        "",
        "Toute autre clé est signalée (« clé inconnue, ignorée — vouliez-vous dire … ? »). "
        "« incompressible seulement » / « compressible seulement » : clé lue par un seul "
        "des deux solveurs (`[physics] compressible = true`) ; « maillage … » : types de "
        "maillage qui l'utilisent. Les valeurs par défaut de `[solver]` sont lues dans le "
        "code (`Settings`, `CompressibleSettings`).",
        "",
        "Exemples complets : `microrans examples` ; tutoriel : `docs/tutoriel.md`.",
    ]

    def table(t: Table, path: str, level: int):
        rows, subs = [], []
        for k, v in t.keys.items():
            if isinstance(v, Table):
                subs.append((k, v))
            else:
                rows.append((k, v))
        if rows:
            lines.append("")
            lines.append("| Clé | Signification | S'applique à |")
            lines.append("|---|---|---|")
            for k, v in rows:
                lines.append(f"| `{k}` | {_cell(v.doc)} | {_applies(v)} |")
        for k, v in subs:
            name = f"{path}.<nom>" if k == "*" else f"{path}.{k}"
            head = f"[[{name}]]" if v.many else f"[{name}]"
            lines.append("")
            doc = v.doc or ("une section par nom (frontière, scalaire…)" if k == "*" else "")
            lines.append(f"{'#' * level} `{head}` — {doc}"
                         + (f" ({_applies(v)})" if _applies(v) else ""))
            if v.free:
                lines.append("")
                lines.append("Contenu libre, sous la forme `nom = valeur`.")
            table(v, name, min(level + 1, 6))

    for sec, t in SCHEMA.keys.items():
        lines.append("")
        head = f"[[{sec}]]" if t.many else f"[{sec}]"
        lines.append(f"## `{head}` — {t.doc}" + (f" ({_applies(t)})" if _applies(t) else ""))
        if sec == "solver":
            lines.append("")
            lines.append("| Clé | Signification | Défaut (incompressible) | "
                         "Défaut (compressible) | S'applique à |")
            lines.append("|---|---|---|---|---|")
            for k, v in sorted(_solver_table().keys.items()):
                di, dc = defaults.get(k, (MISSING, MISSING))
                fi = "" if di is MISSING else _fmt_default(di)
                fc = "" if dc is MISSING else _fmt_default(dc)
                lines.append(f"| `{k}` | {_cell(v.doc)} | {fi} | {fc} | {_applies(v)} |")
            continue
        table(t, sec, 3)
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stdout.write(reference_markdown())
