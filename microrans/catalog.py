"""Exemples fournis, classés : accueil de l'interface et « microrans examples ».

Durées : mesurées un exemple à la fois en ligne de commande (démarrage, maillage, calcul et
figures compris) sur la machine de test (4 cœurs) ; ordre de grandeur, variable selon la
machine. Polaire : les 10 incidences du balayage ([sweep]).
"""
from __future__ import annotations

from pathlib import Path

# (groupe, [(fichier sans .toml, titre, durée)])
GROUPS = [
    ("Commencer ici", [
        ("cavite_re100", "Cavité entraînée, Re = 100 (cas d'école, référence Ghia)", "6 s"),
        ("cylindre_re20", "Cylindre, Re = 20 : traînée Cd, maillage en O", "5 s"),
    ]),
    ("Laminaire", [
        ("sphere_re100_axisym", "Sphère, Re = 100, axisymétrique", "4 s"),
        ("cylindre_re100_urans", "Cylindre, Re = 100 : lâcher de tourbillons (instationnaire)",
         "2,5 min"),
    ]),
    ("Turbulence et aérodynamique", [
        ("plaque_plane_sa", "Plaque plane turbulente, Spalart-Allmaras", "10 s"),
        ("plaque_plane_loi_de_paroi", "Plaque plane avec lois de paroi (maillage grossier)",
         "5 s"),
        ("tuyau_turbulent_sst", "Tuyau turbulent établi, k-ω SST", "10 s"),
        ("naca0012_sa", "Profil NACA 0012 à 4°, Re = 10⁶", "30 s"),
        ("naca0012_polaire", "Polaire Cl(α), Cd(α) du NACA 0012 : « Lancer le balayage »",
         "4 min"),
        ("plaque_plane_transition_t3a", "Transition laminaire → turbulent (T3A)", "40 s"),
    ]),
    ("Thermique", [
        ("convection_naturelle_ra1e5", "Convection naturelle en cavité chauffée, Ra = 10⁵",
         "3 s"),
        ("tuyau_thermique", "Tuyau chauffé à flux uniforme, axisymétrique", "4 s"),
    ]),
    ("Fluides et modèles particuliers", [
        ("sang_carreau_artere", "Sang (non newtonien, Carreau) dans une artère", "4 s"),
        ("melange_deux_courants", "Mélange de deux courants (scalaires transportés)", "9 s"),
        ("filtre_poreux_conduite", "Filtre poreux dans une conduite", "4 s"),
        ("eolienne_disque_actuateur", "Éolienne : disque actuateur, axisymétrique", "5 s"),
    ]),
    ("Compressible", [
        ("compressible_tube_sod", "Tube à choc de Sod", "5 s"),
        ("compressible_rampe_mach2", "Rampe supersonique à Mach 2 : choc oblique", "15 s"),
        ("compressible_plaque_laminaire", "Plaque plane laminaire à Mach 0.2", "13 s"),
        ("compressible_naca0012_transsonique", "NACA 0012 transsonique, Mach 0.8", "50 s"),
    ]),
    ("3D — figures en coupe x, y ou z, champs complets dans ParaView", [
        ("conduite_carree_3d", "Conduite carrée laminaire (solution exacte)", "8 s"),
        ("canal_turbulent_3d", "Canal turbulent Re_τ = 395 extrudé, Spalart-Allmaras",
         "12 s"),
        ("cavite_cubique_re100_3d", "Cavité cubique Re = 100 (démonstration)", "25 s"),
    ]),
    ("Maillage seul (pas de calcul)", [
        ("mesh_cylindre_hybride", "Maillage hybride autour d'un cylindre", "25 s"),
        ("mesh_naca_multi", "Deux profils, dont un contour importé", "40 s"),
    ]),
]
OTHERS = "Autres"


def header(path: Path, commands: bool = True) -> str:
    """Commentaire d'en-tête du fichier (lignes « # » du début), sans les dièses ;
    commands = False : sans les lignes de commande (« microrans run2d … »)."""
    lines = []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        if not ln.startswith("#"):
            break
        ln = ln.lstrip("#").strip()
        if not commands and ("microrans run2d" in ln or "microrans mesh" in ln):
            continue
        lines.append(ln)
    return " ".join(x for x in lines if x)


def catalog(directory: Path) -> list[tuple[str, list[tuple[Path, str, str]]]]:
    """[(groupe, [(fichier, titre, durée)])] ; un exemple absent de GROUPS va dans
    « Autres » sous son nom de fichier (rien n'est caché)."""
    files = {f.stem: f for f in sorted(Path(directory).glob("*.toml"))}
    out = []
    for group, items in GROUPS:
        rows = [(files.pop(stem), title, dur) for stem, title, dur in items if stem in files]
        if rows:
            out.append((group, rows))
    if files:
        out.append((OTHERS, [(f, f.stem, "") for f in files.values()]))
    return out
