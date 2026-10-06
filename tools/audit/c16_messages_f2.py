"""Campagne 16 (lot F2) : messages des points M13 à M22 et L6, en ligne de commande.
Chaque cas : petit fichier de cas (ou commande) et dernière ligne utile de la sortie.

    python tools/audit/c16_messages_f2.py <dossier de sortie>
"""
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")

BASE = """[mesh]
type = "rectangle"
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
nx = {nx}
ny = {ny}
names = {{ left = "walls", right = "walls", bottom = "walls", top = "lid" }}
{mesh_extra}
[physics]
nu = 0.01
{physics_extra}
[boundary.lid]
type = "wall"
U = [1.0, 0.0]
{boundary_walls}
[solver]
max_iter = 3
{solver_extra}
"""
WALLS = '[boundary.walls]\ntype = "wall"'
COMP = """[mesh]
type = "rectangle"
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 0.5
nx = 8
ny = 4
names = { left = "inlet", right = "outlet", bottom = "wall", top = "top" }
[physics]
compressible = true
[flow]
mach = 0.5
pressure = 101325.0
temperature = 300.0
[boundary]
inlet = { type = "farfield" }
outlet = { type = "farfield" }
wall = { type = "slip_wall" }
top = { type = "farfield" }
[solver]
max_iter = 3
"""


def cas(nom, texte, args=(), encodage="utf-8", bom=False):
    f = OUT / f"{nom}.toml"
    data = texte.encode(encodage)
    f.write_bytes((b"\xef\xbb\xbf" if bom else b"") + data)
    return run(nom, ["run2d", str(f), "-o", str(OUT / f"res_{nom}"), "--no-plot", "-q",
                     *args])


def run(nom, args, cwd=None):
    p = subprocess.run([sys.executable, "-m", "microrans", *args], env=env, cwd=cwd or OUT,
                       capture_output=True, text=True)
    lignes = [x for x in (p.stderr + p.stdout).splitlines() if x.strip()]
    utiles = ([x for x in lignes if "Erreur" in x or "Vitesse moyenne" in x]
              + [x for x in lignes if "ATTENTION" in x]) or lignes[-2:]
    print(f"== {nom} (code {p.returncode})")
    for x in utiles[:6]:
        print("   " + x[:300])
    return p


def base(**kw):
    d = dict(nx=4, ny=4, mesh_extra="", physics_extra="", boundary_walls=WALLS,
             solver_extra="")
    d.update(kw)
    return BASE.format(**d)


cas("M13_periodique", base(mesh_extra='periodic = [["walls", "lidd"]]'))
cas("M14_M18_texte_compressible", COMP.replace("max_iter = 3",
                                               'max_iter = 3\nvenkat_k = "0.3"\n'
                                               'linear_iter = "x"\nlinear_tol = "y"'))
cas("M18_negatifs_compressible", COMP.replace(
    "max_iter = 3", "max_iter = 3\nvenkat_k = -1\nlimiter_freeze = -1\nentropy_fix = -1\n"
    "cfl_growth = -1\ncfl_cuts = -1\nlinear_iter = -1\nviscous_factor = -1"))
cas("M18_incompressible", base(solver_extra='cn_theta = -1\nddt_phi_coeff = -1\n'
                                            'nonorth_limit = -1\nthreads = "4"'))
cas("M15_type_patch", base(mesh_extra='patch_types = { walls = "symetrie" }'))
cas("M16_mauvaise_section", base(physics_extra="moment_center = [0.25, 0.0]"))
cas("M17_cl_manquante", base(boundary_walls=""))
cas("M21_une_maille", base(nx=1, ny=1))
cas("M21_une_maille_3d", base(nx=1, ny=1, mesh_extra="extrude = { z0 = 0, z1 = 1, nz = 1 }")
    .replace("U = [1.0, 0.0]", "U = [1.0, 0.0, 0.0]")
    .replace("[solver]", '[boundary.back]\ntype = "symmetry"\n[boundary.front]\n'
             'type = "symmetry"\n[solver]'))
cas("M20_bom", base(), bom=True)
cas("M20_latin1", base().replace("[physics]", "# viscosité\n[physics]"), encodage="latin-1")
(OUT / "fichier_existant").write_text("x", encoding="utf-8")
f = OUT / "M20_sortie_fichier.toml"
f.write_text(base(), encoding="utf-8")
run("M20_sortie_fichier", ["run2d", str(f), "-o", str(OUT / "fichier_existant"), "--no-plot",
                           "-q"])
msh = OUT / "binaire.msh"
msh.write_bytes(b"$MeshFormat\n4.1 1 8\n\x01\x00\x00\x00\n$EndMeshFormat\n")
cas("M20_gmsh_binaire", base().split("[physics]")[0].split("[mesh]")[0]
    + f'[mesh]\ntype = "file"\npath = "{msh.as_posix()}"\n[physics]'
    + base().split("[physics]")[1])
(OUT / "cavite_cubique_re100_3d").mkdir(exist_ok=True)
run("M22_dossier_masque", ["run2d", "cavite_cubique_re100_3d", "-o", str(OUT / "res_M22"),
                           "--no-plot", "-q", "--set", "solver.max_iter=2"])
run("M19_urans_steps", ["urans", "--steps", "50"])
run("M19_trailing_edge", ["run2d", "cylindre_re20", "-o", str(OUT / "res_M19"),
                           "--no-plot", "-q", "--set", "bodies.0.trailing_edge=open",
                           "solver.max_iter=2"])
run("L6_moyenne_3d", ["run2d", "canal_turbulent_3d", "-o", str(OUT / "res_L6"), "--no-plot",
                      "--set", "solver.max_iter=3"])
print("fin")
