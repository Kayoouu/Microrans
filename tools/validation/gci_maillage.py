"""Étude de convergence en maillage (GCI, Celik et al. 2008) sur trois cas de validation du
README : cylindre Re = 20 (C_d), cavité Re = 100 (u minimal sur la verticale x = 0.5,
Ghia et al. 1982), conduite carrée laminaire 3D (vitesse débitante, solution exacte).
Trois maillages par cas, rapport de raffinement 2, résidus serrés (erreur d'itération
négligeable devant l'erreur de maillage). En 3D le résidu plafonne vers 1e-9 : tol = 1e-11
n'est jamais atteint (conduite 32² : 40 000 itérations) ; à 1e-9 la vitesse débitante est à
2e-7 près de celle des 40 000 itérations, devant ~1e-2 d'écart entre maillages.
Le cylindre est aussi calculé avec un domaine 2 fois plus petit et 2 fois plus grand
(rayon du champ lointain 20, 40, 80) : effet de la taille du domaine sur C_d.

    python tools/validation/gci_maillage.py <dossier de sortie>

Écrit gci.json et affiche le tableau du README (§ 6, « Convergence en maillage »).
"""
import copy
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from microrans.cli import examples_dir  # noqa: E402
from microrans.fv2d.case import run_case  # noqa: E402
from microrans.fv2d.sampling import Sampler, line_points  # noqa: E402
from microrans.gci import gci  # noqa: E402
from microrans.mesh2d.builder import load_config  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)


def duct_exact(nu=0.01, f=1.0, a=0.5, b=0.5, terms=2000):
    """Vitesse débitante exacte de la conduite rectangulaire 2a × 2b (White, éq. 3-48)."""
    s = sum(math.tanh(i * math.pi * b / (2 * a)) / i ** 5 for i in range(1, 2 * terms, 2))
    q = 4 * b * a ** 3 * f / (3 * nu) * (1 - 192 * a / (math.pi ** 5 * b) * s)
    return q / (4 * a * b)


def run(cfg, name):
    t0 = time.perf_counter()
    summary, solver = run_case(cfg, out_dir=OUT / name, verbose=False, plot=False,
                               return_solver=True)
    return summary, solver, time.perf_counter() - t0


def cylinder(level):
    c = copy.deepcopy(load_config(examples_dir() / "cylindre_re20.toml"))
    k = 2 ** level                                  # 0 : 48 × 32, 1 : 96 × 64, 2 : 192 × 128
    c["mesh"].update(n_around=48 * k, n_radial=32 * k, first_height=0.02 / k)
    c["solver"].update(tol=1e-10, max_iter=20000)
    s, solver, t = run(c, f"cylindre_{level}")
    return s["cylinder"]["Cd"], solver.mesh.n_cells, s, t


def cavity(level):
    c = copy.deepcopy(load_config(examples_dir() / "cavite_re100.toml"))
    n = 32 * 2 ** level                             # 32², 64², 128²
    c["mesh"].update(nx=n, ny=n)
    c["solver"].update(tol=1e-10, max_iter=40000)
    s, solver, t = run(c, f"cavite_{level}")
    _, pts = line_points((0.5, 0.0), (0.5, 1.0), 4001)
    u = np.asarray(Sampler(solver, pts).sample(["Ux"])["Ux"])
    # point tabulé par Ghia et al. le plus proche du minimum (y = 0.4531)
    s["u_ghia_point"] = float(Sampler(solver, np.array([[0.5, 0.4531]])).sample(["Ux"])["Ux"][0])
    return float(u.min()), solver.mesh.n_cells, s, t


def duct(level):
    c = copy.deepcopy(load_config(examples_dir() / "conduite_carree_3d.toml"))
    n = 16 * 2 ** level                             # 16², 32², 64² (2 mailles selon x)
    c["mesh"].update(ny=n, nz=n)
    c["solver"].update(tol=1e-9, max_iter=40000)
    c["output"].pop("lines", None)
    s, solver, t = run(c, f"conduite_{level}")
    return s["U_mean"][0], solver.mesh.n_cells, s, t


CASES = [
    ("Cylindre Re = 20, O-grid", "C_d", cylinder, 2, 2.045,
     "2.045 (Dennis & Chang 1970)"),
    ("Cavité Re = 100", "u min sur x = 0.5", cavity, 2, -0.21090, "−0.21090 (Ghia et al. 1982)"),
    # raffinée en y et z seulement (x périodique, écoulement invariant) : r = 2, d'où dim = 2
    ("Conduite carrée laminaire 3D", "vitesse débitante", duct, 2, duct_exact(),
     f"{duct_exact():.5f} (série exacte, White)"),
]

results = {}
for title, what, fn, dim, ref, ref_txt in CASES:
    runs = [fn(level) for level in (2, 1, 0)]       # du plus fin au plus grossier
    phi = [r[0] for r in runs]
    cells = [r[1] for r in runs]
    g = gci(phi, cells, dim)
    band = g["gci_fine"] * abs(phi[0])
    results[title] = {
        "grandeur": what, "phi": phi, "cellules": cells, "reference": ref,
        "iterations": [r[2].get("iterations") for r in runs],
        "converge": [r[2].get("converged") for r in runs],
        "temps_s": [round(r[3], 1) for r in runs], **g,
        "u_ghia_point": [r[2].get("u_ghia_point") for r in runs],
        "ecart_fin_reference": (phi[0] - ref) / abs(ref),
        "ecart_extrapole_reference": (g["phi_ext"] - ref) / abs(ref),
        "reference_dans_bande": abs(ref - phi[0]) <= band}
    print(f"== {title} : {what}", flush=True)
    for k, v in results[title].items():
        print(f"   {k} : {v}")

# taille du domaine (cylindre, maillage moyen, même première maille et même progression
# radiale : n_radial ajusté)
dom = {}
for radius, n_radial in ((20.0, 57), (40.0, 64), (80.0, 72)):
    c = copy.deepcopy(load_config(examples_dir() / "cylindre_re20.toml"))
    c["mesh"].update(n_around=96, n_radial=n_radial, first_height=0.01,
                     farfield_radius=radius)
    c["solver"].update(tol=1e-10, max_iter=20000)
    s, solver, t = run(c, f"domaine_{radius:g}")
    dom[radius] = {"Cd": s["cylinder"]["Cd"], "cellules": solver.mesh.n_cells,
                   "iterations": s.get("iterations")}
    print(f"== domaine R = {radius:g} : {dom[radius]}", flush=True)
# même extrapolation que pour le maillage, en 1/R (rapport 2 : « cellules » = R, dim = 1)
dom["extrapolation"] = gci([dom[r]["Cd"] for r in (80.0, 40.0, 20.0)], (80, 40, 20), 1)
print(f"== domaine, extrapolation R → ∞ : {dom['extrapolation']}", flush=True)
results["cylindre_taille_domaine"] = dom

(OUT / "gci.json").write_text(json.dumps(results, indent=1, ensure_ascii=False),
                              encoding="utf-8")
print("\n| Cas | Grandeur | 3 maillages (cellules) | φ₃ / φ₂ / φ₁ | ordre p | extrapolé | "
      "GCI fin | référence |")
print("|---|---|---|---|---|---|---|---|")
for title, r in results.items():
    if "phi" not in r:
        continue
    n3, n2, n1 = r["cellules"][::-1]
    f3, f2, f1 = r["phi"][::-1]
    osc = " (oscillant)" if r["oscillating"] else ""
    print(f"| {title} | {r['grandeur']} | {n3} / {n2} / {n1} | {f3:.5g} / {f2:.5g} / "
          f"{f1:.5g} | {r['p']:.2f}{osc} | {r['phi_ext']:.5g} | {100 * r['gci_fine']:.2g} % | "
          f"{r['reference']:.5g} |")
