"""Campagne 13 : les fichiers wall_<frontière>.csv permettent-ils de retrouver les efforts du
résumé ? Intégration de C_p sur les faces (aires et normales prises dans le maillage, car le
CSV ne les donne pas : constat L7) et comparaison avec summary.json (60 itérations)."""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

from microrans.fv2d.case import run_case  # noqa: E402
from microrans.tomlio import loads  # noqa: E402

EX = REPO / "microrans" / "examples"
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)

for name, patch in (("cylindre_re20", "cylinder"), ("naca0012_sa", "airfoil"),
                    ("conduite_carree_3d", "bottom")):
    cfg = loads((EX / f"{name}.toml").read_text(encoding="utf-8"))
    cfg["solver"]["max_iter"] = 60
    out = OUT / f"o_{name}"
    s, solver = run_case(cfg, base_dir=EX, out_dir=out, verbose=False, plot=False,
                         return_solver=True)
    m = solver.mesh
    p = next(q for q in m.patches if q.name == patch)
    Sf = np.asarray(m.Sf[p.faces])
    A = np.linalg.norm(Sf, axis=1)
    n = Sf / A[:, None]
    csv = np.genfromtxt(out / f"wall_{patch}.csv", delimiter=",", names=True)
    Lref = s["reference_length"]
    Aref = s.get("reference_area")
    if Aref is None:
        Aref = Lref if m.points.shape[1] == 2 else Lref * np.ptp(m.points[:, 2])
    Fp = (csv["Cp"][:, None] * n * A[:, None]).sum(0) / Aref
    print(f"{name} : Cd pression CSV {Fp[0]:+.6f} / résumé {s[patch]['Cd_pressure']:+.6f} ; "
          f"Σ Cf·A / A_ref {np.sum(csv['Cf'] * A) / Aref:+.6f} / résumé "
          f"{s[patch]['Cd_viscous']:+.6f}")
