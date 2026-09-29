"""Incidence de l'écoulement, polaire et balayage de paramètres."""
import copy
import csv

import numpy as np
import pytest

from microrans.cli import main
from microrans.fv2d.case import run_case
from microrans.fv2d.sweep import parse_values, run_sweep


def test_parse_values():
    assert parse_values("-4:12:4") == [-4.0, 0.0, 4.0, 8.0, 12.0]
    assert parse_values("0:1:0.25") == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert parse_values("10, 20; 40") == [10.0, 20.0, 40.0]
    assert parse_values([1, 2]) == [1.0, 2.0]
    with pytest.raises(ValueError):
        parse_values("0:10:-1")


def _ogrid_case(body, re, n_around, n_radial, first, radius=20.0):
    return {"mesh": {"type": "ogrid", "n_around": n_around, "n_radial": n_radial,
                     "farfield_radius": radius, "first_height": first},
            "bodies": [body],
            "physics": {"reynolds": re, "reference_velocity": 1.0, "reference_length": 1.0},
            "initial": {"U": [1.0, 0.0]},
            "boundary": {body["name"]: {"type": "wall"},
                         "farfield": {"type": "farfield", "U": [1.0, 0.0]}},
            "solver": {"max_iter": 3000, "tol": 1e-6},
            "output": {"plots": False, "vtk": False}}


def test_incidence_rotates_flow_and_force_axes(tmp_path):
    """Cylindre (symétrie de révolution) : Cd dans les axes de l'écoulement indépendant de
    l'incidence, portance nulle."""
    cfg = _ogrid_case({"type": "circle", "center": [0, 0], "radius": 0.5, "name": "cyl"},
                      20, 48, 32, 0.02)
    res = {}
    for a in (0.0, 30.0):
        c = copy.deepcopy(cfg)
        c["physics"]["angle_of_attack"] = a
        res[a] = run_case(c, out_dir=tmp_path / str(a), verbose=False, plot=False)["cyl"]
    assert res[30.0]["Cd"] == pytest.approx(res[0.0]["Cd"], rel=1e-3)
    assert abs(res[30.0]["Cl"]) < 2e-3


def test_symmetric_airfoil_polar(tmp_path):
    """Profil symétrique, maillage symétrique : Cl(−α) = −Cl(α), Cd(−α) = Cd(α),
    Cm(−α) = −Cm(α) ; la continuation ne change pas la solution (hors décrochage)."""
    cfg = _ogrid_case({"type": "naca", "code": "0012", "chord": 1.0, "name": "airfoil"},
                      500, 48, 24, 5e-3, radius=10.0)
    cfg["output"]["moment_center"] = [0.25, 0.0]
    cfg["solver"]["tol"] = 1e-5
    rows = run_sweep(cfg, "physics.angle_of_attack", "-3:3:3", out_dir=tmp_path / "c",
                     verbose=False)
    cl, cd, cm = ([r[f"{k}_airfoil"] for r in rows] for k in ("Cl", "Cd", "Cm"))
    assert all(r["converged"] for r in rows)
    assert cl[0] == pytest.approx(-cl[2], rel=1e-3) and cl[2] > 0.05
    assert abs(cl[1]) < 1e-5 * cl[2] + 1e-8
    assert cd[0] == pytest.approx(cd[2], rel=1e-3) and cd[1] < cd[2]
    assert cm[0] == pytest.approx(-cm[2], rel=1e-2, abs=1e-6)
    assert (tmp_path / "c" / "polaire.png").exists()
    cold = run_sweep(cfg, "physics.angle_of_attack", "0:3:3", out_dir=tmp_path / "f",
                     continuation=False, verbose=False, plot=False)
    assert rows[2]["Cl_airfoil"] == pytest.approx(cold[1]["Cl_airfoil"], rel=1e-3)


def test_cli_sweep_reynolds(tmp_path):
    case = tmp_path / "cav.toml"
    case.write_text("""
[mesh]
type = "rectangle"
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
nx = 12
ny = 12
names = { left = "walls", right = "walls", bottom = "walls", top = "lid" }
[physics]
reynolds = 100
reference_velocity = 1.0
[boundary.lid]
type = "wall"
U = [1.0, 0.0]
[boundary.walls]
type = "wall"
[solver]
max_iter = 2000
tol = 1e-6
""")
    out = tmp_path / "out"
    main(["sweep", str(case), "-o", str(out), "-q", "--param", "physics.reynolds",
          "--values", "10", "100"])
    with open(out / "balayage.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert [float(r["physics.reynolds"]) for r in rows] == [10.0, 100.0]
    assert all(r["converged"] == "True" for r in rows)
    # frottement sur le couvercle : décroît avec Re (Cd = force / (½ U² L))
    cd = [abs(float(r["Cd_lid"])) for r in rows]
    assert np.all(np.isfinite(cd)) and cd[0] > 5 * cd[1]
    assert (out / "balayage.png").exists()
