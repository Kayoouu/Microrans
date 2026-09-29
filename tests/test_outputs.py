"""Démarrage multigrille, sondes, profils, moyennes temporelles."""
import copy
import csv

import numpy as np
import pytest

from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.case import run_case
from microrans.fv2d.fmg import coarsen_config
from microrans.fv2d.sampling import Sampler, line_points, locate, parse_points
from microrans.mesh2d import Circle, Rectangle, channel_mesh, triangulate


def _cavity(n=32, **output):
    return {"mesh": {"type": "rectangle", "x0": 0.0, "x1": 1.0, "y0": 0.0, "y1": 1.0,
                     "nx": n, "ny": n,
                     "names": {"left": "walls", "right": "walls", "bottom": "walls",
                               "top": "lid"}},
            "physics": {"nu": 0.01, "reference_velocity": 1.0},
            "boundary": {"lid": {"type": "wall", "U": [1.0, 0.0]}, "walls": {"type": "wall"}},
            "solver": {"max_iter": 3000, "tol": 1e-6},
            "output": {"vtk": False, "plots": False, **output}}


def test_coarsen_config():
    c = coarsen_config(_cavity(64), 2)
    assert (c["mesh"]["nx"], c["mesh"]["ny"]) == (16, 16)
    og = {"mesh": {"type": "ogrid", "n_around": 130, "n_radial": 64, "first_height": 1e-3}}
    c = coarsen_config(og, 1)["mesh"]
    assert c["n_around"] % 2 == 0 and c["n_radial"] == 32 and c["first_height"] == 2e-3
    hy = {"mesh": {"type": "hybrid", "h_max": 1.0, "h_surface": 0.05,
                   "layers": {"n": 12, "first_height": 1e-4}}}
    c = coarsen_config(hy, 1)["mesh"]
    assert c["h_max"] == 2.0 and c["layers"]["n"] == 6 and c["layers"]["first_height"] == 2e-4
    bl = {"mesh": {"type": "blocks", "blocks": [{"cells": [40, 21]}]}}
    assert coarsen_config(bl, 1)["mesh"]["blocks"][0]["cells"] == [20, 10]
    assert coarsen_config({"mesh": {"type": "file", "path": "x.msh"}}, 1) is None


def test_fmg_same_solution_fewer_fine_iterations(tmp_path):
    ref = run_case(_cavity(64), out_dir=tmp_path / "a", verbose=False, plot=False)
    cfg = _cavity(64)
    cfg["solver"]["fmg_levels"] = 2
    s = run_case(cfg, out_dir=tmp_path / "b", verbose=False, plot=False)
    assert [lv["n_cells"] for lv in s["fmg"]["levels"]] == [256, 1024]
    assert s["restart"]["mode"] == "interpolé"
    assert s["iterations"] < 0.7 * ref["iterations"]
    assert s["lid"]["Cd"] == pytest.approx(ref["lid"]["Cd"], rel=1e-4)


def test_locate_matches_matplotlib_on_triangles():
    from matplotlib.path import Path
    m = triangulate(Rectangle(0, 0, 2, 1) - Circle((1.0, 0.5), 0.25), 0.15, max_iter=80)
    rng = np.random.default_rng(0)
    pts = rng.uniform([0, 0], [2, 1], size=(300, 2))
    cell = locate(m, pts)
    inside_body = np.hypot(pts[:, 0] - 1.0, pts[:, 1] - 0.5) < 0.24
    assert np.all(cell[inside_body] == -1)
    for p, c in zip(pts[cell >= 0], cell[cell >= 0]):
        poly = m.points[m.cell_nodes[c, :m.cell_nv[c]]]
        assert Path(poly).contains_point(p, radius=1e-9) or Path(poly).contains_point(
            p, radius=-1e-9)


def test_line_profile_poiseuille_and_wall_values():
    """Reconstruction φ_P + ∇φ·d : profil parabolique à 0.8 % près avec 16 mailles (erreur
    d'une reconstruction linéaire dans la maille), valeurs de paroi ≈ 0 ; points hors
    domaine = NaN."""
    m = channel_mesh(1.0, 2.0, 2, 16)
    s = Solver2D(m, 0.1, {"bottom": {"type": "wall"}, "top": {"type": "wall"}},
                 body_force=(1.0, 0.0), settings=Settings(relax_U=1.0))
    s.run_steady(max_iter=200, tol=1e-10)
    _, pts = line_points((0.5, 0.0), (0.5, 2.0), 41)
    v = Sampler(s, pts).sample(["Ux", "p"])
    y = pts[:, 1]
    ex = y * (2 - y) / 0.2
    assert np.max(np.abs(v["Ux"] - ex)) < 1e-2 * ex.max()
    assert abs(v["Ux"][0]) < 1e-2 * ex.max() and abs(v["Ux"][-1]) < 1e-2 * ex.max()
    assert np.isnan(Sampler(s, [[0.5, 3.0]]).sample(["Ux"])["Ux"][0])


def test_probes_and_lines_in_case_outputs(tmp_path):
    cfg = _cavity(24, probes="0.5 0.5; 5 5",
                  lines=[{"name": "v", "start": [0.5, 0.0], "end": [0.5, 1.0], "n": 50}])
    cfg["output"]["plots"] = True
    s = run_case(cfg, out_dir=tmp_path, verbose=False, plot=True)
    assert np.isfinite(s["probes"][0]["Ux"]) and np.isnan(s["probes"][1]["Ux"])
    with open(tmp_path / "history.csv", encoding="utf-8") as fh:
        row = list(csv.DictReader(fh))[-1]
    assert float(row["probe1_Ux"]) == pytest.approx(s["probes"][0]["Ux"], rel=1e-6)
    d = np.genfromtxt(tmp_path / "line_v.csv", delimiter=",", names=True)
    assert len(d) == 50 and d["Ux"][-1] == pytest.approx(1.0, abs=0.02)
    assert (tmp_path / "line_v.png").exists()
    assert parse_points("1, 2; 3 4").tolist() == [[1, 2], [3, 4]]


def test_time_averages_and_exact_restart(tmp_path):
    cfg = _cavity(16, average_from=0.5, checkpoint_minutes=0)
    cfg["solver"].update(mode="transient", dt=0.05, t_end=2.0)
    _, a = run_case(copy.deepcopy(cfg), out_dir=tmp_path / "a", verbose=False, plot=False,
                    return_solver=True)
    c1 = copy.deepcopy(cfg)
    c1["solver"]["t_end"] = 1.2
    run_case(c1, out_dir=tmp_path / "b", verbose=False, plot=False)
    c2 = copy.deepcopy(cfg)
    c2["initial"] = {"restart": str(tmp_path / "b" / "checkpoint.npz")}
    _, b = run_case(c2, out_dir=tmp_path / "c", verbose=False, plot=False, return_solver=True)
    ma, mb = a.mean_fields(), b.mean_fields()
    assert a.averager.weight == pytest.approx(1.5)
    assert set(ma) == {"Ux_mean", "Ux_rms", "Uy_mean", "Uy_rms", "p_mean", "p_rms"}
    for k in ma:
        assert np.array_equal(ma[k], mb[k])
    # démarrage de la cavité : le champ évolue, donc les fluctuations ne sont pas nulles
    assert ma["Ux_rms"].max() > 1e-3
