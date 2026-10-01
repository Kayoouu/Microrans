"""Solveur volumes finis 2D : vérification (solutions exactes) et validation (références publiées)."""
import json

import numpy as np
import pytest
from scipy.interpolate import LinearNDInterpolator

from microrans.cases import run_rans_channel
from microrans.cli import main
from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.benchmarks import GHIA_RE100_U, GHIA_RE100_V
from microrans.mesh2d import Circle, cavity_mesh, channel_mesh, o_grid, rectangle_mesh
from microrans.reference import womersley_channel

WALLS = {"bottom": {"type": "wall"}, "top": {"type": "wall"}}


def poiseuille_error(ny):
    m = channel_mesh(1.0, 2.0, 2, ny)
    s = Solver2D(m, 0.1, WALLS, body_force=(1.0, 0.0), settings=Settings(relax_U=1.0))
    assert s.run_steady(max_iter=200, tol=1e-10)
    y = m.cell_centers[:, 1]
    ex = y * (2 - y) / (2 * 0.1)
    return np.max(np.abs(s.U[:, 0] - ex)) / ex.max()


def test_periodic_poiseuille_second_order():
    e1, e2 = poiseuille_error(16), poiseuille_error(32)
    assert e2 < 2e-3
    assert np.log2(e1 / e2) == pytest.approx(2.0, abs=0.2)


def test_developing_channel_reaches_parabolic_profile():
    m = rectangle_mesh(0, 10, 0, 1, 80, 20, names={"left": "inlet", "right": "outlet",
                                                    "bottom": "wall", "top": "wall"})
    s = Solver2D(m, 0.01, {"inlet": {"type": "inlet", "U": [1, 0]},
                           "outlet": {"type": "outlet"}, "wall": {"type": "wall"}})
    assert s.run_steady(max_iter=1000, tol=1e-7)
    C = m.cell_centers
    sel = C[:, 0] > 9.8
    ex = 6 * C[sel, 1] * (1 - C[sel, 1])
    assert np.max(np.abs(s.U[sel, 0] - ex)) < 0.01
    # conservation de la masse : débit de sortie = débit d'entrée
    fo = s.F_b[s.patch_slices["outlet"]].sum()
    fi = s.F_b[s.patch_slices["inlet"]].sum()
    assert fo == pytest.approx(-fi, rel=1e-8)


def test_lid_driven_cavity_ghia_re100():
    m = cavity_mesh(48)
    s = Solver2D(m, 0.01, {"lid": {"type": "wall", "U": [1, 0]}, "walls": {"type": "wall"}},
                 reference_velocity=1.0)
    assert s.run_steady(max_iter=1500, tol=1e-6)
    C = m.cell_centers
    iu = LinearNDInterpolator(C, s.U[:, 0])
    iv = LinearNDInterpolator(C, s.U[:, 1])
    yu, xv = GHIA_RE100_U[1:-1, 0], GHIA_RE100_V[1:-1, 0]
    eu = iu(np.column_stack([np.full_like(yu, 0.5), yu])) - GHIA_RE100_U[1:-1, 1]
    ev = iv(np.column_stack([xv, np.full_like(xv, 0.5)])) - GHIA_RE100_V[1:-1, 1]
    assert np.nanmax(np.abs(eu)) < 0.015
    assert np.nanmax(np.abs(ev)) < 0.015


def _womersley(steps, scheme):
    nu, A, om, T = 0.02, 5.0, 2 * np.pi, 1.0
    m = channel_mesh(1.0, 2.0, 2, 128)
    y = m.cell_centers[:, 1]
    s = Solver2D(m, nu, WALLS, settings=Settings(time_scheme=scheme), reference_velocity=1.0)
    s.U[:, 0] = womersley_channel(y, 0.0, 1.0, A, om, nu)
    s.F_i = np.sum(s.fvm.interp(s.U) * s.fvm.Si, axis=1)
    dt = T / steps
    s.body_force = np.array([1.0 + A * np.sin(om * dt), 0.0])

    def cb(sv, n):
        sv.body_force = np.array([1.0 + A * np.sin(om * (sv.time + dt)), 0.0])
    s.run_transient(dt, T, callback=cb)
    ex = womersley_channel(y, T, 1.0, A, om, nu)
    return np.max(np.abs(s.U[:, 0] - ex)) / np.max(np.abs(ex))


def test_pimple_backward_is_second_order_in_time():
    e1, e2 = _womersley(20, "backward"), _womersley(40, "backward")
    assert np.log2(e1 / e2) > 1.8


def test_turbulent_channel_2d_matches_1d_sa():
    re_tau = 395.0
    m = channel_mesh(1.0, 2.0, 2, 96, first_height=0.4 / re_tau)
    s = Solver2D(m, 1 / re_tau, WALLS, model="sa", body_force=(1.0, 0.0), initial_U=(15.0, 0.0),
                 turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0})
    assert s.run_steady(max_iter=4000, tol=1e-7)
    ub = np.sum(s.U[:, 0] * m.cell_volumes) / np.sum(m.cell_volumes)
    r1 = run_rans_channel("sa", re_tau, n_cells=256, y1_plus=0.2).summary
    assert ub == pytest.approx(r1["Ub_plus"] * r1["u_tau"], rel=1e-3)
    _, tau, _ = s.wall_shear("bottom")
    assert tau.mean() == pytest.approx(1.0, abs=1e-3)


@pytest.mark.parametrize("model", ["sst", "kw", "ke"])
def test_turbulent_channel_2d_other_models_converge(model):
    re_tau = 395.0
    m = channel_mesh(1.0, 2.0, 2, 96, first_height=0.4 / re_tau)
    s = Solver2D(m, 1 / re_tau, WALLS, model=model, body_force=(1.0, 0.0), initial_U=(15.0, 0.0),
                 turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0})
    assert s.run_steady(max_iter=4000, tol=1e-6)
    ub = np.sum(s.U[:, 0] * m.cell_volumes) / np.sum(m.cell_volumes)
    # même équations que le 1D : écart < 2 % à ce maillage (sensibilité à y1+, voir README)
    r1 = run_rans_channel(model, re_tau, n_cells=1024, y1_plus=0.05, max_iter=50000).summary
    assert ub == pytest.approx(r1["Ub_plus"] * r1["u_tau"], rel=0.02)


def test_cylinder_re20_drag():
    m = o_grid(Circle((0, 0), 0.5, "cylinder"), 80, 56, 40.0, 0.01)
    s = Solver2D(m, 1 / 20, {"cylinder": {"type": "wall"},
                             "farfield": {"type": "farfield", "U": [1, 0]}}, initial_U=(1, 0))
    assert s.run_steady(max_iter=1500, tol=1e-6)
    f = s.forces()["cylinder"]["total"]
    cd = f[0] / 0.5
    assert 2.0 < cd < 2.1            # Dennis & Chang (1970) : 2.045
    assert abs(f[1] / 0.5) < 1e-6    # portance nulle (symétrie)


def test_hybrid_mesh_does_not_diverge():
    """Non-régression : SIMPLEC divergeait sur maillage hybride (non-orthogonalité ~60°)
    tant que le terme (rAtU − rAU)∂p/∂n ignorait la correction non orthogonale."""
    import warnings

    from microrans.mesh2d import Rectangle, hybrid_mesh
    dom = Rectangle(-6, -6, 12, 6, names={"left": "inlet", "right": "outlet",
                                          "bottom": "side", "top": "side"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = hybrid_mesh(dom, [Circle((0, 0), 0.5, "cylinder").as_wall()], 1.5, 0.1,
                        n_layers=5, first_height=0.02, ratio=1.2, max_iter=150)
    bcs = {"cylinder": {"type": "wall"}, "inlet": {"type": "inlet", "U": [1, 0]},
           "outlet": {"type": "outlet"}, "side": {"type": "symmetry"}}
    s = Solver2D(m, 1 / 20, bcs, initial_U=(1, 0), settings=Settings(relax_U=0.9))
    s.run_steady(max_iter=400, tol=1e-5)          # ne doit pas lever FloatingPointError
    cd = s.forces()["cylinder"]["total"][0] / 0.5
    assert 2.0 < cd < 2.8                          # domaine confiné : Cd > valeur non confinée
    assert s.history[-1]["Ux"] < 1e-3


def test_cli_run2d(tmp_path):
    case = tmp_path / "cav.toml"
    case.write_text("""
[mesh]
type = "rectangle"
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
nx = 16
ny = 16
names = { left = "walls", right = "walls", bottom = "walls", top = "lid" }
[physics]
nu = 0.01
reference_velocity = 1.0
[boundary.lid]
type = "wall"
U = [1.0, 0.0]
[boundary.walls]
type = "wall"
[solver]
max_iter = 50
""")
    out = tmp_path / "out"
    main(["run2d", str(case), "-o", str(out), "-q", "--set", "solver.max_iter=400"])
    s = json.loads((out / "summary.json").read_text())
    assert s["converged"] and s["iterations"] > 50
    for f in ("fields.vtk", "history.csv", "U.png", "convergence.png"):
        assert (out / f).exists()


@pytest.mark.parametrize("model,tol", [("sa", 0.04), ("sst", 0.06), ("kw", 0.02)])
def test_wall_functions_coarse_channel(model, tol):
    """Lois de paroi (Spalding) : 1re cellule à y⁺ ≈ 50, Re_τ = 2000 ; débit comparé au même
    modèle résolu jusqu'à la paroi (1D, y⁺ = 0.2). L'écart vient surtout de la loi log
    universelle (κ = 0.41, B = 5.2) imposée par la loi de paroi."""
    re_tau = 2000.0
    m = channel_mesh(1.0, 2.0, 2, 24, first_height=100.0 / re_tau)
    s = Solver2D(m, 1 / re_tau, WALLS, model=model, body_force=(1.0, 0.0), initial_U=(20.0, 0.0),
                 turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0},
                 settings=Settings(wall_treatment="wall_function"))
    assert s.run_steady(max_iter=3000, tol=1e-7)
    _, tau, yp = s.wall_shear("bottom")
    assert tau.mean() == pytest.approx(1.0, abs=1e-3)
    assert 40 < yp.mean() < 60
    ub = np.sum(s.U[:, 0] * m.cell_volumes) / np.sum(m.cell_volumes)
    r1 = run_rans_channel(model, re_tau, n_cells=256, y1_plus=0.2).summary
    assert ub == pytest.approx(r1["Ub_plus"] * r1["u_tau"], rel=tol)


def test_wall_functions_rejected_for_low_re_k_epsilon():
    with pytest.raises(ValueError, match="bas-Reynolds"):
        Solver2D(channel_mesh(1.0, 2.0, 2, 16), 1e-3, WALLS, model="ke",
                 settings=Settings(wall_treatment="wall_function"))


def test_pseudo_transient_fine_stretched_channel():
    """Maillage fin et étiré (384 cellules, y⁺ = 0.1) : la sous-relaxation implicite de
    SIMPLEC y demande O(N²) itérations (> 8 000) ; le pseudo-pas local convectif converge
    en ~100 itérations vers la même solution."""
    re_tau = 395.0
    m = channel_mesh(1.0, 2.0, 2, 384, first_height=0.1 / re_tau)
    s = Solver2D(m, 1 / re_tau, WALLS, model="sa", body_force=(1.0, 0.0), initial_U=(15.0, 0.0),
                 turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0},
                 settings=Settings(pseudo_cfl=200.0, relax_turb=1.0))
    assert s.run_steady(max_iter=400, tol=1e-7)
    assert s.iterations < 250
    ub = np.sum(s.U[:, 0] * m.cell_volumes) / np.sum(m.cell_volumes)
    r1 = run_rans_channel("sa", re_tau, n_cells=256, y1_plus=0.2).summary
    assert ub == pytest.approx(r1["Ub_plus"] * r1["u_tau"], rel=2e-3)


def test_divergence_reported_with_tips(tmp_path):
    """Audit : une divergence finissait en « Factor is exactly singular » (matrice de
    pression) ; les vitesses démesurées mais finies (~1e50) n'étaient pas détectées."""
    import tomllib
    from pathlib import Path

    from microrans.fv2d.case import run_case
    ex = Path(__file__).resolve().parent.parent / "microrans" / "examples" / "cavite_re100.toml"
    cfg = tomllib.loads(ex.read_text(encoding="utf-8"))
    cfg["mesh"].update(nx=12, ny=12)
    cfg["physics"]["nu"] = 1e-6
    cfg["solver"].update(max_iter=200, algorithm="SIMPLE", relax_U=1.0, relax_p=1.0)
    cfg["output"] = {"plots": False, "vtk": False}
    with pytest.raises(FloatingPointError, match=r"Le calcul a divergé à l'itération \d+ .*"
                                                 r"Pistes : démarrer en convection_U = \"upwind\""):
        run_case(cfg, out_dir=tmp_path, verbose=False, plot=False)
