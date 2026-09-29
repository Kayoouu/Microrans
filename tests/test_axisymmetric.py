"""Écoulements axisymétriques : solutions exactes et références publiées."""
import numpy as np
import pytest
from scipy.optimize import brentq

from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.case import run_case
from microrans.mesh2d import Circle, o_grid, rectangle_mesh
from microrans.mesh2d.blocks import grading_for_first_cell


def _pipe(ny, nx=2, length=1.0, radius=1.0, grading=1.0, periodic=True):
    return rectangle_mesh(0, length, 0, radius, nx, ny, grading=(1.0, grading),
                          names={"left": "inlet", "right": "outlet", "bottom": "axis",
                                 "top": "wall"}, types={"top": "wall"},
                          periodic=[("inlet", "outlet")] if periodic else None)


def test_gradient_of_constant_and_linear_fields():
    """Gradient de Green-Gauss pondéré par r : exact pour les champs linéaires (la
    contribution des faces latérales du secteur est incluse)."""
    m = _pipe(8, 8, periodic=False)
    s = Solver2D(m, 0.1, {"axis": {"type": "axis"}, "wall": {"type": "wall"},
                          "inlet": {"type": "wall"}, "outlet": {"type": "wall"}},
                 axisymmetric=True)
    x, y = m.cell_centers.T
    fb = m.face_centers[m.n_internal:]
    for f, fbv, exact in [(np.ones_like(x), np.ones(len(fb)), (0.0, 0.0)),
                          (2 * x + 3 * y, 2 * fb[:, 0] + 3 * fb[:, 1], (2.0, 3.0))]:
        g = s.fvm.grad(f, fbv)
        interior = (y > 0.2) & (y < 0.8) & (x > 0.2) & (x < 0.8)    # hors frontières
        assert np.allclose(g[interior], exact, atol=1e-10)


def test_hagen_poiseuille_second_order():
    """Tuyau : u(r) = f (R² − r²) / (4ν), débit π f R⁴ / (8ν)."""
    nu, f = 0.1, 1.0
    err = []
    for ny in (8, 16, 32):
        m = _pipe(ny)
        s = Solver2D(m, nu, {"axis": {"type": "axis"}, "wall": {"type": "wall"}},
                     body_force=(f, 0.0), settings=Settings(relax_U=1.0), axisymmetric=True)
        assert s.run_steady(max_iter=300, tol=1e-11)
        r = m.cell_centers[:, 1]
        ex = f / (4 * nu) * (1 - r ** 2)
        err.append(np.max(np.abs(s.U[:, 0] - ex)) / ex.max())
        assert np.abs(s.U[:, 1]).max() < 1e-12
    order = np.log2(np.array(err[:-1]) / err[1:])
    assert np.all(order > 1.95), order
    q = 2 * np.pi * np.sum(s.U[:, 0] * s.fvm.V)             # longueur 1
    assert q == pytest.approx(np.pi * f / (8 * nu), rel=2e-3)


def test_radial_source_flow_exercises_hoop_terms():
    """u_r = C/r, u_x = 0, p = −C²/(2r²) : solution exacte de Navier-Stokes où le terme
    circonférentiel −2ν u_r/r² compense exactement le reste de ∇·τ (ν grand : un terme
    circonférentiel faux ou absent donne une vitesse fausse)."""
    C, r1, r2, nu = 0.5, 0.5, 2.0, 1.0
    eu, ep = [], []
    for ny in (16, 32):
        m = rectangle_mesh(0, 0.5, r1, r2, 2, ny, names={"left": "a", "right": "b",
                                                         "bottom": "inner", "top": "outer"},
                           periodic=[("a", "b")])
        s = Solver2D(m, nu, {"inner": {"type": "inlet", "U": [0.0, C / r1]},
                             "outer": {"type": "inlet", "U": [0.0, C / r2]}},
                     settings=Settings(convection_U="linearUpwind"), axisymmetric=True,
                     reference_velocity=C / r1)
        assert s.run_steady(max_iter=3000, tol=1e-10)
        r = m.cell_centers[:, 1]
        eu.append(np.max(np.abs(s.U[:, 1] - C / r)) / (C / r1))
        k = (r > 0.8) & (r < 1.7)            # loin des bords (gradient de p nul imposé)
        pex = -C ** 2 / (2 * r ** 2)
        d = (s.p - s.p[k].mean()) - (pex - pex[k].mean())
        ep.append(np.max(np.abs(d[k])) / np.ptp(pex))
    assert eu[1] < 5e-3 and eu[0] / eu[1] > 3.3
    assert ep[1] < 2e-3 and ep[0] / ep[1] > 3.0


def _sphere_case():
    return {"mesh": {"type": "ogrid", "n_around": 128, "n_radial": 64, "farfield_radius": 30.0,
                     "first_height": 0.005, "cut_axis": True},
            "bodies": [{"type": "circle", "center": [0, 0], "radius": 0.5, "name": "sphere"}],
            "physics": {"reynolds": 100, "reference_velocity": 1.0, "reference_length": 1.0,
                        "axisymmetric": True},
            "initial": {"U": [1.0, 0.0]},
            "boundary": {"sphere": {"type": "wall"}, "axis": {"type": "axis"},
                         "farfield": {"type": "farfield", "U": [1.0, 0.0]}},
            "solver": {"max_iter": 2000, "tol": 1e-6},
            "output": {"plots": False, "vtk": False, "checkpoint": False}}


def test_sphere_drag_re100(tmp_path):
    """Sphère, Re = 100 : Cd = 1.085 (Fornberg 1988) ; corrélation de Clift et al. 1.087.
    Maillage en O coupé à l'axe ([mesh] cut_axis)."""
    s = run_case(_sphere_case(), out_dir=tmp_path, verbose=False, plot=False)
    assert s["converged"] and s["axisymmetric"]
    assert s["sphere"]["Cd"] == pytest.approx(1.085, rel=0.015)
    assert s["sphere"]["Cl"] == 0.0


def test_laminar_pipe_nusselt_uniform_heat_flux():
    """Écoulement établi, flux pariétal uniforme : Nu = 48/11 = 4.364 (Shah & London)."""
    R, L, nu = 0.5, 10.0, 1.0 / 50
    m = _pipe(24, 120, L, R, periodic=False)
    bcs = {"inlet": {"type": "inlet", "U": ["2*(1-(y/0.5)**2)", 0.0], "T": 0.0},
           "outlet": {"type": "outlet"}, "axis": {"type": "axis"},
           "wall": {"type": "wall", "q": 1.0}}
    s = Solver2D(m, nu, bcs, energy={"Pr": 1.0}, axisymmetric=True, reference_velocity=1.0,
                 settings=Settings(relax_T=1.0))
    assert s.run_steady(max_iter=2000, tol=1e-8)
    alpha = nu
    C, V = m.cell_centers, s.fvm.V
    col = np.isclose(C[:, 0], C[np.argmin(np.abs(C[:, 0] - 8.0)), 0])
    Tb = np.sum(s.U[col, 0] * s.T[col] * V[col]) / np.sum(s.U[col, 0] * V[col])
    Tw, q = s.wall_heat_flux("wall")
    xf = m.face_centers[m.n_internal:][s.patch_slices["wall"]][:, 0]
    nu_x = q.mean() * 2 * R / (alpha * (np.interp(C[col, 0][0], xf, Tw) - Tb))
    assert nu_x == pytest.approx(48 / 11, rel=3e-3)
    # bilan d'énergie : puissance pariétale = débit × élévation de température moyenne
    heat = 2 * np.pi * np.sum(q * s.fvm.magSb[s.patch_slices["wall"]])
    assert heat == pytest.approx(2 * np.pi * R * L * 1.0, rel=1e-9)


def test_turbulent_pipe_friction_sa():
    """Tuyau turbulent lisse, Re_τ = 550 (Re_D ≈ 19 000), SA : λ à ±4 % de la loi de
    Prandtl 1/√λ = 2 log10(Re √λ) − 0.8 (+2.5 % mesuré)."""
    re_tau, ny = 550.0, 48
    r = grading_for_first_cell(ny, 0.5 / re_tau)
    m = _pipe(ny, grading=1.0 / r)
    s = Solver2D(m, 1 / re_tau, {"axis": {"type": "axis"}, "wall": {"type": "wall"}},
                 model="sa", body_force=(2.0, 0.0), initial_U=(15.0, 0.0),
                 reference_velocity=15.0, axisymmetric=True,
                 turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0})
    assert s.run_steady(max_iter=6000, tol=1e-7)
    ub = np.sum(s.U[:, 0] * s.fvm.V) / np.sum(s.fvm.V)        # u_τ = 1
    re = ub * 2 * re_tau
    lam = 8 / ub ** 2
    ref = brentq(lambda x: 1 / x - (2.0 * np.log10(re * x) - 0.8), 1e-3, 1.0) ** 2
    assert lam == pytest.approx(ref, rel=0.04)


def test_axisymmetric_input_checks():
    m = rectangle_mesh(0, 1, -1, 1, 4, 4, names={"left": "w", "right": "w", "bottom": "w",
                                                 "top": "w"})
    with pytest.raises(ValueError, match="demi-plan"):
        Solver2D(m, 0.1, {"w": {"type": "wall"}}, axisymmetric=True)
    m = _pipe(4)
    with pytest.raises(ValueError, match="axe x"):
        Solver2D(m, 0.1, {"axis": {"type": "axis"}, "wall": {"type": "wall"}},
                 body_force=(0.0, 1.0), axisymmetric=True)
    # coupe impossible : pas d'arêtes sur y = 0 (nombre impair de points autour du corps)
    with pytest.raises(ValueError, match="cut_at_axis"):
        o_grid(Circle((0, 0), 0.5, "c"), 33, 8, 5.0, 0.05).cut_at_axis()
    half = o_grid(Circle((0, 0), 0.5, "c"), 32, 8, 5.0, 0.05).cut_at_axis()
    assert half.n_cells == 16 * 8 and half.patch("axis").type == "symmetry"
    cfg = _sphere_case()
    cfg["physics"]["angle_of_attack"] = 5.0
    with pytest.raises(ValueError, match="incidence"):
        run_case(cfg, out_dir=None, verbose=False, plot=False)
