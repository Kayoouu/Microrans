"""Solveur couplé pression-vitesse ([solver] algorithm = "coupled") : même physique que
SIMPLE(C), solution indépendante de la sous-relaxation, convergence en peu d'itérations."""
import numpy as np
import pytest
from scipy.interpolate import LinearNDInterpolator

from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.benchmarks import GHIA_RE100_U, GHIA_RE100_V
from microrans.mesh2d import Circle, cavity_mesh, channel_mesh, o_grid, rectangle_mesh

CAVITY = {"lid": {"type": "wall", "U": [1, 0]}, "walls": {"type": "wall"}}


def _cavity(n, **kw):
    return Solver2D(cavity_mesh(n), 0.01, CAVITY, reference_velocity=1.0,
                    settings=Settings(**kw))


def test_cavity_ghia_re100_in_few_iterations():
    s = _cavity(48, algorithm="coupled", relax_U=1.0)
    assert s.run_steady(max_iter=100, tol=1e-6)
    assert s.iterations < 40                      # SIMPLEC : plusieurs centaines
    C = s.mesh.cell_centers
    iu = LinearNDInterpolator(C, s.U[:, 0])
    iv = LinearNDInterpolator(C, s.U[:, 1])
    yu, xv = GHIA_RE100_U[1:-1, 0], GHIA_RE100_V[1:-1, 0]
    eu = iu(np.column_stack([np.full_like(yu, 0.5), yu])) - GHIA_RE100_U[1:-1, 1]
    ev = iv(np.column_stack([xv, np.full_like(xv, 0.5)])) - GHIA_RE100_V[1:-1, 1]
    assert np.nanmax(np.abs(eu)) < 0.015
    assert np.nanmax(np.abs(ev)) < 0.015


def test_solution_independent_of_relaxation_and_close_to_simplec():
    a = _cavity(24, algorithm="coupled", relax_U=1.0)
    b = _cavity(24, algorithm="coupled", relax_U=0.7)
    d = _cavity(24, algorithm="coupled", pseudo_cfl=50.0)
    c = _cavity(24, algorithm="SIMPLEC", relax_U=0.9)
    for s in (a, b, c, d):
        assert s.run_steady(max_iter=2000, tol=1e-9)
    # flux de Rhie-Chow avec la diagonale physique (non relaxée, sans pseudo-temps) : même
    # solution à la tolérance près
    assert np.max(np.abs(a.U - b.U)) < 1e-6
    assert np.max(np.abs(a.U - d.U)) < 1e-6
    # SIMPLEC (diagonale relaxée dans Rhie-Chow) : écart de l'ordre de l'erreur de
    # discrétisation seulement (0,5 % de la vitesse du couvercle sur 24 × 24, dans les coins)
    assert np.max(np.abs(a.U - c.U)) < 1e-2


def test_channel_mass_conservation_and_parabolic_profile():
    m = rectangle_mesh(0, 10, 0, 1, 80, 20, names={"left": "inlet", "right": "outlet",
                                                    "bottom": "wall", "top": "wall"})
    s = Solver2D(m, 0.01, {"inlet": {"type": "inlet", "U": [1, 0]},
                           "outlet": {"type": "outlet"}, "wall": {"type": "wall"}},
                 settings=Settings(algorithm="coupled", relax_U=1.0))
    assert s.run_steady(max_iter=200, tol=1e-7)
    C = m.cell_centers
    sel = C[:, 0] > 9.8
    assert np.max(np.abs(s.U[sel, 0] - 6 * C[sel, 1] * (1 - C[sel, 1]))) < 0.01
    fo = s.F_b[s.patch_slices["outlet"]].sum()
    fi = s.F_b[s.patch_slices["inlet"]].sum()
    assert fo == pytest.approx(-fi, rel=1e-8)
    # flux conservatifs : divergence nulle dans chaque cellule
    assert np.max(np.abs(s.fvm.div(s.F_i, s.F_b))) < 1e-8


def test_periodic_poiseuille_body_force():
    m = channel_mesh(1.0, 2.0, 2, 32)
    s = Solver2D(m, 0.1, {"bottom": {"type": "wall"}, "top": {"type": "wall"}},
                 body_force=(1.0, 0.0), settings=Settings(algorithm="coupled", relax_U=1.0))
    assert s.run_steady(max_iter=100, tol=1e-10)
    y = m.cell_centers[:, 1]
    ex = y * (2 - y) / (2 * 0.1)
    assert np.max(np.abs(s.U[:, 0] - ex)) / ex.max() < 2e-3


def test_cylinder_re20_drag():
    m = o_grid(Circle((0, 0), 0.5, "cylinder"), 80, 56, 40.0, 0.01)
    s = Solver2D(m, 1 / 20, {"cylinder": {"type": "wall"},
                             "farfield": {"type": "farfield", "U": [1, 0]}}, initial_U=(1, 0),
                 settings=Settings(algorithm="coupled", relax_U=1.0))
    assert s.run_steady(max_iter=100, tol=1e-6)
    assert s.iterations < 40
    cd = s.forces()["cylinder"]["total"][0] / 0.5
    assert 2.0 < cd < 2.1            # Dennis & Chang (1970) : 2.045
    # factorisation LU réutilisée d'une itération à l'autre (préconditionneur GMRES)
    assert s._coupled_lin["factorizations"] < s.iterations


def test_hagen_poiseuille_axisymmetric():
    nu, f = 0.1, 1.0
    m = rectangle_mesh(0, 1, 0, 1, 2, 24, names={"left": "inlet", "right": "outlet",
                                                 "bottom": "axis", "top": "wall"},
                       types={"top": "wall"}, periodic=[("inlet", "outlet")])
    s = Solver2D(m, nu, {"axis": {"type": "axis"}, "wall": {"type": "wall"}},
                 body_force=(f, 0.0), axisymmetric=True,
                 settings=Settings(algorithm="coupled", relax_U=1.0))
    assert s.run_steady(max_iter=100, tol=1e-10)
    r = m.cell_centers[:, 1]
    ex = f * (1 - r ** 2) / (4 * nu)
    assert np.max(np.abs(s.U[:, 0] - ex)) / ex.max() < 5e-3


def _cavity_thermal(n, Ra, names=None, **kw):
    names = names or {"left": "hot", "right": "cold", "bottom": "adiab", "top": "adiab"}
    m = rectangle_mesh(0, 1, 0, 1, n, n, names=names)
    walls = {k: {"type": "wall"} for k in set(names.values())}
    walls["hot"] = {"type": "wall", "T": 0.5}
    walls["cold"] = {"type": "wall", "T": -0.5}
    return Solver2D(m, 0.71, walls, energy={"Pr": 0.71, "beta": Ra * 0.71, "gravity": (0.0, -1.0)},
                    reference_velocity=20.0, settings=Settings(algorithm="coupled", **kw))


def test_natural_convection_de_vahl_davis_ra1e4():
    s = _cavity_thermal(32, 1e4, relax_T=1.0)
    assert s.settings.relax_U == 0.9                 # auto : flottabilité → 0.9
    assert s.run_steady(max_iter=500, tol=1e-6)
    _, q = s.wall_heat_flux("hot")
    Nu = float(np.sum(q * s.fvm.magSb[s.patch_slices["hot"]]))
    assert Nu == pytest.approx(2.243, rel=0.015)     # de Vahl Davis (1983)


def test_stable_stratification_stays_at_rest():
    s = _cavity_thermal(20, 1e5, names={"left": "side", "right": "side", "bottom": "cold",
                                        "top": "hot"})
    s.run_steady(max_iter=300, tol=1e-9)
    assert np.abs(s.U).max() < 1e-6
