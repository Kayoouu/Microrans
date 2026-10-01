"""Solveur incompressible en 3D (même Solver2D, 3 composantes) : solution exacte de la
conduite carrée, écoulement ABC instationnaire exact, canal turbulent périodique, équivalence
exacte avec le 2D sur un cas extrudé, options refusées."""
import numpy as np
import pytest

from microrans.cases import run_rans_channel
from microrans.fv2d.solver import Settings, Solver2D
from microrans.mesh2d import Circle, cavity_mesh, channel_mesh, o_grid
from microrans.mesh3d import box_mesh, extrude


def _duct_exact(y, z, a, f, nu, nmax=101):
    """Conduite carrée |y|, |z| ≤ a, établie (White, « Viscous Fluid Flow », éq. 3-48)."""
    u = np.zeros_like(y)
    for n in range(1, nmax, 2):
        k = n * np.pi / (2 * a)
        u += (-1) ** ((n - 1) // 2) / n ** 3 * (1 - np.cosh(k * z) / np.cosh(k * a)) * np.cos(
            k * y)
    return 16 * a * a * f / (nu * np.pi ** 3) * u


def _duct(n, a=0.5, nu=0.01, f=1.0):
    m = box_mesh(0, 0.5, -a, a, -a, a, 2, n, n, names={"left": "in", "right": "out"},
                 periodic=[("in", "out")])
    bcs = {k: {"type": "wall"} for k in ("bottom", "top", "back", "front")}
    S = Solver2D(m, nu, bcs, body_force=(f, 0.0, 0.0), settings=Settings(max_iter=2000,
                                                                          tol=1e-8))
    assert S.run_steady()
    C, V = m.cell_centers, m.cell_volumes
    ue = _duct_exact(C[:, 1], C[:, 2], a, f, nu)
    err = np.sqrt(np.sum((S.U[:, 0] - ue) ** 2 * V) / V.sum()) / ue.max()
    return S, err


def test_square_duct_matches_exact_series_second_order():
    S8, e8 = _duct(8)
    S16, e16 = _duct(16)
    assert e16 < 6e-3                                  # mesuré : 2.0e-2, 5.0e-3, 1.3e-3 (32²)
    assert np.log2(e8 / e16) > 1.9
    assert np.abs(S16.U[:, 1:]).max() < 1e-12           # pas d'écoulement secondaire
    # débit : U_moyen exact = (4 a⁴ f / 3ν)(1 − (192/π⁵) Σ tanh(nπ/2)/n⁵) / (4 a²)
    s = sum(np.tanh(n * np.pi / 2) / n ** 5 for n in range(1, 200, 2))
    u_mean = 4 * 0.5 ** 4 / (3 * 0.01) * (1 - 192 / np.pi ** 5 * s) / (4 * 0.5 ** 2)
    # débit trop fort de 5.9 %, 1.5 %, 0.38 % (8², 16², 32²) : ordre 2
    d = [np.sum(S.U[:, 0] * S.mesh.cell_volumes) / S.mesh.cell_volumes.sum() / u_mean - 1
         for S in (S8, S16)]
    assert 0 < d[1] < 0.016 and d[0] / d[1] > 3.8


def _abc(n, scheme, nu=0.1, dt=0.05, t_end=1.0):
    """Écoulement ABC (Arnold-Beltrami-Childress) dans le cube périodique [0, 2π]³ : solution
    exacte de Navier-Stokes 3D, u = u₀ e^{−νt}, p = −|u|²/2 (champ de Beltrami, ω = u)."""
    L = 2 * np.pi
    m = box_mesh(0, L, 0, L, 0, L, n, n, n,
                 periodic=[("left", "right"), ("bottom", "top"), ("back", "front")])

    def u0(C):
        x, y, z = C.T
        return np.column_stack([np.sin(z) + np.cos(y), np.sin(x) + np.cos(z),
                                np.sin(y) + np.cos(x)])

    S = Solver2D(m, nu, {}, settings=Settings(time_scheme=scheme))
    S.U = u0(m.cell_centers)
    S.F_i = np.sum(S.fvm.interp(S.U) * S.fvm.Si, axis=1)
    S.run_transient(dt, t_end)
    ue = u0(m.cell_centers) * np.exp(-nu * t_end)
    eu = np.sqrt(np.mean(np.sum((S.U - ue) ** 2, axis=1)) / np.mean(np.sum(ue ** 2, axis=1)))
    pe = -0.5 * np.sum(ue ** 2, axis=1)
    ep = np.sqrt(np.mean(((S.p - S.p.mean()) - (pe - pe.mean())) ** 2)) / np.ptp(pe)
    return eu, ep


def test_abc_flow_pimple_matches_exact_solution():
    """PIMPLE 3D instationnaire : erreurs mesurées (backward, Δt = 0.05, t = 1) U 5.8e-2,
    6.8e-3, 5.7e-4 et p 2.6e-2, 7.7e-3, 2.5e-3 sur 8³, 16³, 32³ ; Δt = 0.025 ne change pas
    16³ (erreur spatiale dominante)."""
    e8, _ = _abc(8, "backward")
    e16, p16 = _abc(16, "backward")
    assert e16 < 8e-3 and p16 < 1e-2
    assert e8 / e16 > 4.0
    e16_rk3, _ = _abc(16, "rk3")                        # mesuré : 5.2e-3
    assert e16_rk3 < 7e-3


@pytest.mark.parametrize("model", ["sa", "sst_gamma"])
def test_periodic_turbulent_channel_3d_identical_to_2d(model):
    """Canal Re_τ = 395 extrudé sur 2 couches, périodique en x et en z : même solution que le
    2D (sst_gamma passe par le gradient normal 3D de la transition) ; SA comparé au 1D."""
    re_tau = 395.0
    m2 = channel_mesh(1.0, 2.0, 2, 96, first_height=0.4 / re_tau)
    m3 = extrude(m2, 0.0, 0.5, 2, names={"back": "z0", "front": "z1"},
                 periodic=[("z0", "z1")])
    walls = {"bottom": {"type": "wall"}, "top": {"type": "wall"}}
    ub = []
    for m in (m2, m3):
        z = (0.0,) * (m.dim - 2)
        S = Solver2D(m, 1 / re_tau, walls, model=model, body_force=(1.0, 0.0) + z,
                     initial_U=(15.0, 0.0) + z,
                     turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0})
        assert S.run_steady(max_iter=4000, tol=1e-6)
        ub.append(np.sum(S.U[:, 0] * m.cell_volumes) / np.sum(m.cell_volumes))
    assert ub[1] == pytest.approx(ub[0], rel=1e-5)      # mesuré : 1.9e-7 (sa), 3.0e-6
    assert np.abs(S.U[:, 2]).max() < 1e-12
    for p in walls:
        assert S.wall_shear(p)[1].mean() == pytest.approx(1.0, abs=2e-3)
    if model == "sa":
        r1 = run_rans_channel("sa", re_tau, n_cells=256, y1_plus=0.2).summary
        assert ub[1] == pytest.approx(r1["Ub_plus"] * r1["u_tau"], rel=1e-3)


def test_turbulent_square_duct_converges():
    """Conduite carrée turbulente établie (4 parois, coins), SST, Re_τ = 180 : converge, pas
    d'écoulement secondaire (attendu : un modèle à viscosité turbulente linéaire ne prédit pas
    les écoulements secondaires de 2e espèce), coefficient de perte de charge proche de Blasius
    (contrôle d'ordre de grandeur seulement : mesuré λ = 0.0419 contre 0.0410)."""
    from microrans.mesh2d.blocks import grading_for_first_cell
    re_tau, n = 180.0, 32
    r = grading_for_first_cell(n // 2, 0.5 / re_tau)
    g = [(0.5, 0.5, r), (0.5, 0.5, 1 / r)]
    m = box_mesh(0, 0.5, -1, 1, -1, 1, 2, n, n, grading=(1.0, g, g),
                 names={"left": "in", "right": "out"}, periodic=[("in", "out")])
    bcs = {k: {"type": "wall"} for k in ("bottom", "top", "back", "front")}
    S = Solver2D(m, 1 / re_tau, bcs, model="sst", body_force=(1.0, 0.0, 0.0),
                 initial_U=(15.0, 0.0, 0.0),
                 turbulence_inflow={"intensity": 0.05, "viscosity_ratio": 50.0})
    assert S.run_steady(max_iter=4000, tol=1e-6)
    ub = np.sum(S.U[:, 0] * m.cell_volumes) / m.cell_volumes.sum()
    lam, re = 8 * 0.5 / ub ** 2, ub * 2.0 * re_tau     # τ_w moyen = f·A/P = 0.5
    assert lam == pytest.approx(0.316 * re ** -0.25, rel=0.1)
    assert np.abs(S.U[:, 1:]).max() < 1e-10


def _pair(m2, bcs2, nu):
    """Même cas en 2D et extrudé sur une couche, faces z en symétrie."""
    out = []
    for m in (m2, extrude(m2, 0.0, 0.1, 1, types={"back": "symmetry", "front": "symmetry"})):
        bcs = {k: dict(v) for k, v in bcs2.items()}
        if m.dim == 3:
            for v in bcs.values():
                if "U" in v:
                    v["U"] = list(v["U"]) + [0.0]
            bcs["back"] = bcs["front"] = {"type": "symmetry"}
        S = Solver2D(m, nu, bcs, settings=Settings(max_iter=3000, tol=1e-7))
        S.run_steady()
        out.append(S)
    return out


def test_extruded_cavity_identical_to_2d():
    a, b = _pair(cavity_mesh(16), {"lid": {"type": "wall", "U": [1.0, 0.0]},
                                   "walls": {"type": "wall"}}, 0.01)
    assert a.iterations == b.iterations
    assert np.abs(a.U - b.U[:, :2]).max() < 1e-12 and np.abs(b.U[:, 2]).max() < 1e-14
    assert np.abs((a.p - a.p.mean()) - (b.p - b.p.mean())).max() < 1e-12


def test_extruded_thin_cylinder_identical_to_2d():
    """Cellules minces (épaisseur 0.1, mailles jusqu'à ~2.6 dans le plan) contre les plans
    de symétrie : divergeait à la 4e itération tant que la diffusion normale de ces faces
    entrait dans la diagonale commune (voir Solver2D._common_diagonal)."""
    m2 = o_grid(Circle((0, 0), 0.5), 32, 16, 20.0, 0.02, wall_name="cyl")
    a, b = _pair(m2, {"cyl": {"type": "wall"}, "farfield": {"type": "farfield",
                                                             "U": [1.0, 0.0]}}, 0.05)
    assert np.abs(a.U - b.U[:, :2]).max() < 1e-10
    fa = a.forces(["cyl"])["cyl"]["total"]
    fb = b.forces(["cyl"])["cyl"]["total"]
    assert fb[0] / 0.1 == pytest.approx(fa[0], rel=1e-10)
    assert abs(fb[2]) < 1e-12


def test_options_unavailable_in_3d_refused():
    m = box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2, types={k: "wall" for k in
                                                   ("left", "right", "bottom", "top",
                                                    "back", "front")})
    bcs = {k: {"type": "wall"} for k in ("left", "right", "bottom", "top", "back", "front")}
    with pytest.raises(ValueError, match="zones poreuses"):
        Solver2D(m, 0.01, bcs, porous=[{"region": "rectangle"}])
    with pytest.raises(ValueError, match="solveur couplé"):
        Solver2D(m, 0.01, bcs, settings=Settings(algorithm="coupled"))
    with pytest.raises(ValueError, match="axisymétrique disponible"):
        Solver2D(m, 0.01, bcs, axisymmetric=True)
    bad = {**bcs, "top": {"type": "wall", "U": [1.0, 0.0]}}
    with pytest.raises(ValueError, match=r"\[boundary.top\] U = \[1.0, 0.0\] : 3 composantes"):
        Solver2D(m, 0.01, bad)
    with pytest.raises(ValueError, match="body_force"):
        Solver2D(m, 0.01, bcs, body_force=(1.0, 0.0))
