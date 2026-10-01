"""Solveur compressible 2D (compressible.py) : solutions exactes, conservation, conditions
aux limites, schéma implicite, cas de calcul. Tests rapides (maillages grossiers)."""
import json

import numpy as np
import pytest

from microrans.fv2d.compressible import (CompressibleSettings, CompressibleSolver2D, Gas,
                                         euler_flux, hllc_flux, make_state, roe_flux)
from microrans.fv2d.gasdynamics import oblique_shock, riemann_exact, riemann_star
from microrans.mesh2d.blocks import block_mesh, rectangle_mesh


# ---------------------------------------------------------------- références exactes
def test_riemann_exact_sod_matches_toro():
    # Toro 2009, tableau 4.3 (test 1) : p* = 0.30313, u* = 0.92745
    p, u = riemann_star((1.0, 0.0, 1.0), (0.125, 0.0, 0.1))
    assert p == pytest.approx(0.30313, abs=1e-5)
    assert u == pytest.approx(0.92745, abs=1e-5)
    r, _, _, w = riemann_exact((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), np.array([0.0]))
    assert w["rho_star_left"] == pytest.approx(0.42632, abs=1e-5)
    assert w["rho_star_right"] == pytest.approx(0.26557, abs=1e-5)
    assert w["right_shock"] == pytest.approx(1.75216, abs=1e-4)


def test_oblique_shock_theory():
    # NACA Report 1135 : M = 2, θ = 10° → β = 39.31°, p2/p1 = 1.7066, M2 = 1.6405
    d = oblique_shock(2.0, 10.0)
    assert d["beta_deg"] == pytest.approx(39.314, abs=1e-3)
    assert d["p_ratio"] == pytest.approx(1.7066, abs=1e-4)
    assert d["M2"] == pytest.approx(1.6405, abs=1e-4)
    with pytest.raises(ValueError):
        oblique_shock(2.0, 30.0)                  # au-delà de θ_max ≈ 22.97° : détaché


# ---------------------------------------------------------------- flux
def test_fluxes_consistent_and_jacobians():
    from microrans.fv2d.compressible_implicit import flux_jacobian, roe_dissipation
    g = 1.4
    rng = np.random.default_rng(0)
    n = 20
    WL = np.vstack([rng.uniform(0.5, 2, n), rng.uniform(-300, 300, n),
                    rng.uniform(-300, 300, n), rng.uniform(5e4, 2e5, n)])
    WR = WL * rng.uniform(0.8, 1.2, (4, n))
    a = rng.uniform(0, 2 * np.pi, n)
    nx, ny = np.cos(a), np.sin(a)
    F = np.array(euler_flux(*WL, nx, ny, g))
    for f in (roe_flux, hllc_flux):          # consistance : F(W, W) = F(W)
        assert np.allclose(np.array(f(WL, WL, nx, ny, g)), F, rtol=1e-12, atol=1e-9)
    # dissipation de Roe : |Ã| (Q_R − Q_L) = F_L + F_R − 2 F_roe (propriété de Roe)

    def cons(W):
        r, u, v, p = W
        return np.vstack([r, r * u, r * v, p / (g - 1) + 0.5 * r * (u * u + v * v)])
    D = roe_dissipation(WL, WR, nx, ny, g, 0.1)
    FR = np.array(euler_flux(*WR, nx, ny, g))
    Froe = np.array(roe_flux(WL, WR, nx, ny, g, 0.1))
    lhs = np.einsum("nij,jn->in", D, cons(WR) - cons(WL))
    assert np.allclose(lhs, F + FR - 2 * Froe, rtol=1e-9, atol=1e-6 * np.abs(F).max())
    # jacobienne du flux : différences finies
    A = flux_jacobian(WL[:, :1], nx[:1], ny[:1], g)[0]
    Q = cons(WL[:, :1])[:, 0]

    def flux_of(Qv):
        r = Qv[0]
        u, v = Qv[1] / r, Qv[2] / r
        p = (g - 1) * (Qv[3] - 0.5 * r * (u * u + v * v))
        return np.array(euler_flux(r, u, v, p, nx[0], ny[0], g))
    J = np.column_stack([(flux_of(Q + h * e) - flux_of(Q - h * e)) / (2 * h)
                         for e, h in zip(np.eye(4), 1e-6 * np.abs(Q))])
    assert np.allclose(J, A, rtol=1e-6, atol=1e-6 * np.abs(A).max())


# ---------------------------------------------------------------- tube à choc
def _sod(nx, **kw):
    mesh = rectangle_mesh(0, 1, 0, 1.0 / nx, nx, 1,
                          names={"left": "left", "right": "right", "bottom": "w", "top": "w"})
    gas = Gas()
    fs = make_state(gas, 0.0, 1e5, density=1.0)
    bcs = {"left": {"type": "supersonic_outlet"}, "right": {"type": "supersonic_outlet"},
           "w": {"type": "symmetry"}}
    init = lambda x, y: (np.where(x < 0.5, 1.0, 0.125), 0.0, 0.0,  # noqa: E731
                         np.where(x < 0.5, 1e5, 1e4))
    S = CompressibleSolver2D(mesh, gas, fs, bcs, CompressibleSettings(cfl=0.8, **kw), init)
    t_end = 0.2 / np.sqrt(1e5)
    S.run_transient(t_end)
    x = mesh.cell_centers[:, 0]
    r, u, p, w = riemann_exact((1.0, 0, 1e5), (0.125, 0, 1e4), (x - 0.5) / t_end)
    return S, x, r, u, p, w


@pytest.mark.parametrize("flux", ["roe", "hllc"])
def test_sod_shock_tube_against_exact_solution(flux):
    S, x, r, u, p, w = _sod(100, flux=flux)
    W = S.W
    err = np.mean(np.abs(W[:, 0] - r))
    # mesuré : 6.0e-3 (Roe), 5.6e-3 (HLLC) ; ordre 1 : 2.0e-2
    assert err < 7.5e-3
    assert np.mean(np.abs(W[:, 3] - p)) / 1e5 < 6e-3
    # position du choc (ρ à mi-chemin entre ρ*_R et ρ_R) à ± 1.5 maille
    rho_mid = 0.5 * (w["rho_star_right"] + 0.125)
    k = np.nonzero((W[:-1, 0] > rho_mid) & (W[1:, 0] <= rho_mid) & (x[:-1] > 0.75))[0][0]
    xs = x[k] + (W[k, 0] - rho_mid) / (W[k, 0] - W[k + 1, 0]) * (x[k + 1] - x[k])
    x_exact = 0.5 + w["right_shock"] * 0.2 / np.sqrt(1e5)
    assert abs(xs - x_exact) < 0.015
    # plateau étoile (entre contact et choc)
    plateau = (x > 0.72) & (x < 0.82)
    assert np.allclose(W[plateau, 3], w["p_star"], rtol=0.01)


def test_sod_first_order_is_less_accurate():
    S1, x, r, *_ = _sod(100, order=1)
    S2, *_ = _sod(100, order=2)
    e1 = np.mean(np.abs(S1.W[:, 0] - r))
    e2 = np.mean(np.abs(S2.W[:, 0] - r))
    assert e2 < 0.5 * e1


# ---------------------------------------------------------------- conservation
def test_conservation_closed_box_and_periodic():
    # boîte fermée (parois glissantes) : masse et énergie conservées au bit près ;
    # boîte doublement périodique : quantité de mouvement aussi
    gas = Gas()
    fs = make_state(gas, 0.3, 1e5, 300.0)

    def bump(x, y):
        r2 = (x - 0.4) ** 2 + (y - 0.5) ** 2
        return (1.2 * (1 + 0.3 * np.exp(-r2 / 0.02)), 30.0, -20.0,
                1e5 * (1 + 0.5 * np.exp(-r2 / 0.02)))
    m = rectangle_mesh(0, 1, 0, 1, 16, 16, names={"left": "w", "right": "w", "bottom": "w",
                                                   "top": "w"})
    S = CompressibleSolver2D(m, gas, fs, {"w": {"type": "slip_wall"}},
                             CompressibleSettings(cfl=0.8), bump)
    t0 = S.totals()
    S.run_transient(1e-3)
    t1 = S.totals()
    assert S.iterations > 20
    assert abs(t1["mass"] / t0["mass"] - 1) < 1e-13
    assert abs(t1["energy"] / t0["energy"] - 1) < 1e-13
    mp = rectangle_mesh(0, 1, 0, 1, 16, 16, names={"left": "l", "right": "r", "bottom": "b",
                                                    "top": "t"},
                        periodic=[("l", "r"), ("b", "t")])
    S = CompressibleSolver2D(mp, gas, fs, {}, CompressibleSettings(cfl=0.8, flux="hllc"),
                             bump)
    t0 = S.totals()
    S.run_transient(1e-3)
    t1 = S.totals()
    for k in t0:
        assert abs(t1[k] - t0[k]) <= 1e-12 * abs(t0[k]) + 1e-12 * t0["mass"] * 400.0


def test_freestream_preserved_on_triangles():
    from microrans.mesh2d.geometry import Rectangle
    from microrans.mesh2d.unstructured import triangulate
    dom = Rectangle(0, 0, 2, 1, names={"left": "in", "right": "out", "bottom": "b",
                                       "top": "t"})
    mesh = triangulate(dom, 0.15)
    gas = Gas()
    fs = make_state(gas, 0.7, 8e4, 250.0, angle_deg=20.0)
    S = CompressibleSolver2D(mesh, gas, fs, {n: {"type": "farfield"} for n in
                                             ("in", "out", "b", "t")})
    R = S.residual(S.Q)
    scale = np.abs(np.array(euler_flux(fs.rho, fs.u, fs.v, fs.p, 1.0, 0.0, 1.4))).max()
    assert np.abs(R).max() < 1e-10 * scale


# ---------------------------------------------------------------- conditions aux limites
def test_boundary_ghost_states():
    mesh = rectangle_mesh(0, 1, 0, 1, 4, 4, names={"left": "in", "right": "out",
                                                    "bottom": "wall", "top": "far"})
    gas = Gas(viscosity="constant", mu=1e-3)
    fs = make_state(gas, 0.5, 1e5, 300.0)
    bcs = {"in": {"type": "inlet"}, "out": {"type": "outlet", "p": 9e4},
           "wall": {"type": "wall", "T": 350.0}, "far": {"type": "farfield"}}
    S = CompressibleSolver2D(mesh, gas, fs, bcs)
    Wi = np.tile(np.array([fs.rho, fs.u, fs.v, fs.p])[:, None], (1, S.nb))
    G = S.ghost(Wi)
    g = gas.gamma
    sl = S.patch_slices
    # champ lointain avec l'état amont à l'intérieur : état fantôme = état amont
    assert np.allclose(G[:, sl["far"]], Wi[:, sl["far"]], rtol=1e-12)
    # entrée subsonique : p0, T0 amont retrouvées, écoulement dans la direction imposée
    r, u, v, p = G[:, sl["in"]]
    T = p / (r * gas.R)
    M2 = (u * u + v * v) / (g * p / r)
    assert np.allclose(p * (1 + 0.5 * (g - 1) * M2) ** (g / (g - 1)), fs.p0, rtol=1e-12)
    assert np.allclose(T + 0.5 * (u * u + v * v) / gas.cp, fs.T0, rtol=1e-12)
    assert np.allclose(u, fs.speed * 1.0, rtol=1e-9) and np.allclose(v, 0.0, atol=1e-9)
    # sortie subsonique : pression imposée
    assert np.allclose(G[3, sl["out"]], 9e4)
    # paroi adhérente : vitesse moyenne nulle à la face, T_w imposée pour les gradients
    Wb = S.boundary_values(np.tile(np.array([fs.rho, fs.u, fs.v, fs.p])[:, None],
                                   (1, S.nc)))
    assert np.allclose(Wb[1:3, sl["wall"]], 0.0)
    assert np.allclose(Wb[4, sl["wall"]], 350.0)
    # paroi glissante (Euler) : flux de masse nul, quantité de mouvement normale
    S2 = CompressibleSolver2D(mesh, Gas(), fs, {**bcs, "wall": {"type": "slip_wall"}})
    Wf = Wi[:, sl["wall"]]
    nb = S2.nb_hat[sl["wall"]].T
    for f in (roe_flux, hllc_flux):
        F = np.array(f(Wf, S2.ghost(Wi)[:, sl["wall"]], nb[0], nb[1], g))
        assert np.allclose(F[0], 0.0, atol=1e-9) and np.allclose(F[3], 0.0, atol=1e-6)
        assert np.allclose(F[1] * nb[1] - F[2] * nb[0], 0.0, atol=1e-6)   # pas tangentiel


# ---------------------------------------------------------------- stationnaire
def _ramp(n, scheme):
    th = np.radians(10.0)
    V = [[0, 0], [0.5, 0], [1.5, np.tan(th)], [1.5, 1.0], [0.5, 1.0], [0, 1.0]]
    blocks = [{"vertices": [0, 1, 4, 5], "cells": [n // 2, n]},
              {"vertices": [1, 2, 3, 4], "cells": [n, n]}]
    patches = {"inlet": {"type": "patch", "faces": [[5, 0]]},
               "outlet": {"type": "patch", "faces": [[2, 3]]},
               "top": {"type": "patch", "faces": [[3, 4], [4, 5]]},
               "ramp": {"type": "wall", "faces": [[0, 1], [1, 2]]}}
    mesh = block_mesh(V, blocks, patches=patches)
    gas = Gas()
    fs = make_state(gas, 2.0, 101325.0, 288.15)
    bcs = {"inlet": {"type": "supersonic_inlet"}, "outlet": {"type": "supersonic_outlet"},
           "top": {"type": "farfield"}, "ramp": {"type": "wall"}}
    S = CompressibleSolver2D(mesh, gas, fs, bcs, CompressibleSettings(
        cfl=5.0 if scheme == "implicit" else 2.0, steady_scheme=scheme, max_iter=2000,
        tol=1e-8))
    assert S.run_steady()
    return S, fs


def test_supersonic_ramp_oblique_shock():
    S, fs = _ramp(24, "implicit")
    ref = oblique_shock(2.0, 10.0)
    C = S.mesh.cell_centers
    p = S.p / fs.p
    # zone entre la rampe et le choc, loin du coin et des frontières
    xi, yi = C[:, 0], C[:, 1]
    yr = (xi - 0.5) * np.tan(np.radians(10.0))
    ys = (xi - 0.5) * np.tan(np.radians(ref["beta_deg"]))
    zone = (xi > 0.9) & (xi < 1.45) & (yi > yr + 0.05) & (yi < ys - 0.12)
    assert zone.sum() > 10
    assert np.mean(p[zone]) == pytest.approx(ref["p_ratio"], rel=0.01)     # mesuré : 0.1 %
    assert np.mean(S.mach[zone]) == pytest.approx(ref["M2"], rel=0.01)
    # angle du choc : abscisse du saut de pression (mi-hauteur) sur deux lignes y = cte
    from scipy.interpolate import LinearNDInterpolator
    f = LinearNDInterpolator(C, p)
    pm = 0.5 * (1 + ref["p_ratio"])
    xs = []
    for y in (0.3, 0.6):
        x = np.linspace(0.6, 1.45, 800)
        pv = f(x, np.full_like(x, y))
        k = np.nonzero((pv[:-1] < pm) & (pv[1:] >= pm))[0][0]
        xs.append(x[k] + (pm - pv[k]) / (pv[k + 1] - pv[k]) * (x[k + 1] - x[k]))
    beta = np.degrees(np.arctan(0.3 / (xs[1] - xs[0])))
    assert beta == pytest.approx(ref["beta_deg"], abs=1.0)
    # effort sur la rampe : (p2 − p1) × projection de la rampe
    F = S.forces(["ramp"])["ramp"]["total"]
    assert F[0] == pytest.approx((ref["p_ratio"] - 1) * fs.p * np.tan(np.radians(10)),
                                 rel=0.01)


def test_explicit_and_implicit_reach_same_solution():
    Si, _ = _ramp(12, "implicit")
    Se, _ = _ramp(12, "rk3")
    assert Si.iterations < Se.iterations
    assert np.allclose(Si.W, Se.W, rtol=1e-5, atol=1e-5 * np.abs(Se.W).max(axis=0))


def test_couette_flow_viscous_heating():
    # Couette compressible, μ constant : u = U y/h, T = T_w + Pr U²/(2 c_p) η(1 − η)
    # (White, « Viscous Fluid Flow », 3e éd., § 3-2 : Couette avec dissipation visqueuse)
    h, U, Tw = 1e-3, 200.0, 300.0
    mesh = rectangle_mesh(0, h / 6, 0, h, 4, 24,
                          names={"left": "l", "right": "r", "bottom": "b", "top": "t"},
                          periodic=[("l", "r")])
    gas = Gas(viscosity="constant", mu=1e-2)
    fs = make_state(gas, 0.0, 101325.0, Tw)
    S = CompressibleSolver2D(mesh, gas, fs, {"b": {"type": "wall", "T": Tw},
                                             "t": {"type": "wall", "T": Tw, "U": [U, 0]}},
                             CompressibleSettings(steady_scheme="implicit", cfl=5,
                                                  max_iter=600, tol=1e-9))
    assert S.run_steady()
    eta = mesh.cell_centers[:, 1] / h
    dT = gas.Pr * U ** 2 / (2 * gas.cp)
    assert np.abs(S.U[:, 0] - U * eta).max() < 1e-6 * U
    assert np.abs(S.T - (Tw + dT * eta * (1 - eta))).max() < 0.005 * dT / 4
    _, tau, _ = S.wall_shear("b")
    assert np.allclose(tau, gas.mu * U / h, rtol=1e-5)
    _, q = S.wall_heat_flux("b")
    assert np.allclose(q, -gas.mu * U ** 2 / (2 * h), rtol=1e-4)


# ---------------------------------------------------------------- cas de calcul
def test_case_file_run_and_restart(tmp_path):
    from microrans.cli import examples_dir
    from microrans.fv2d.case import run_case
    from microrans.mesh2d.builder import load_config
    cfg = load_config(examples_dir() / "compressible_rampe_mach2.toml")
    for b in cfg["mesh"]["blocks"]:
        b["cells"] = [max(c // 4, 4) for c in b["cells"]]
    cfg["solver"].update(max_iter=60, steady_scheme="implicit", cfl=5.0)
    out = tmp_path / "rampe"
    s1 = run_case(cfg, base_dir=examples_dir(), out_dir=out, verbose=False, plot=True)
    for f in ("summary.json", "fields.vtk", "history.csv", "Mach.png", "p.png",
              "wall_ramp.csv", "checkpoint.npz", "convergence.png"):
        assert (out / f).is_file(), f
    d = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert d["solver"] == "compressible" and d["freestream"]["mach"] == pytest.approx(2.0)
    assert d["ramp"]["Cd"] > 0 and d["ramp"]["Cl"] < 0
    # reprise exacte : la numérotation des itérations continue
    cfg2 = {**cfg, "initial": {"restart": str(out / "checkpoint.npz")}}
    s2 = run_case(cfg2, base_dir=examples_dir(), out_dir=tmp_path / "suite", verbose=False,
                  plot=False)
    assert s2["restart"]["mode"] == "exact"
    assert s2["iterations"] == s1["iterations"] + s2["iterations_this_run"]


def test_case_file_transient_sod(tmp_path):
    from microrans.cli import examples_dir
    from microrans.fv2d.case import run_case
    from microrans.mesh2d.builder import load_config
    cfg = load_config(examples_dir() / "compressible_tube_sod.toml")
    cfg["mesh"].update(nx=80, y1=1.0 / 80)
    cfg["output"]["probes"] = [[0.75, 0.5 / 80]]
    cfg["output"]["lines"] = [{"name": "axe", "start": [0, 0.5 / 80], "end": [1, 0.5 / 80],
                               "n": 50}]
    s = run_case(cfg, base_dir=examples_dir(), out_dir=tmp_path, verbose=False, plot=False)
    assert s["mode"] == "transient" and s["time"] == pytest.approx(6.32455532e-4)
    # bords ouverts (extrapolation) : quelques 1e-9 de la masse sortent par les extrémités
    assert s["totals_final"]["mass"] == pytest.approx(s["totals_initial"]["mass"], rel=1e-7)
    # sonde dans le plateau étoile droit : ρ*_R = 0.26557
    assert s["probes"][0]["rho"] == pytest.approx(0.26557, rel=0.03)
    assert (tmp_path / "line_axe.csv").is_file()


def test_bad_boundary_type_and_missing_patch():
    mesh = rectangle_mesh(0, 1, 0, 1, 2, 2)
    gas = Gas()
    fs = make_state(gas, 0.5, 1e5, 300.0)
    with pytest.raises(ValueError, match="inconnue"):
        CompressibleSolver2D(mesh, gas, fs, {n: {"type": "velocity_inlet"} for n in
                                             ("left", "right", "bottom", "top")})
    with pytest.raises(ValueError, match="manquantes"):
        CompressibleSolver2D(mesh, gas, fs, {"left": {"type": "farfield"}})


def test_farfield_vortex_correction_removes_domain_size_effect(tmp_path):
    """NACA 0012, M = 0.5, α = 1.25°, maillage grossier : sans correction, C_l dépend de la
    distance du champ lointain (mesuré : 0.1590 à 10 cordes, 0.1651 à 100) ; avec le
    tourbillon ponctuel, à 0.2 % près (0.1655 / 0.1658)."""
    import math

    from microrans.cli import examples_dir
    from microrans.fv2d.case import run_case
    from microrans.mesh2d.builder import load_config
    cl = {}
    for R in (10.0, 100.0):
        for vortex in (False, True):
            cfg = load_config(examples_dir() / "compressible_naca0012_transsonique.toml")
            nr = int(round(24 * math.log(R / 0.002) / math.log(30 / 0.002)))
            cfg["mesh"].update(n_around=64, n_radial=nr, farfield_radius=R, first_height=6e-3)
            cfg["flow"]["mach"] = 0.5
            cfg["solver"].update(max_iter=600, tol=1e-7, monitor_tol=1e-5)
            if vortex:
                cfg["boundary"]["farfield"]["vortex"] = [0.25, 0.0]
            else:
                cfg["boundary"]["farfield"].pop("vortex", None)
            cfg["output"] = {"vtk": False, "checkpoint": False, "forces": ["airfoil"]}
            s = run_case(cfg, examples_dir(), out_dir=tmp_path, verbose=False, plot=False)
            assert s["converged"]
            cl[R, vortex] = s["airfoil"]["Cl"]
    assert cl[100.0, True] == pytest.approx(cl[10.0, True], rel=5e-3)
    assert cl[100.0, False] > 1.02 * cl[10.0, False]


def test_implicit_cfl_cap_kept_during_transonic_startup():
    """NACA 0012, M = 0.8, 192 × 64 : au démarrage les résidus stagnent pendant que le choc
    se déplace. L'ancienne règle divisait alors le plafond de CFL par 2 (dès l'itération
    ~60, jusqu'à 3 fois, sans retour : calcul non convergé en 3000 itérations) ; il n'est
    réduit que si la solution oscille (cycle limite), pas quand elle dérive."""
    from microrans.cli import examples_dir
    from microrans.fv2d.compressible_case import build_compressible_solver
    from microrans.mesh2d.builder import load_config
    cfg = load_config(examples_dir() / "compressible_naca0012_transsonique.toml")
    cfg["solver"].update(max_iter=120, tol=1e-12, monitor_tol=None)
    s = build_compressible_solver(cfg, examples_dir())
    s.run_steady()
    assert s.cfl_cuts_done == 0
    assert max(s.history[-1][k] for k in ("rho", "rhoU", "rhoV", "rhoE")) < 1e-2


def test_unphysical_state_gives_measured_tips(tmp_path):
    # audit D2 : rampe Mach 2 en RK3 à cfl = 4 → arrêt à la 9e itération ; l'arrêt
    # propose le schéma implicite (même solution à cfl 8, 6 fois plus rapide)
    import tomllib
    from pathlib import Path

    from microrans.fv2d.compressible_case import run_compressible_case
    ex = Path(__file__).parents[1] / "microrans" / "examples"
    cfg = tomllib.loads((ex / "compressible_rampe_mach2.toml").read_text(encoding="utf-8"))
    cfg["solver"].update(cfl=4.0, max_iter=50)
    with pytest.raises(FloatingPointError) as exc:
        run_compressible_case(cfg, base_dir=ex, out_dir=tmp_path, verbose=False, plot=False)
    assert "État non physique" in str(exc.value)
    assert "steady_scheme = \"implicit\"" in str(exc.value)
