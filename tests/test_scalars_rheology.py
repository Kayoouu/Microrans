"""Scalaires passifs (transport, limiteur, bilans, reprise) et fluides non newtoniens."""
import copy

import numpy as np
import pytest
from scipy.optimize import brentq

from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.case import run_case
from microrans.fv2d.rheology import Rheology
from microrans.mesh2d import channel_mesh, rectangle_mesh

PLUG = {"inlet": {"type": "inlet", "U": [1.0, 0.0]}, "outlet": {"type": "outlet"},
        "sym": {"type": "symmetry"}}


def _plug_mesh(nx):
    return rectangle_mesh(0, 1.0, 0, 0.1, nx, 2, names={"left": "inlet", "right": "outlet",
                                                        "bottom": "sym", "top": "sym"})


def test_scalar_convection_diffusion_source_exact_order2():
    """U c' = D c'' + S, c(0) = 0, c'(L) = 0 (sortie à gradient nul) : solution exacte ;
    ordre 2 et bilan (flux sortants = source) à la précision machine."""
    D, S = 0.02, 1.0
    errs = []
    for nx in (20, 40, 80):
        m = _plug_mesh(nx)
        s = Solver2D(m, 0.01, PLUG, initial_U=(1.0, 0.0),
                     scalars={"c": {"diffusivity": D, "source": S}})
        s.run_steady(max_iter=300, tol=1e-9)
        x = m.cell_centers[:, 0]
        ex = S * x - S * D * (np.exp((x - 1.0) / D) - np.exp(-1.0 / D))
        errs.append(np.sqrt(np.mean((s.scalars["c"] - ex) ** 2)))
        f = s.scalar_fluxes("c")
        assert f["inlet"]["flux"] + f["outlet"]["flux"] + f["sym"]["flux"] == pytest.approx(
            f["source_total"], rel=1e-10)
    order = np.log2(np.array(errs[:-1]) / errs[1:])
    assert np.all(order > 1.9), order


def test_scalar_equals_temperature_with_same_diffusivity():
    """Même équation que T (D = ν/Pr, mêmes conditions) : mêmes champs."""
    m = rectangle_mesh(0, 4.0, 0, 1.0, 40, 12, names={"left": "inlet", "right": "outlet",
                                                     "bottom": "wall", "top": "wall"})
    bc = {"inlet": {"type": "inlet", "U": [1.0, 0.0], "T": 0.0, "scalars": {"c": 0.0}},
          "outlet": {"type": "outlet"},
          "wall": {"type": "wall", "T": 1.0, "scalars": {"c": 1.0}}}
    s = Solver2D(m, 0.02, bc, energy={"Pr": 2.0}, settings=Settings(relax_T=1.0),
                 scalars={"c": {"schmidt": 2.0, "scheme": "linearUpwind"}})
    s.run_steady(max_iter=2000, tol=1e-10)
    assert np.allclose(s.scalars["c"], s.T, rtol=0, atol=1e-8)


def test_scalar_mass_conserved_and_bounded_in_closed_cavity():
    """Cavité fermée, créneau de concentration : masse conservée exactement ; le limiteur
    réduit les dépassements de 0 ≤ c ≤ 1 d'un facteur > 100 par rapport à linearUpwind
    (≈ 10⁻¹ → ≈ 2·10⁻⁴ en PIMPLE : correction différée non itérée à convergence)."""
    m = rectangle_mesh(0, 1, 0, 1, 24, 24, names={"left": "w", "right": "w", "bottom": "w",
                                                   "top": "lid"})
    bc = {"lid": {"type": "wall", "U": [1.0, 0.0]}, "w": {"type": "wall"}}
    blob = "where(abs(x - 0.5) < 0.2, where(abs(y - 0.7) < 0.15, 1.0, 0.0), 0.0)"
    over = {}
    for scheme in ("linearUpwindLimited", "linearUpwind"):
        s = Solver2D(m, 0.01, bc, reference_velocity=1.0,
                     scalars={"c": {"diffusivity": 0.0, "initial": blob, "scheme": scheme}})
        m0 = float(np.sum(s.scalars["c"] * s.fvm.V))
        s.run_transient(0.01, 1.0)
        c = s.scalars["c"]
        assert float(np.sum(c * s.fvm.V)) == pytest.approx(m0, rel=1e-9)
        over[scheme] = max(-c.min(), c.max() - 1.0)
    assert over["linearUpwindLimited"] < 1e-3
    assert over["linearUpwind"] > 100 * over["linearUpwindLimited"]


def test_scalar_case_outputs_and_exact_restart(tmp_path):
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 3, "y0": 0, "y1": 1, "nx": 30, "ny": 10,
                    "names": {"left": "inlet", "right": "outlet", "bottom": "wall",
                              "top": "wall"}},
           "physics": {"nu": 0.02},
           "scalars": {"c": {"diffusivity": 0.01}, "age": {"diffusivity": 0.001,
                                                          "source": 1.0}},
           "boundary": {"inlet": {"type": "inlet", "U": [1.0, 0.0], "scalars": {"c": 1.0}},
                        "outlet": {"type": "outlet"}, "wall": {"type": "wall"}},
           "solver": {"mode": "transient", "dt": 0.05, "t_end": 1.0,
                      "time_scheme": "backward"},
           "output": {"vtk": False, "plots": False, "checkpoint_minutes": 0,
                      "probes": [[1.5, 0.5]], "animate": "c", "animate_every": 5}}
    _, a = run_case(copy.deepcopy(cfg), out_dir=tmp_path / "a", verbose=False, plot=False,
                    return_solver=True)
    c1 = copy.deepcopy(cfg)
    c1["solver"]["t_end"] = 0.6
    run_case(c1, out_dir=tmp_path / "b", verbose=False, plot=False)
    c2 = copy.deepcopy(cfg)
    c2["initial"] = {"restart": str(tmp_path / "b" / "checkpoint.npz")}
    s2, b = run_case(c2, out_dir=tmp_path / "c", verbose=False, plot=False,
                     return_solver=True)
    for k in ("c", "age"):
        assert np.array_equal(a.scalars[k], b.scalars[k])
    assert "probe1_c" in open(tmp_path / "c" / "history.csv", encoding="utf-8").readline()
    assert s2["scalars"]["c"]["inlet"]["flux"] < 0              # entre par l'entrée
    assert s2["animation"].endswith("animation_c.gif")


def test_scalar_name_and_parameter_validation():
    m = _plug_mesh(10)
    with pytest.raises(ValueError, match="réservé"):
        Solver2D(m, 0.01, PLUG, scalars={"p": {"diffusivity": 1.0}})
    with pytest.raises(ValueError, match="diffusivity"):
        Solver2D(m, 0.01, PLUG, scalars={"c": {}})


# ------------------------------------------------------------------ non newtonien
def _channel_run(spec, ny, uref):
    m = channel_mesh(1.0, 2.0, 2, ny)
    s = Solver2D(m, 0.1, {"bottom": {"type": "wall"}, "top": {"type": "wall"}},
                 body_force=(1.0, 0.0), viscosity=spec, reference_velocity=uref,
                 settings=Settings(relax_U=1.0))
    assert s.run_steady(max_iter=3000, tol=1e-10)
    return s, np.abs(m.cell_centers[:, 1] - 1.0)


@pytest.mark.parametrize("n, order", [(0.5, 1.95), (1.5, 1.6)])
def test_power_law_channel_analytic(n, order):
    """Loi puissance, canal plan entraîné par une force G = 1 : profil exact
    u = n/(n+1) (G/K)^(1/n) (h^((n+1)/n) − η^((n+1)/n)). Ordre 2 pour n = 0.5 ; 5/3 pour
    n = 1.5 (dérivée seconde du profil exact infinie sur l'axe)."""
    K = 0.1
    umax = n / (n + 1) * (1.0 / K) ** (1 / n)
    errs = []
    for ny in (16, 32, 64):
        s, eta = _channel_run({"model": "power_law", "K": K, "n": n}, ny, umax)
        ex = umax * (1 - eta ** ((n + 1) / n))
        errs.append(np.max(np.abs(s.U[:, 0] - ex)) / umax)
    assert errs[-1] < 1e-3
    assert np.all(np.log2(np.array(errs[:-1]) / errs[1:]) > order)


def test_carreau_channel_against_integrated_profile():
    """Carreau : profil de référence intégré numériquement (τ = G η, γ̇ = f⁻¹(τ))."""
    spec = {"model": "carreau", "nu0": 1.0, "nu_inf": 0.01, "lambda": 2.0, "n": 0.4}
    rh = Rheology(spec)
    s, eta = _channel_run(spec, 64, 1.0)
    sel = np.argsort(eta)[::4]

    def gam(tau):
        return brentq(lambda g: float(rh(np.array([g]))[0]) * g - tau, 0.0, 1e8, xtol=1e-14)
    ref = []
    for e in eta[sel]:
        xs = np.linspace(e, 1.0, 801)
        ref.append(np.trapezoid([gam(x) for x in xs], xs))
    ref = np.array(ref)
    assert np.max(np.abs(s.U[sel, 0] - ref)) < 1.5e-3 * ref.max()


def test_bingham_plug_flow():
    """Bingham (τ_y = 0.3, K = 0.1) : bouchon |η| < 0.3, vitesse du bouchon 2.45 ;
    régularisation bi-visqueuse (nu_max) + surface d'écoulement → ordre ~1.3."""
    spec = {"model": "bingham", "tau_y": 0.3, "K": 0.1, "nu_max": 100.0}
    s, eta = _channel_run(spec, 128, 2.0)
    ex = np.where(eta < 0.3, 5.0 * 0.49, 5.0 * (1 - eta ** 2) - 3.0 * (1 - eta))
    assert np.max(np.abs(s.U[:, 0] - ex)) < 6e-3 * ex.max()
    plug = s.nu_lam >= 0.999 * 100.0
    assert np.all(eta[plug] < 0.32) and np.all(plug[eta < 0.25])


def test_power_law_pipe_axisymmetric():
    n, K = 0.5, 0.1
    umax = n / (n + 1) * (1.0 / (2 * K)) ** (1 / n)
    m = rectangle_mesh(0, 1.0, 0, 1.0, 2, 32, names={"left": "i", "right": "o",
                                                      "bottom": "axis", "top": "wall"},
                       periodic=[("i", "o")])
    s = Solver2D(m, K, {"axis": {"type": "axis"}, "wall": {"type": "wall"}},
                 body_force=(1.0, 0.0), axisymmetric=True, reference_velocity=umax,
                 viscosity={"model": "power_law", "K": K, "n": n},
                 settings=Settings(relax_U=1.0))
    assert s.run_steady(max_iter=2000, tol=1e-10)
    r = m.cell_centers[:, 1]
    ex = umax * (1 - r ** ((n + 1) / n))
    assert np.max(np.abs(s.U[:, 0] - ex)) < 2e-3 * umax


def test_rheology_models_and_validation():
    g = np.array([1e-3, 1.0, 1e3])
    assert np.allclose(Rheology({"model": "power_law", "K": 2.0, "n": 1.0,
                                 "nu_min": 0, "nu_max": 1e9})(g), 2.0)
    cr = Rheology({"model": "cross", "nu0": 1.0, "nu_inf": 0.1, "m": 1.0, "n": 1.0})
    assert cr(np.array([1.0]))[0] == pytest.approx(0.55)
    cas = Rheology({"model": "casson", "tau_y": 1.0, "nu_inf": 1.0, "nu_max": 1e9})
    assert cas(np.array([1.0]))[0] == pytest.approx(4.0)
    with pytest.raises(ValueError, match="manquants"):
        Rheology({"model": "carreau", "nu0": 1.0})
    m = _plug_mesh(10)
    with pytest.raises(ValueError, match="laminaire"):
        Solver2D(m, 0.01, PLUG, model="sa", viscosity={"model": "power_law", "K": 1, "n": 1})


def test_non_newtonian_case_reports_viscosity(tmp_path):
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 4, "y0": 0, "y1": 1, "nx": 40,
                    "ny": 12, "names": {"left": "inlet", "right": "outlet",
                                       "bottom": "wall", "top": "wall"}},
           "physics": {"reference_velocity": 1.0, "reference_length": 1.0,
                       "viscosity": {"model": "power_law", "K": 0.05, "n": 0.6}},
           "boundary": {"inlet": {"type": "inlet", "flow_rate": 1.0},
                        "outlet": {"type": "outlet"}, "wall": {"type": "wall"}},
           "solver": {"max_iter": 3000, "tol": 1e-7},
           "output": {"vtk": False, "plots": True}}
    s = run_case(cfg, out_dir=tmp_path, verbose=False, plot=True)
    assert s["converged"]
    assert s["nu"] == pytest.approx(0.05)                   # ν(γ̇_ref = 1)
    v = s["viscosity"]
    assert v["nu_min"] < 0.05 < v["nu_max"]                 # rhéofluidifiant
    assert (tmp_path / "viscosity.png").exists()
