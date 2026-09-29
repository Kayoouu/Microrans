"""Schémas en temps 2D : ordre de convergence, précision, pas adaptatif, choix automatique."""
import warnings

import numpy as np
import pytest

from microrans.fv2d import Settings, Solver2D
from microrans.mesh2d import channel_mesh, rectangle_mesh
from microrans.reference import womersley_channel

NU, AMP, OM = 0.02, 2.0, 2 * np.pi
WALLS = {"bottom": {"type": "wall"}, "top": {"type": "wall"}}


def _womersley(scheme, steps, ny=16):
    """Canal pulsé (Womersley) : pression uniforme → erreur de pure intégration en temps."""
    m = channel_mesh(4.0, 2.0, 2, ny)          # Δx grand : Courant axial faible
    s = Solver2D(m, NU, WALLS, settings=Settings(time_scheme=scheme), reference_velocity=1.0,
                 body_force=lambda t: (1.0 + AMP * np.sin(OM * t), 0.0))
    s.U[:, 0] = womersley_channel(m.cell_centers[:, 1], 0.0, 1.0, AMP, OM, NU)
    s.F_i = np.sum(s.fvm.interp(s.U) * s.fvm.Si, axis=1)
    s.run_transient(1.0 / steps, 1.0)
    return s.U[:, 0].copy()


@pytest.mark.parametrize("scheme,order", [("euler", 1), ("backward", 2), ("crankNicolson", 2),
                                          ("rk2", 2), ("rk3", 3), ("rk4", 4)])
def test_temporal_order_womersley(scheme, order):
    ref = _womersley(scheme, 640)
    e1 = np.max(np.abs(_womersley(scheme, 40) - ref))
    e2 = np.max(np.abs(_womersley(scheme, 80) - ref))
    assert np.log2(e1 / e2) == pytest.approx(order, abs=0.3)


L2PI = 2 * np.pi
U0 = np.array([1.0, 0.5])


def _tgv_exact(C, t, nu):
    F = np.exp(-2 * nu * t)
    x, y = C[:, 0] - U0[0] * t, C[:, 1] - U0[1] * t
    return np.column_stack([U0[0] - np.cos(x) * np.sin(y) * F, U0[1] + np.sin(x) * np.cos(y) * F])


def _tgv_solver(n, nu=0.01, **kw):
    m = rectangle_mesh(0, L2PI, 0, L2PI, n, n, names={"left": "L", "right": "R",
                                                      "bottom": "B", "top": "T"},
                       periodic=[("L", "R"), ("B", "T")])
    s = Solver2D(m, nu, {}, settings=Settings(**kw), reference_velocity=1.0)
    s.U = _tgv_exact(m.cell_centers, 0.0, nu)
    s.F_i = np.sum(s.fvm.interp(s.U) * s.fvm.Si, axis=1)
    return s


@pytest.mark.parametrize("scheme,tol", [("backward", 0.02), ("crankNicolson", 0.02),
                                        ("rk3", 0.011), ("rk4", 0.011)])
def test_convected_taylor_green_vortex(scheme, tol):
    """Tourbillon advecté (solution exacte par invariance galiléenne), Courant ≈ 1."""
    s = _tgv_solver(32, time_scheme=scheme)
    s.run_transient(0.08, 2.0)
    err = np.max(np.abs(s.U - _tgv_exact(s.mesh.cell_centers, 2.0, 0.01)))
    assert err < tol
    # flux aux faces à divergence nulle (projection / PISO)
    assert np.max(np.abs(s.fvm.div(s.F_i, s.F_b))) < 1e-8


def test_adaptive_time_step_respects_courant():
    s = _tgv_solver(24, time_scheme="rk3", adjust_dt=True, max_co=0.8)
    hist = s.run_transient(1.0, 1.0)                    # Δt initial démesuré : corrigé
    assert s.time == pytest.approx(1.0)
    assert max(h["Co"] for h in hist) < 0.8 * 1.05
    assert all(np.isfinite(h["Co"]) for h in hist)


def test_explicit_scheme_warns_beyond_stability_limit():
    s = _tgv_solver(16, time_scheme="rk4")
    with pytest.warns(UserWarning, match="stabilité"):
        with warnings.catch_warnings():
            warnings.simplefilter("default")
            try:
                s.run_transient(0.5, 0.5)
            except FloatingPointError:
                pass


def test_auto_scheme_choice():
    # tourbillon laminaire : limite convective → RK3 explicite
    s = _tgv_solver(16, time_scheme="auto")
    s.run_transient(0.05, 0.1)
    assert s.settings.time_scheme == "rk3"
    # paroi très raffinée (y⁺ ~ 0.1) : limite de diffusion → BDF2 implicite
    m = channel_mesh(1.0, 2.0, 2, 64, first_height=1e-4)
    s = Solver2D(m, 1e-3, WALLS, body_force=(1.0, 0.0), initial_U=(1.0, 0.0),
                 settings=Settings(time_scheme="auto"))
    s.run_transient(0.05, 0.1)
    assert s.settings.time_scheme == "backward"
