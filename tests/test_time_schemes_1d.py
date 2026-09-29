"""Schémas en temps du solveur 1D : ordres (Womersley exact), stabilité explicite, URANS."""
import warnings

import numpy as np
import pytest

from microrans.cases import PulsatingForcing
from microrans.grid import channel_grid, tanh_grid
from microrans.models import get_model
from microrans.reference import womersley_channel
from microrans.solver import (Solution, explicit_dt_limit, model_diffusivities, solve_steady,
                              solve_unsteady)

NU, AMP, OM = 0.02, 5.0, 2 * np.pi


def _womersley(scheme, steps, n=32):
    g = tanh_grid(n, 2.0)
    m = get_model("laminar", g, NU)
    U0 = womersley_channel(g.y, 0.0, 1.0, AMP, OM, NU)
    r = solve_unsteady(m, g, NU, lambda t: 1.0 + AMP * np.sin(OM * t), Solution(m, g, NU, U0, {}),
                       t_end=1.0, dt=1.0 / steps, scheme=scheme)
    return r.final.U


@pytest.mark.parametrize("scheme,order,steps", [
    ("euler", 1, 40), ("bdf2", 2, 40), ("cn", 2, 40), ("sdirk2", 2, 40), ("sdirk3", 3, 40),
    ("rk2", 2, 400), ("rk3", 3, 400), ("rk4", 4, 400), ("ab2", 2, 400)])
def test_1d_temporal_orders(scheme, order, steps):
    ref = _womersley(scheme, 16 * steps)
    e1 = np.max(np.abs(_womersley(scheme, steps) - ref))
    e2 = np.max(np.abs(_womersley(scheme, 2 * steps) - ref))
    assert np.log2(e1 / e2) == pytest.approx(order, abs=0.35)


@pytest.fixture(scope="module")
def sa_channel():
    re_tau = 395.0
    g = channel_grid(128, re_tau, 0.3)
    m = get_model("sa", g, 1 / re_tau)
    return g, m, solve_steady(m, g, 1 / re_tau, forcing=1.0)


def test_explicit_limit_is_a_real_stability_bound(sa_channel):
    g, m, st = sa_channel
    nu = 1 / 395.0
    lim = explicit_dt_limit(g, *model_diffusivities(m, g, nu, st.U, st.state), scheme="rk3")
    f = PulsatingForcing(1.0, 10.0, 3.95)
    solve_unsteady(m, g, nu, f, st, t_end=300 * lim, dt=0.95 * lim, scheme="rk3")   # stable
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(FloatingPointError):
            solve_unsteady(m, g, nu, f, st, t_end=300 * 4 * lim, dt=4 * lim, scheme="rk3")


def test_turbulent_urans_high_order_dirk_beats_bdf2(sa_channel):
    """À coût comparable, SDIRK3 (16 pas/période) est bien plus précis que BDF2 (64 pas)."""
    g, m, st = sa_channel
    nu, T = 1 / 395.0, 2 * np.pi / 3.95
    f = PulsatingForcing(1.0, 10.0, 3.95)

    def run(scheme, steps):
        return solve_unsteady(m, g, nu, f, st, t_end=T, dt=T / steps, scheme=scheme).final.U

    ref = run("sdirk3", 512)
    e_bdf2 = np.max(np.abs(run("bdf2", 64) - ref))
    e_sd3 = np.max(np.abs(run("sdirk3", 16) - ref))
    e_cn = np.max(np.abs(run("cn", 32) - ref))
    assert e_sd3 < 0.2 * e_bdf2
    assert e_cn < 0.5 * e_bdf2
