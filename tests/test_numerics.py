"""Vérification du noyau numérique contre des solutions exactes."""
import numpy as np
import pytest

from microrans import channel_grid, tanh_grid
from microrans.numerics import d2dy2, ddy, integrate
from microrans.verification import (diffusion_mms, poiseuille_error, womersley_space_order,
                                    womersley_time_order)


def test_channel_grid_hits_target_y1plus_and_is_symmetric():
    g = channel_grid(128, re_tau=395.0, y1_plus=0.5)
    assert g.dy[0] * 395.0 == pytest.approx(0.5, rel=1e-8)
    assert np.allclose(g.y + g.y[::-1], 2.0, atol=1e-12)
    assert np.all(np.diff(g.y) > 0)


def test_channel_grid_falls_back_to_uniform():
    g = channel_grid(64, re_tau=10.0, y1_plus=5.0)
    assert np.allclose(np.diff(g.y), 2.0 / 64)


def test_derivatives_exact_for_quadratics_on_stretched_grid():
    g = tanh_grid(40, 2.5)
    f = 3.0 * g.y ** 2 - 2.0 * g.y + 1.0
    assert np.allclose(ddy(g, f), 6.0 * g.y - 2.0, atol=1e-9)
    assert np.allclose(d2dy2(g, f)[1:-1], 6.0, atol=1e-7)


def test_trapezoid_integral():
    g = tanh_grid(400, 2.0)
    assert integrate(g, np.sin(0.5 * np.pi * g.y)) == pytest.approx(4.0 / np.pi, rel=1e-4)


def test_diffusion_operator_second_order():
    study = diffusion_mms()
    assert study.passed, study.report()


def test_laminar_poiseuille_exact():
    assert poiseuille_error() < 1e-10


@pytest.mark.parametrize("scheme,order", [("euler", 1.0), ("bdf2", 2.0)])
def test_womersley_temporal_order(scheme, order):
    study = womersley_time_order(scheme)
    assert study.passed, study.report()
    assert study.errors[-1] < 2e-3


def test_womersley_spatial_order():
    study = womersley_space_order()
    assert study.passed, study.report()
