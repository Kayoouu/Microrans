"""Conditions aux limites : débit imposé, pression totale imposée."""
import numpy as np
import pytest

from microrans.fv2d import Settings, Solver2D
from microrans.mesh2d import rectangle_mesh


def _channel(nx=60, ny=16, length=6.0, height=1.0, names=None):
    return rectangle_mesh(0, length, 0, height, nx, ny, names=names or {
        "left": "inlet", "right": "outlet", "bottom": "wall", "top": "wall"})


@pytest.mark.parametrize("profile", ["uniform", "parabolic"])
def test_flow_rate_inlet_planar(profile):
    """Débit par unité de profondeur conservé exactement ; profil parabolique : entrée déjà
    établie (profil de sortie identique à Poiseuille)."""
    m = _channel()
    s = Solver2D(m, 0.02, {"inlet": {"type": "inlet", "flow_rate": 0.5, "profile": profile},
                           "outlet": {"type": "outlet"}, "wall": {"type": "wall"}})
    assert s.run_steady(max_iter=3000, tol=1e-8)
    assert np.sum(s.F_b[s.patch_slices["outlet"]]) == pytest.approx(0.5, rel=1e-9)
    Ub = s.boundary_U(s.U)[s.patch_slices["inlet"]]
    assert np.all(Ub[:, 0] > 0) and np.allclose(Ub[:, 1], 0.0)
    if profile == "parabolic":
        y = m.face_centers[m.n_internal:][s.patch_slices["inlet"]][:, 1]
        assert np.allclose(Ub[:, 0] / Ub[:, 0].max(), y * (1 - y) / (y * (1 - y)).max())


def test_flow_rate_inlet_axisymmetric_total_over_360_degrees():
    m = _channel(names={"left": "inlet", "right": "outlet", "bottom": "axis", "top": "wall"},
                 height=0.5)
    q = np.pi * 0.25                                    # U_b = 1 dans un tuyau de rayon 0.5
    s = Solver2D(m, 0.02, {"inlet": {"type": "inlet", "flow_rate": q, "profile": "parabolic"},
                           "outlet": {"type": "outlet"}, "wall": {"type": "wall"},
                           "axis": {"type": "axis"}}, axisymmetric=True)
    assert s.run_steady(max_iter=3000, tol=1e-8)
    assert 2 * np.pi * np.sum(s.F_b[s.patch_slices["outlet"]]) == pytest.approx(q, rel=1e-9)
    assert s.U[:, 0].max() == pytest.approx(2.0, rel=0.01)   # Poiseuille : 2 U_b sur l'axe


def test_total_pressure_inlet_poiseuille():
    """Écoulement entraîné par p0 − p_sortie = 0.1 : p + ½|U|² = p0 sur l'entrée ; débit à
    1 % du débit de Poiseuille (perte d'entrée et pression dynamique ~1 %)."""
    nu, dp, L = 0.02, 0.1, 10.0
    m = _channel(100, 20, L)
    s = Solver2D(m, nu, {"inlet": {"type": "pressure_inlet", "p0": dp},
                         "outlet": {"type": "outlet", "p": 0.0}, "wall": {"type": "wall"}},
                 reference_velocity=0.05)
    assert s.run_steady(max_iter=5000, tol=1e-8)
    sl = s.patch_slices["inlet"]
    pb, Ub = s.boundary_p(s.p)[sl], s.boundary_U(s.U)[sl]
    assert np.allclose(pb + 0.5 * np.sum(Ub ** 2, axis=1), dp, rtol=1e-12)
    q = np.sum(s.F_b[s.patch_slices["outlet"]])
    assert q == pytest.approx(dp / (12 * nu * L), rel=0.012)


def test_total_pressure_inlet_explicit_scheme():
    """RK3 explicite : la pression imposée variable est prise en compte à chaque projection
    (le second membre de Poisson est recalculé)."""
    m = _channel(40, 10, 4.0)
    s = Solver2D(m, 0.05, {"inlet": {"type": "pressure_inlet", "p0": 0.05},
                           "outlet": {"type": "outlet", "p": 0.0}, "wall": {"type": "wall"}},
                 reference_velocity=0.1,
                 settings=Settings(time_scheme="rk3", adjust_dt=True, max_co=0.5))
    s.run_transient(0.01, 20.0)
    sl = s.patch_slices["inlet"]
    pb, Ub = s.boundary_p(s.p)[sl], s.boundary_U(s.U)[sl]
    assert np.allclose(pb + 0.5 * np.sum(Ub ** 2, axis=1), 0.05, rtol=1e-3)
    assert np.sum(s.F_b[sl]) < 0                         # le fluide entre
