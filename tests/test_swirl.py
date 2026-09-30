"""Rotation propre (swirl, axisymétrique) et disques actuateurs."""
import numpy as np
import pytest

from microrans.fv2d import Solver2D
from microrans.mesh2d import rectangle_mesh


def _couette(ny, R1=1.0, R2=2.0):
    return rectangle_mesh(0, 0.2, R1, R2, 2, ny, names={"left": "i", "right": "o",
                                                        "bottom": "inner", "top": "outer"},
                          periodic=[("i", "o")])


def test_taylor_couette_velocity_pressure_torque():
    """Cylindres coaxiaux infinis, intérieur tournant : u_θ = A r + B/r, ∂p/∂r = u_θ²/r,
    couple 4πνΩR1²R2²/(R2² − R1²) par unité de longueur ; ordre 2, couples équilibrés."""
    R1, R2, Om, nu = 1.0, 2.0, 1.0, 0.5
    A = -Om * R1 ** 2 / (R2 ** 2 - R1 ** 2)
    B = Om * R1 ** 2 * R2 ** 2 / (R2 ** 2 - R1 ** 2)
    Tex = -4 * np.pi * nu * Om * R1 ** 2 * R2 ** 2 / (R2 ** 2 - R1 ** 2) * 0.2
    errs, terr = [], []
    for ny in (8, 16, 32):
        m = _couette(ny)
        s = Solver2D(m, nu, {"inner": {"type": "wall", "omega": Om},
                             "outer": {"type": "wall"}}, axisymmetric=True, swirl=True,
                     reference_velocity=Om * R1)
        s.run_steady(max_iter=1500, tol=1e-9)
        r = m.cell_centers[:, 1]
        errs.append(np.max(np.abs(s.scalars["U_theta"] - (A * r + B / r))) / (Om * R1))
        terr.append(abs(s.torque("inner") / Tex - 1))
        assert abs(s.torque("inner") + s.torque("outer")) < 1e-3 * abs(Tex)
    assert np.all(np.log2(np.array(errs[:-1]) / errs[1:]) > 1.8)
    assert np.all(np.log2(np.array(terr[:-1]) / terr[1:]) > 1.8) and terr[-1] < 3e-4
    # équilibre radial de la pression (dernier maillage)
    rr = np.linspace(R1, R2, 4001)
    w2r = (A * rr + B / rr) ** 2 / rr
    pex = np.concatenate([[0], np.cumsum(0.5 * (w2r[1:] + w2r[:-1]) * np.diff(rr))])
    pe = np.interp(r, rr, pex)
    assert np.max(np.abs((s.p - s.p.mean()) - (pe - pe.mean()))) < 5e-3 * np.ptp(pe)


def test_rotating_pipe_solid_body_rotation():
    """Tuyau tournant (Ω = 2) : rotation solide u_θ = Ω r, p(R) − p(0) = Ω²R²/2."""
    m = rectangle_mesh(0, 0.2, 0, 1.0, 2, 20, names={"left": "i", "right": "o",
                                                      "bottom": "axis", "top": "wall"},
                       periodic=[("i", "o")])
    s = Solver2D(m, 0.2, {"axis": {"type": "axis"}, "wall": {"type": "wall", "omega": 2.0}},
                 axisymmetric=True, swirl=True, reference_velocity=2.0)
    s.run_steady(max_iter=1500, tol=1e-9)
    r = m.cell_centers[:, 1]
    assert np.max(np.abs(s.scalars["U_theta"] - 2.0 * r)) < 1e-4       # itératif : 2e-5
    col = np.argsort(r)
    dp = np.polyfit(r[col] ** 2, s.p[col], 1)[0]              # p = Ω² r² / 2 + c
    assert dp == pytest.approx(2.0, rel=1e-2)


def _disk_case(nx, ny, **disk):
    gx = [(0.36, 0.25, 1 / 8.0), (0.08, 0.35, 1.0), (0.56, 0.40, 10.0)]
    gy = [(0.15, 0.5, 1.0), (0.85, 0.5, 15.0)]
    m = rectangle_mesh(-20, 30, 0, 10, nx, ny, grading=(gx, gy),
                       names={"left": "inlet", "right": "outlet", "bottom": "axis",
                              "top": "top"})
    return m, Solver2D(m, 0.01, {"inlet": {"type": "inlet", "U": [1.0, 0.0]},
                                 "outlet": {"type": "outlet"}, "axis": {"type": "axis"},
                                 "top": {"type": "symmetry"}},
                       axisymmetric=True, initial_U=(1.0, 0.0), reference_velocity=1.0,
                       swirl="torque" in disk,
                       actuator_disks=[{"name": "rotor", "x0": -0.075, "x1": 0.075,
                                        "radius": 1.0, **disk}])


def test_actuator_disk_froude_momentum_theory():
    """Éolienne C_T = 0.5 (disque uniforme, Re = U R/ν = 100) : vitesse au disque et
    puissance à < 1.5 % de la théorie de Froude (a = 0.146) ; mesuré −0.8 %."""
    _, s = _disk_case(120, 60, thrust_coefficient=0.5, mode="turbine")
    assert s.run_steady(max_iter=2000, tol=1e-7)
    a = (1 - np.sqrt(0.5)) / 2
    rep = s.disk_report()[0]
    assert rep["thrust"] == pytest.approx(-0.5 * 0.5 * np.pi)
    assert rep["U_disk"] == pytest.approx(1 - a, rel=0.015)
    assert rep["power"] == pytest.approx(4 * a * (1 - a) ** 2 * 0.5 * np.pi, rel=0.015)


def test_actuator_disk_torque_gives_angular_momentum_flux():
    """Couple Q imposé : le flux de moment cinétique r u_θ qui sort du domaine vaut Q
    (conservation ; parois glissantes, pas de couple pariétal)."""
    m, s = _disk_case(120, 60, thrust=-0.3, torque=0.2)
    s.run_steady(max_iter=2000, tol=1e-7)
    sl = s.patch_slices["outlet"]
    rb = m.face_centers[m.n_internal:][sl][:, 1]
    Wb = s.boundary_scalar("U_theta")[sl]
    flux = 2 * np.pi * float(np.sum(s.F_b[sl] * rb * Wb))
    assert flux == pytest.approx(0.2, rel=2e-2)


def test_swirl_requires_axisymmetric():
    m = _couette(4)
    with pytest.raises(ValueError, match="axisymétrique"):
        Solver2D(m, 0.1, {"inner": {"type": "wall"}, "outer": {"type": "wall"}}, swirl=True)
