"""Thermique : conduction exacte, flux imposé, convection naturelle (de Vahl Davis 1983)."""
import numpy as np
import pytest

from microrans.fv2d import Settings, Solver2D
from microrans.mesh2d import rectangle_mesh


def _box(n, gr=1.0, names=None):
    g = [(0.5, 0.5, gr), (0.5, 0.5, 1 / gr)] if gr != 1.0 else 1.0
    return rectangle_mesh(0, 1, 0, 1, n, n, grading=(g, g), names=names or {
        "left": "hot", "right": "cold", "bottom": "adiab", "top": "adiab"})


def test_conduction_linear_profile_and_heat_balance():
    m = _box(16)
    s = Solver2D(m, 0.1, {"hot": {"type": "wall", "T": 1.0}, "cold": {"type": "wall", "T": 0.0},
                          "adiab": {"type": "wall"}}, energy={"Pr": 1.0},
                 settings=Settings(relax_T=1.0))
    s.run_steady(max_iter=50, tol=1e-12)
    x = m.cell_centers[:, 0]
    assert np.max(np.abs(s.T - (1.0 - x))) < 1e-8
    _, qh = s.wall_heat_flux("hot")
    _, qc = s.wall_heat_flux("cold")
    # flux conductif exact α ΔT / L = 0.1 ; conservation : ce qui entre ressort
    assert np.allclose(qh, 0.1, rtol=1e-6)
    assert np.sum(qh) == pytest.approx(-np.sum(qc), rel=1e-8)


def test_imposed_heat_flux():
    m = _box(16)
    q = 2.0
    s = Solver2D(m, 0.1, {"hot": {"type": "wall", "q": q}, "cold": {"type": "wall", "T": 0.0},
                          "adiab": {"type": "wall"}}, energy={"Pr": 0.5},
                 settings=Settings(relax_T=1.0))
    s.run_steady(max_iter=50, tol=1e-12)
    alpha = 0.1 / 0.5
    Tw, qh = s.wall_heat_flux("hot")
    assert np.allclose(qh, q, rtol=1e-8)
    assert np.allclose(Tw, q / alpha, rtol=1e-6)            # T(0) = q L / α


@pytest.mark.parametrize("Ra,nu_ref,umax_ref,vmax_ref", [(1e4, 2.243, 16.178, 19.617)])
def test_natural_convection_de_vahl_davis(Ra, nu_ref, umax_ref, vmax_ref):
    Pr = 0.71
    m = _box(32, 3.0)
    s = Solver2D(m, Pr, {"hot": {"type": "wall", "T": 0.5}, "cold": {"type": "wall", "T": -0.5},
                         "adiab": {"type": "wall"}},
                 energy={"Pr": Pr, "beta": Ra * Pr, "gravity": (0.0, -1.0)},
                 reference_velocity=20.0)
    assert s.run_steady(max_iter=2000, tol=1e-6)
    _, q = s.wall_heat_flux("hot")
    Nu = float(np.sum(q * s.fvm.magSb[s.patch_slices["hot"]]))
    assert Nu == pytest.approx(nu_ref, rel=0.015)
    # fluide immobile ailleurs que dans la couche : vitesses max proches de la référence
    assert np.abs(s.U[:, 0]).max() == pytest.approx(umax_ref, rel=0.05)
    assert np.abs(s.U[:, 1]).max() == pytest.approx(vmax_ref, rel=0.05)


def test_stable_stratification_stays_at_rest():
    """Chaud en haut, froid en bas : équilibre hydrostatique exact (pas de courants
    parasites grâce à la force évaluée aux faces)."""
    m = _box(20, names={"left": "side", "right": "side", "bottom": "cold", "top": "hot"})
    s = Solver2D(m, 0.71, {"hot": {"type": "wall", "T": 0.5}, "cold": {"type": "wall", "T": -0.5},
                           "side": {"type": "wall"}},
                 energy={"Pr": 0.71, "beta": 1e5 * 0.71, "gravity": (0.0, -1.0)},
                 reference_velocity=1.0)
    s.run_steady(max_iter=300, tol=1e-9)
    assert np.abs(s.U).max() < 1e-6
