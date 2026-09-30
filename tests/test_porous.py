"""Zones poreuses (Darcy-Forchheimer) : Brinkman, perte de charge, anisotropie tournée."""
import numpy as np
import pytest

from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.case import run_case
from microrans.mesh2d import channel_mesh, rectangle_mesh

SLIP = {"inlet": {"type": "inlet", "U": [0.5, 0.0]}, "outlet": {"type": "outlet"},
        "sym": {"type": "symmetry"}}


def _duct(nx=240, length=12.0):
    return rectangle_mesh(0, length, 0, 1, nx, 10, names={"left": "inlet", "right": "outlet",
                                                           "bottom": "sym", "top": "sym"})


def test_brinkman_channel_order():
    """Canal entièrement poreux entre deux parois (équation de Brinkman) : solution exacte
    u = G/(ν d) (1 − cosh(√d (y − h)) / cosh(√d h)) ; ordre → 2."""
    nu, G, d = 0.1, 1.0, 25.0
    errs = []
    for ny in (16, 32, 64):
        m = channel_mesh(1.0, 2.0, 2, ny)
        s = Solver2D(m, nu, {"bottom": {"type": "wall"}, "top": {"type": "wall"}},
                     body_force=(G, 0.0), reference_velocity=G / (nu * d),
                     porous=[{"region": "expression", "expression": "x > -1", "darcy": d}],
                     settings=Settings(relax_U=1.0))
        assert s.run_steady(max_iter=500, tol=1e-11)
        y = m.cell_centers[:, 1]
        ex = G / (nu * d) * (1 - np.cosh(np.sqrt(d) * (y - 1)) / np.cosh(np.sqrt(d)))
        errs.append(np.max(np.abs(s.U[:, 0] - ex)) / ex.max())
    order = np.log2(np.array(errs[:-1]) / errs[1:])
    assert errs[-1] < 3e-3 and order[-1] > 1.7


@pytest.mark.parametrize("angle", [0.0, 30.0])
def test_anisotropic_rotated_plug_pressure_gradients(angle):
    """Bouchon poreux anisotrope (d = 100 / 400, f = 2 / 8) tourné de `angle` dans une
    conduite à parois glissantes : au cœur de la zone u = (U, 0) et ∇p = −K·U exactement
    (K = ν D + ½|U| F), y compris la composante transversale."""
    nu, U, dd, ff = 0.01, 0.5, [100.0, 400.0], [2.0, 8.0]
    m = _duct()
    s = Solver2D(m, nu, SLIP, initial_U=(U, 0.0),
                 porous=[{"name": "plug", "region": "rectangle", "x0": 2, "x1": 10,
                          "y0": -1, "y1": 2, "darcy": dd, "forchheimer": ff,
                          "angle": angle}])
    s.run_steady(max_iter=3000, tol=1e-8)
    a = np.radians(angle)
    c2, s2, cs = np.cos(a) ** 2, np.sin(a) ** 2, np.sin(a) * np.cos(a)
    kxx = nu * (dd[0] * c2 + dd[1] * s2) + 0.5 * U * (ff[0] * c2 + ff[1] * s2)
    kyx = nu * (dd[0] - dd[1]) * cs + 0.5 * U * (ff[0] - ff[1]) * cs
    C = m.cell_centers
    col = np.abs(C[:, 0] - 6.025) < 1e-6
    row = np.abs(C[:, 1] - 0.45) < 1e-6
    sel = row & (C[:, 0] > 4) & (C[:, 0] < 8)
    assert np.polyfit(C[sel, 0], s.p[sel], 1)[0] == pytest.approx(-kxx * U, rel=1e-5)
    assert np.polyfit(C[col, 1], s.p[col], 1)[0] == pytest.approx(-kyx * U, rel=1e-5,
                                                                  abs=1e-10)
    assert np.allclose(s.U[col, 0], U, rtol=1e-5) and np.abs(s.U[col, 1]).max() < 1e-5


def test_porous_transient_and_case_report(tmp_path):
    """Instationnaire (schéma auto → BDF2 implicite) et bilan par zone : puissance dissipée
    ≈ perte de charge × débit dans le bouchon (Darcy seul, zone de 2 × 1). Écart 0.5 % :
    oscillation de la vitesse reconstruite (±2.6 %) dans les cellules voisines des
    interfaces fluide / poreux (Rhie-Chow au saut de résistance, limite documentée)."""
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 6, "y0": 0, "y1": 1, "nx": 60,
                    "ny": 6, "names": {"left": "inlet", "right": "outlet", "bottom": "sym",
                                       "top": "sym"}},
           "physics": {"nu": 0.01},
           "initial": {"U": [0.5, 0.0]},
           "porous": [{"name": "filtre", "region": "rectangle", "x0": 2, "x1": 4, "y0": -1,
                       "y1": 2, "permeability": 0.01}],
           "boundary": SLIP,
           "solver": {"mode": "transient", "dt": 0.05, "t_end": 2.0},
           "output": {"vtk": False, "plots": False, "checkpoint": False}}
    s, solver = run_case(cfg, out_dir=tmp_path, verbose=False, plot=False,
                         return_solver=True)
    assert solver.settings.time_scheme == "backward"
    z = s["porous"][0]
    assert z["name"] == "filtre" and z["cells"] == 120
    assert z["dissipation"] == pytest.approx(0.01 * 100 * 0.25 * 2.0, rel=1e-2)
    assert float(np.sum(solver.F_b[solver.patch_slices["outlet"]])) == pytest.approx(0.5)


def test_porous_validation_errors():
    m = _duct(24, 6.0)
    with pytest.raises(ValueError, match="aucune cellule"):
        Solver2D(m, 0.01, SLIP, porous=[{"x0": 20, "x1": 30, "y0": 0, "y1": 1,
                                         "darcy": 1.0}])
    with pytest.raises(ValueError, match="darcy"):
        Solver2D(m, 0.01, SLIP, porous=[{"x0": 0, "x1": 3, "y0": 0, "y1": 1}])
    with pytest.raises(ValueError, match="recouvre"):
        Solver2D(m, 0.01, SLIP, porous=[{"x0": 0, "x1": 3, "y0": 0, "y1": 1, "darcy": 1.0},
                                        {"region": "circle", "center": [2, 0.5],
                                         "radius": 0.5, "darcy": 1.0}])
