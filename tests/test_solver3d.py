"""Solveur incompressible en 3D (même Solver2D, 3 composantes) : solution exacte de la
conduite carrée, équivalence exacte avec le 2D sur un cas extrudé, options refusées."""
import numpy as np
import pytest

from microrans.fv2d.solver import Settings, Solver2D
from microrans.mesh2d import Circle, cavity_mesh, o_grid
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
