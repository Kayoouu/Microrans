"""Noyaux Numba facultatifs ([solver] numba = true) : mêmes résultats que NumPy."""
import numpy as np
import pytest

pytest.importorskip("numba")

from microrans.fv2d import Settings, Solver2D                        # noqa: E402
from microrans.fv2d.fvm import FVM                                    # noqa: E402
from microrans.mesh2d import Circle, Rectangle, channel_mesh, rectangle_mesh, triangulate  # noqa: E402


@pytest.mark.parametrize("axi", [False, True])
def test_operators_match_numpy(axi):
    m = triangulate(Rectangle(0, 0.2, 2, 1.2) - Circle((1.0, 0.7), 0.2), 0.15, max_iter=60)
    a, b = FVM(m, axisymmetric=axi), FVM(m, axisymmetric=axi, numba=True)
    assert b.fast and not a.fast
    rng = np.random.default_rng(0)
    phi, pb = rng.normal(size=a.nc), rng.normal(size=a.nb)
    fi, fb = rng.normal(size=a.ni), rng.normal(size=a.nb)
    for arr in (phi, rng.normal(size=(a.nc, 2)), rng.normal(size=(a.nc, 2, 2))):
        assert np.allclose(a.interp(arr), b.interp(arr), rtol=0, atol=1e-13)
    assert np.allclose(a.sum_faces(fi, fb), b.sum_faces(fi, fb), rtol=0, atol=1e-12)
    assert np.allclose(a.grad(phi, pb), b.grad(phi, pb), rtol=1e-12, atol=1e-12)
    for idx in ("P", "N", "Pb"):
        w = rng.normal(size=len(getattr(a, idx)))
        assert np.allclose(a._sum(getattr(a, idx), w), b._sum(getattr(b, idx), w), atol=1e-12)
    d, up, lo = rng.normal(size=a.nc), rng.normal(size=a.ni), rng.normal(size=a.ni)
    assert abs(a.matrix(d, up, lo) - b.matrix(d, up, lo)).max() < 1e-13


def test_solver_runs_match_numpy_and_thread_count():
    def run(**kw):
        m = channel_mesh(1.0, 2.0, 2, 24, first_height=2e-3)
        s = Solver2D(m, 1 / 2000, {"bottom": {"type": "wall"}, "top": {"type": "wall"}},
                     model="sst", body_force=(1.0, 0.0), initial_U=(20.0, 0.0),
                     settings=Settings(**kw))
        s.run_steady(max_iter=15, tol=1e-30)
        return s
    ref, one, two = run(), run(numba=True), run(numba=True, threads=2)
    assert one.fvm.fast and not ref.fvm.fast
    assert np.max(np.abs(ref.U - one.U)) < 1e-8
    assert np.array_equal(one.U, two.U)          # sommes par cellule : ordre indépendant des fils


def test_cavity_same_converged_solution():
    def run(**kw):
        m = rectangle_mesh(0, 1, 0, 1, 24, 24, names={"left": "w", "right": "w",
                                                       "bottom": "w", "top": "lid"})
        s = Solver2D(m, 0.01, {"lid": {"type": "wall", "U": [1.0, 0.0]}, "w": {"type": "wall"}},
                     settings=Settings(**kw))
        assert s.run_steady(max_iter=2000, tol=1e-8)
        return s
    a, b = run(), run(numba=True)
    assert np.max(np.abs(a.U - b.U)) < 1e-6
