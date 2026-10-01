"""Opérateurs volumes finis en 3D (fvm.py, dimension lue sur le maillage) : gradient,
limiteur, laplacien ; mêmes fonctions qu'en 2D."""
import numpy as np
import pytest
import scipy.sparse.linalg as spla

from microrans.fv2d.fvm import FVM
from microrans.mesh3d import box_mesh


def _lin(X):
    return 1.0 + 2.0 * X[:, 0] - 3.0 * X[:, 1] + 0.5 * X[:, 2]


@pytest.mark.parametrize("numba", [False, True])
def test_gradient_of_linear_field_exact(numba):
    pytest.importorskip("numba") if numba else None
    m = box_mesh(0, 1, 0, 2, 0, 1, 5, 6, 7, grading=(2.0, 1.0, 0.5))
    f = FVM(m, numba=numba)
    assert f.dim == 3
    g = f.grad(_lin(m.cell_centers), _lin(m.face_centers[f.ni:]))
    assert g.shape == (m.n_cells, 3)
    assert np.abs(g - [2.0, -3.0, 0.5]).max() < 1e-12


def test_limiter_keeps_face_values_within_neighbours():
    m = box_mesh(0, 1, 0, 1, 0, 1, 6, 6, 6)
    f = FVM(m)
    phi = np.where(m.cell_centers[:, 0] < 0.5, 1.0, 0.0)          # saut
    phib = phi[m.owner[f.ni:]]
    g = f.limit_grad(phi, f.grad(phi, phib), phib)
    ni = f.ni
    P, N = m.owner[:ni], m.neighbour
    for cells, other in ((P, N), (N, P)):
        face = phi[cells] + np.sum(g[cells] * (m.face_centers[:ni] - m.cell_centers[cells]), 1)
        lo = np.minimum(phi[cells], phi[other]) - 1e-12
        assert np.all(face >= np.minimum(lo, 0.0)) and np.all(face <= 1.0 + 1e-12)


def _exact(X):
    x, y, z = X.T
    return np.exp(x) * np.sin(2 * y) * (1 + z * z)


def _laplacian_exact(X):
    x, y, z = X.T
    return np.exp(x) * np.sin(2 * y) * (-3.0 * (1 + z * z) + 2.0)


def test_laplacian_second_order_on_graded_hexahedra():
    errs = []
    for n in (8, 16):
        m = box_mesh(0, 1, 0, 1, 0, 1, n, n, n, grading=(2.0, 0.5, 3.0))
        f = FVM(m)
        phib = _exact(m.face_centers[f.ni:])
        bc = (np.zeros(f.nb), phib, -1.0 / f.dperp, phib / f.dperp)
        d, u, lo, r = f.assemble(np.zeros(f.ni), np.zeros(f.nb), np.ones(f.ni), np.ones(f.nb),
                                 bc)
        phi = spla.spsolve(f.matrix(d, u, lo).tocsc(), r - _laplacian_exact(m.cell_centers) * f.V)
        errs.append(np.sqrt(np.sum((phi - _exact(m.cell_centers)) ** 2 * f.V) / f.V.sum()))
    assert np.log2(errs[0] / errs[1]) > 1.8                       # mesuré : 1.86


def test_axisymmetric_refused_in_3d():
    with pytest.raises(ValueError, match="maillage 2D seulement"):
        FVM(box_mesh(0, 1, 0, 1, 0, 1, 2, 2, 2), axisymmetric=True)
