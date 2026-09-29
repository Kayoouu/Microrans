"""Solveurs linéaires : AMG par agrégation, Krylov, assemblage CSR à structure figée."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from microrans.fv2d.fvm import FVM
from microrans.linalg import AggregationAMG, JacobiPrecond, LinearSolver, fcg, pbicgstab
from microrans.mesh2d import Circle, channel_mesh, o_grid


def _laplacian(mesh):
    f = FVM(mesh)
    nb = f.nb
    a, b, g, d = np.zeros(nb), np.zeros(nb), -1.0 / f.dperp, np.zeros(nb)   # Dirichlet 0
    diag, up, lo, _ = f.assemble(np.zeros(f.ni), np.zeros(nb), np.ones(f.ni), np.ones(nb),
                                 (a, b, g, d))
    return f, f.matrix(diag, up, lo)


def test_csr_assembly_matches_coo_with_duplicate_faces():
    # canal périodique de 2 cellules en x : deux faces relient la même paire de cellules
    m = channel_mesh(1.0, 2.0, 2, 8)
    f = FVM(m)
    rng = np.random.default_rng(1)
    diag, up, lo = rng.random(f.nc), rng.random(f.ni), rng.random(f.ni)
    ref = sp.coo_matrix((np.concatenate([diag, up, lo]), (f.rows, f.cols)),
                        shape=(f.nc, f.nc)).toarray()
    assert np.allclose(f.matrix(diag, up, lo).toarray(), ref)


def test_amg_coarsens_and_converges_uniformly():
    its = []
    for n in (40, 80):
        m = o_grid(Circle((0, 0), 0.5, "c"), 2 * n, n, 20.0, 0.5 / n)
        f, A = _laplacian(m)
        b = np.random.default_rng(0).standard_normal(f.nc) * f.V
        H = AggregationAMG(f.nc, f.P, f.N, f.g)
        ratios = [nc / n0 for (_, nc), n0 in zip(H.levels, [f.nc] + [nc for _, nc in H.levels])]
        assert max(ratios) < 0.35                       # agrégats de ~4 cellules
        x, k, ok = fcg(A, b, np.zeros_like(b), H.setup(A), rtol=1e-8, maxiter=200)
        assert ok
        ref = spsolve(A.tocsc(), b)
        assert np.linalg.norm(x - ref) < 1e-6 * np.linalg.norm(ref)
        its.append(k)
    assert its[1] < 2.0 * its[0] + 5                    # quasi indépendant de la taille


def test_chebyshev_smoother_and_bicgstab():
    m = o_grid(Circle((0, 0), 0.5, "c"), 64, 32, 20.0, 0.01)
    f, A = _laplacian(m)
    b = np.ones(f.nc) * f.V
    H = AggregationAMG(f.nc, f.P, f.N, f.g, smoother="chebyshev")
    x, k, ok = fcg(A, b, np.zeros_like(b), H.setup(A), rtol=1e-8)
    assert ok and k < 60
    # matrice non symétrique (convection upwind) : BiCGStab + Jacobi
    F = np.random.default_rng(2).standard_normal(f.ni) * 0.1
    nb = f.nb
    diag, up, lo, _ = f.assemble(F, np.zeros(nb), np.ones(f.ni), np.ones(nb),
                                 (np.zeros(nb), np.zeros(nb), -1.0 / f.dperp, np.zeros(nb)))
    An = f.matrix(diag, up, lo)
    x, k, ok = pbicgstab(An, b, np.zeros_like(b), JacobiPrecond(An), rtol=1e-8, maxiter=2000)
    assert ok
    assert np.linalg.norm(An @ x - b) <= 1e-7 * np.linalg.norm(b)


def test_linear_solver_facade_methods_agree():
    m = o_grid(Circle((0, 0), 0.5, "c"), 48, 24, 20.0, 0.01)
    f, A = _laplacian(m)
    b = np.random.default_rng(3).standard_normal(f.nc)
    ls = LinearSolver(f.nc, f.P, f.N, f.g)
    ref = spsolve(A.tocsc(), b)
    for method in ("direct", "amg", "bicgstab", "auto"):
        x = ls.solve(A, b, None, method, rtol=1e-10, maxiter=5000, symmetric=True)
        assert np.linalg.norm(x - ref) < 1e-6 * np.linalg.norm(ref), method
