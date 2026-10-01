"""Chemin GPU (CuPy) vérifié sans carte graphique grâce au « faux GPU » de fake_device.py :
tout mélange implicite CPU/GPU lève une erreur, et les résultats doivent égaler le CPU."""
import warnings

import numpy as np
import pytest

from fake_device import DeviceArray, MixingError, register_fake_gpu
from microrans.backend import get_backend
from microrans.fv2d import Settings, Solver2D
from microrans.mesh2d import Circle, Rectangle, cavity_mesh, channel_mesh, rectangle_mesh
from microrans.mesh2d.unstructured import triangulate

register_fake_gpu()


def _pair(build, run, **settings):
    out = {}
    for be in ("cpu", "fakegpu"):
        s = build(Settings(backend=be, **settings))
        if be == "fakegpu":
            assert isinstance(s.U, DeviceArray) and isinstance(s.fvm.V, DeviceArray)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run(s)
        s.to_cpu()
        assert isinstance(s.U, np.ndarray)
        out[be] = s
    return out["cpu"], out["fakegpu"]


def test_fake_device_rejects_mixing():
    xp = get_backend("fakegpu").xp
    a = xp.zeros(3)
    with pytest.raises(MixingError):
        a + np.ones(3)
    with pytest.raises(MixingError):
        np.asarray(a)
    with pytest.raises(MixingError):
        a[np.array([0, 1])]


@pytest.mark.parametrize("solver_p", ["auto", "amg"])
def test_steady_laminar_cavity(solver_p):
    def build(st):
        return Solver2D(cavity_mesh(20), 0.01, {"lid": {"type": "wall", "U": [1, 0]},
                                                 "walls": {"type": "wall"}},
                        reference_velocity=1.0, settings=st)
    c, g = _pair(build, lambda s: s.run_steady(max_iter=200, tol=1e-7), solver_p=solver_p)
    assert c.iterations == g.iterations or abs(c.iterations - g.iterations) < 5
    assert np.max(np.abs(c.U - g.U)) < 1e-5


@pytest.mark.parametrize("model", ["sa", "sst", "ke", "kw"])
def test_turbulent_nonorthogonal_open_domain(model):
    """Triangles (correction non orthogonale), entrée / sortie / symétrie / paroi."""
    dom = Rectangle(-2, -1.5, 4, 1.5, names={"left": "inlet", "right": "outlet",
                                             "bottom": "side", "top": "side"})
    m = triangulate(dom - Circle((0, 0), 0.3, "body").as_wall(), 0.35,
                    [{"shape": Circle((0, 0), 0.3), "h": 0.1, "growth": 0.3}], max_iter=60)

    def build(st):
        return Solver2D(m, 1e-3, {"inlet": {"type": "inlet", "U": [1, 0]},
                                  "outlet": {"type": "outlet"}, "side": {"type": "symmetry"},
                                  "body": {"type": "wall"}}, model=model, initial_U=(1, 0),
                        settings=st)
    # solveurs directs : même arithmétique des deux côtés (l'AMG utilise Gauss-Seidel sur
    # CPU et Chebyshev sur GPU, d'où des itérés intermédiaires légèrement différents)
    c, g = _pair(build, lambda s: s.run_steady(max_iter=25, tol=1e-12), solver_p="direct")
    assert np.max(np.abs(c.U - g.U)) < 1e-6
    for k in c.state:
        assert np.max(np.abs(c.state[k] - g.state[k])) <= 1e-6 * np.max(np.abs(c.state[k]))


@pytest.mark.parametrize("scheme", ["backward", "crankNicolson", "rk3", "ab2"])
def test_transient_schemes_periodic(scheme):
    L = 2 * np.pi

    def build(st):
        m = rectangle_mesh(0, L, 0, L, 12, 12, names={"left": "L", "right": "R",
                                                       "bottom": "B", "top": "T"},
                           periodic=[("L", "R"), ("B", "T")])
        s = Solver2D(m, 0.05, {}, settings=st, reference_velocity=1.0)
        C = m.cell_centers
        U0 = np.column_stack([-np.cos(C[:, 0]) * np.sin(C[:, 1]),
                              np.sin(C[:, 0]) * np.cos(C[:, 1])])
        s.U = s.backend.asarray(U0)
        s.F_i = s.xp.sum(s.fvm.interp(s.U) * s.fvm.Si, axis=1)
        return s
    c, g = _pair(build, lambda s: s.run_transient(0.05, 0.5), time_scheme=scheme)
    assert np.max(np.abs(c.U - g.U)) < 1e-7


def test_channel_body_force_and_forces_on_gpu():
    def build(st):
        return Solver2D(channel_mesh(1.0, 2.0, 2, 24), 0.1,
                        {"bottom": {"type": "wall"}, "top": {"type": "wall"}},
                        body_force=(1.0, 0.0), settings=st)
    c, g = _pair(build, lambda s: s.run_steady(max_iter=100, tol=1e-9), monitor_tol=1e-9)
    assert np.max(np.abs(c.U - g.U)) < 1e-8
    fc, fg = c.forces()["bottom"]["total"], g.forces()["bottom"]["total"]
    assert np.allclose(fc, fg)


def test_wall_functions_on_gpu():
    def build(st):
        return Solver2D(channel_mesh(1.0, 2.0, 2, 16, first_height=0.05), 1 / 2000.0,
                        {"bottom": {"type": "wall"}, "top": {"type": "wall"}}, model="sst",
                        body_force=(1.0, 0.0), initial_U=(20.0, 0.0), settings=st)
    c, g = _pair(build, lambda s: s.run_steady(max_iter=30, tol=1e-12),
                 wall_treatment="wall_function", solver_p="direct")
    assert np.max(np.abs(c.U - g.U)) < 1e-6


def test_scalars_and_non_newtonian_on_gpu():
    def build(st):
        m = channel_mesh(4.0, 1.0, 20, 10, periodic=False)
        return Solver2D(m, 0.05, {"inlet": {"type": "inlet", "flow_rate": 1.0,
                                            "scalars": {"c": "where(y < 0.5, 1.0, 0.0)"}},
                                  "outlet": {"type": "outlet"}, "bottom": {"type": "wall"},
                                  "top": {"type": "wall"}},
                        viscosity={"model": "power_law", "K": 0.05, "n": 0.6},
                        scalars={"c": {"diffusivity": 0.01}}, settings=st)
    c, g = _pair(build, lambda s: s.run_steady(max_iter=40, tol=1e-12))
    assert np.max(np.abs(c.U - g.U)) < 1e-8
    assert np.max(np.abs(c.scalars["c"] - g.scalars["c"])) < 1e-8
    assert np.allclose(c.nu_lam, g.nu_lam, rtol=1e-8, atol=0)


def test_porous_zone_on_gpu():
    def build(st):
        m = channel_mesh(4.0, 1.0, 20, 6, periodic=False)
        return Solver2D(m, 0.05, {"inlet": {"type": "inlet", "U": [0.5, 0.0]},
                                  "outlet": {"type": "outlet"}, "bottom": {"type": "symmetry"},
                                  "top": {"type": "symmetry"}},
                        porous=[{"region": "rectangle", "x0": 1.5, "x1": 2.5, "y0": -1,
                                 "y1": 2, "darcy": [50.0, 200.0], "forchheimer": 3.0,
                                 "angle": 30.0}], settings=st)
    c, g = _pair(build, lambda s: s.run_steady(max_iter=40, tol=1e-12))
    assert np.max(np.abs(c.U - g.U)) < 1e-8 and np.max(np.abs(c.p - g.p)) < 1e-8


def test_generic_csr_matches_scipy():
    """Matrices creuses « génériques » (bincount, lectures indexées) utilisées avec dpnp :
    vérifiées ici avec NumPy comme module de tableaux."""
    import scipy.sparse as sp

    from microrans.sparse_generic import GenericSparseModule
    rng = np.random.default_rng(1)
    A = sp.random(40, 40, density=0.2, random_state=2, format="csr") + sp.eye(40)
    A.sort_indices()
    G = GenericSparseModule(np).csr_matrix((A.data, A.indices, A.indptr), shape=A.shape)
    x = rng.normal(size=40)
    assert np.allclose(G @ x, A @ x) and np.allclose(G.diagonal(), A.diagonal())
    assert np.allclose(G.toarray(), A.toarray())
    assert abs(G.to_scipy(np.asarray) - A).max() == 0


def test_intel_backend_matches_cpu_if_available():
    """Backend intel (dpnp) : exécuté seulement si dpnp et un matériel oneAPI en double
    précision sont présents (pas en intégration continue)."""
    pytest.importorskip("dpnp")
    try:
        get_backend("intel:cpu")
    except RuntimeError as exc:
        pytest.skip(str(exc))

    def build(st):
        return Solver2D(channel_mesh(4.0, 1.0, 12, 6, periodic=False), 0.05,
                        {"inlet": {"type": "inlet", "U": [1.0, 0.0]},
                         "outlet": {"type": "outlet"}, "bottom": {"type": "wall"},
                         "top": {"type": "wall"}}, settings=st)
    res = {}
    for be in ("cpu", "intel:cpu"):
        s = build(Settings(backend=be, solver_p="amg", solver_U="bicgstab"))
        if be == "cpu":
            s.fvm.lin.amg.smoother = "chebyshev"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s.run_steady(max_iter=10, tol=1e-30)
        res[be] = s.to_cpu()
    assert np.max(np.abs(res["cpu"].U - res["intel:cpu"].U)) < 1e-10


def test_coupled_rejected_on_gpu():
    m = cavity_mesh(8)
    with pytest.raises(ValueError, match="CPU seulement"):
        Solver2D(m, 0.01, {"lid": {"type": "wall", "U": [1, 0]}, "walls": {"type": "wall"}},
                 settings=Settings(backend="fakegpu", algorithm="coupled"))


def test_gpu_backends_explained_in_executable(monkeypatch):
    """Audit M6 : dans l'exécutable (PyInstaller), les cartes graphiques ne peuvent pas
    fonctionner ; le message ne conseille plus un « pip install » impossible."""
    import sys

    from microrans.backend import available
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    assert available("cpu") and not available("cuda") and not available("intel")
    for name in ("cuda", "rocm", "intel"):
        with pytest.raises(RuntimeError, match="exécutable calcule uniquement sur le processeur"):
            get_backend(name)
    with pytest.raises(ValueError, match="Backend inconnu"):
        get_backend("gpu2")
