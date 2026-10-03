"""Coupes planes des maillages 3D pour les figures (lot E4) : x = cte, y = cte (et z = cte)
sur hexaèdres, prismes, tétraèdres et pyramides. Contrôle : la somme des aires des
polygones de coupe est l'aire exacte de la section, chaque polygone est orienté (aire > 0),
la coupe z redonne ZSlice."""
import itertools
import warnings

import numpy as np
import pytest

from microrans.mesh2d import Circle, Rectangle, triangulate
from microrans.mesh3d import Mesh3D, box_mesh, extrude
from microrans.mesh3d.slice import PlaneSlice, ZSlice, slice_mesh


def _areas(sl):
    out = []
    for row, k in zip(sl.cell_nodes, sl.cell_nv):
        x, y = sl.points[row[:k]].T
        out.append(0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1)))
    return np.array(out)


def _check(sl, exact):
    a = _areas(sl)
    assert np.all(a > 0)                                  # sommets rangés, sans croisement
    assert a.sum() == pytest.approx(exact, rel=1e-12)
    assert np.all(np.diff(sl.cells) > 0)                  # cellules distinctes, triées
    assert sl.cell_centers.shape == (sl.n_cells, 2)


def test_box_cuts_through_cells_nodes_and_boundaries():
    m = box_mesh(0, 2, 0, 1, 0, 3, 8, 4, 6, grading=(3, 1, 0.5))
    for axis, values, exact, labels in (("x", (None, 0.37, 0.0, 2.0), 3.0, ("y", "z")),
                                        ("y", (0.5, 1.0, 0.25), 6.0, ("x", "z"))):
        for v in values:
            sl = PlaneSlice(m, axis, v)
            _check(sl, exact)
            assert sl.labels == labels and set(sl.cell_nv) == {4}
    # plan sur une face entre deux couches (y = 0.5, 4 couches) : cellules du dessus
    sl = PlaneSlice(m, "y", 0.5)
    assert np.all(m.cell_centers[sl.cells, 1] > 0.5) and sl.n_cells == 48
    with pytest.raises(ValueError, match="plan x = 5 hors du domaine"):
        PlaneSlice(m, "x", 5.0)


def test_extruded_mesh_exact_section_and_same_z_cut_as_zslice():
    """Aire exacte : épaisseur × longueur de la droite x = c (ou y = c) dans le maillage 2D,
    trou du cylindre compris."""
    m2 = triangulate(Rectangle(0, 0, 2, 1) - Circle((1, .5), .2).as_wall(), 0.1, [])
    m = extrude(m2, 0, 1.5, 3)

    def length(a, c):
        L = 0.0
        for row, k in zip(m2.cell_nodes, m2.cell_nv):
            q = m2.points[row[:k]]
            s = q[:, a] - c
            p = [q[i] + s[i] / (s[i] - s[(i + 1) % k]) * (q[(i + 1) % k] - q[i])
                 for i in range(k) if (s[i] > 0) != (s[(i + 1) % k] > 0)]
            L += np.linalg.norm(p[0] - p[1]) if len(p) == 2 else 0.0
        return L

    for axis, c in (("x", 1.0), ("x", 0.83), ("y", 0.5), ("y", 0.31)):
        _check(PlaneSlice(m, axis, c), 1.5 * length("xy".index(axis), c))
    for z in (None, 0.0, 0.5, 1.5):
        a, b = ZSlice(m, z), PlaneSlice(m, "z", z)
        assert np.array_equal(a.cells, b.cells)
        assert np.allclose(_areas(a), _areas(b), rtol=1e-12)
    assert isinstance(slice_mesh(m, "z"), ZSlice) and slice_mesh(m, "y").axis == "y"


def _kuhn_cube(n, jitter):
    """Cube unité : n³ cubes en 6 tétraèdres (Kuhn), nœuds intérieurs déplacés."""
    g = np.linspace(0, 1, n + 1)
    P = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    inner = np.all((P > 1e-12) & (P < 1 - 1e-12), axis=1)
    P[inner] += np.random.default_rng(0).uniform(-jitter, jitter, (inner.sum(), 3)) / n
    T = []
    for ijk in itertools.product(range(n), repeat=3):
        for perm in itertools.permutations(range(3)):
            c = np.array(ijk)
            path = [c.copy()]
            for ax in perm:
                c[ax] += 1
                path.append(c.copy())
            T.append([(i * (n + 1) + j) * (n + 1) + k for i, j, k in path])
    T = np.array(T)
    a, b, c, d = (P[T[:, q]] for q in range(4))
    neg = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) < 0
    T[neg] = T[neg][:, [0, 2, 1, 3]]
    return Mesh3D(P, list(T))


def test_tetrahedra_and_pyramids():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tm = _kuhn_cube(5, 0.25)
        P = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1],
                      [1, 1, 1], [0, 1, 1], [.5, .5, .5]], float)
        faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6],
                 [3, 0, 4, 7]]
        pm = Mesh3D(P, [f + [8] for f in faces])         # 6 pyramides de sommet le centre
    for axis, v in (("x", 0.5), ("y", 0.37), ("z", 0.21), ("x", 0.0), ("x", 1.0)):
        sl = PlaneSlice(tm, axis, v)
        _check(sl, 1.0)
        assert set(sl.cell_nv) <= {3, 4}
    for axis, v in (("x", 0.5), ("x", 0.2), ("y", 0.75)):   # x = 0.5 : par le sommet commun
        _check(PlaneSlice(pm, axis, v), 1.0)
