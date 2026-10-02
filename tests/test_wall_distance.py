"""Distance à la paroi (lot E2) : la recherche accélérée (arbre k-d + arbre de boîtes + bornes
inférieures) donne exactement le résultat de la comparaison à toutes les faces de paroi, en
2D et en 3D, y compris sur les cas qui la mettent en défaut : faces gauches, faces de tailles
très différentes sur une même paroi, prismes, maillages très étirés."""
import warnings

import numpy as np

from microrans.mesh2d import Circle, Rectangle, bvh, o_grid, rectangle_mesh, triangulate
from microrans.mesh2d.mesh import Mesh2D
from microrans.mesh3d import Mesh3D, box_mesh, extrude


def _brute_3d(m, names):
    """Toutes les faces pour toutes les cellules, dans l'ordre (égalité : premier numéro)."""
    fn = np.vstack([m.patch_face_nodes(n) for n in names])
    kk = np.sum(fn >= 0, axis=1)
    X = m.points[np.where(fn >= 0, fn, fn[:, :1])]
    valid = np.arange(4)[None, :] < kk[:, None]
    xbar = np.sum(np.where(valid[..., None], X, 0.0), axis=1) / kk[:, None]
    cand = np.broadcast_to(np.arange(len(fn)), (m.n_cells, len(fn)))
    return Mesh3D._dist_to_faces(m.cell_centers, cand, X, xbar, kk)


def _brute_2d(m, names):
    """L'ancien calcul (tous les segments, argmin) : référence au bit près."""
    seg = np.vstack([m.face_nodes[m.patch(n).faces] for n in names])
    a, b = m.points[seg[:, 0]], m.points[seg[:, 1]]
    ab = b - a
    ab2 = np.maximum(np.sum(ab * ab, axis=1), 1e-300)
    ap = m.cell_centers[:, None, :] - a[None]
    t = np.clip(np.sum(ap * ab[None], axis=2) / ab2[None], 0.0, 1.0)
    r = ap - t[..., None] * ab[None]
    d2 = np.sum(r ** 2, axis=2)
    j = d2.argmin(axis=1)
    rows = np.arange(len(j))
    return np.sqrt(d2[rows, j]), r[rows, j]


def _check_3d(m):
    d, v = m._nearest_on_patches(m.wall_patches)
    d0, v0 = _brute_3d(m, m.wall_patches)
    assert np.array_equal(d, d0)
    assert np.allclose(v, v0, rtol=0, atol=1e-14 * np.abs(v0).max())


def test_3d_warped_wall_faces():
    """Sommets déplacés au hasard : faces de paroi gauches (écart au plan moyen > 0)."""
    walls = {n: "wall" for n in ("left", "right", "bottom", "top", "back", "front")}
    m = box_mesh(0, 1, 0, 1, 0, 1, 7, 6, 5, types=walls)
    m.points = m.points + np.random.default_rng(1).uniform(-0.04, 0.04, m.points.shape)
    _check_3d(m)


def test_3d_wall_faces_of_very_different_sizes():
    """Rapport 200 entre la première et la dernière maille sur les parois : la recherche par
    sphère de rayon « d + plus grand rayon de face » examinait presque toutes les petites
    faces (64 000 cellules : 207 s avant, 2.3 s après)."""
    m = box_mesh(0, 10, 0, 1, 0, 2, 30, 6, 8, grading=(200, 1, 30),
                 types={"bottom": "wall", "back": "wall"})
    _check_3d(m)


def test_3d_prisms_and_curved_wall():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m2 = triangulate(Rectangle(0, 0, 2, 1) - Circle((0.7, 0.5), 0.25), h_max=0.15)
    m = extrude(m2, 0, 0.3, 2, types={"circle": "wall", "back": "wall"})
    assert m.quality()["cell_types"] == {"prismes": m.n_cells}
    _check_3d(m)


def test_2d_identical_to_all_segments():
    """Ancien calcul : chaque cellule contre tous les segments (490 000 cellules, 2 800
    segments : 91 s ; nouveau : 2.7 s). Même résultat au bit près, égalités comprises
    (segment de plus petit numéro), sur un rectangle à 4 parois (égalités sur les
    diagonales), un maillage très resserré et un maillage en O étiré."""
    m = rectangle_mesh(0, 1, 0, 1, 40, 40)
    m.set_patch_types({n: "wall" for n in ("left", "right", "bottom", "top")})
    meshes = [m, rectangle_mesh(0, 30, 0, 1, 60, 30, grading=(200, 50),
                                types={"bottom": "wall", "top": "wall"}),
              o_grid(Circle((0, 0), 0.5), n_around=64, n_radial=40, farfield_radius=50,
                     first_height=1e-3)]
    for m in meshes:
        d, v = Mesh2D._nearest_on_patches(m, m.wall_patches)
        d0, v0 = _brute_2d(m, m.wall_patches)
        assert np.array_equal(d, d0) and np.array_equal(v, v0)


def test_bvh_structure():
    """Chaque face dans une seule feuille, boîtes des nœuds contenant celles des faces."""
    rng = np.random.default_rng(0)
    lo = rng.uniform(0, 1, (1000, 3))
    hi = lo + rng.uniform(0, 0.05, (1000, 3))
    order, start, end, left, right, blo, bhi = bvh.build(lo, hi, leaf=8)
    leaves = np.nonzero(left < 0)[0]
    assert np.array_equal(np.sort(np.concatenate([order[start[i]:end[i]] for i in leaves])),
                          np.arange(1000))
    assert (end[leaves] - start[leaves]).max() <= 8
    for i in range(len(start)):
        f = order[start[i]:end[i]]
        assert np.all(blo[i] <= lo[f]) and np.all(bhi[i] >= hi[f])
