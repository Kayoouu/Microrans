"""Maillage 3D (mesh3d) : géométrie des cellules et des faces, pavé, extrusion d'un maillage
2D, périodicité, distance à la paroi exacte."""
import numpy as np
import pytest

from microrans.mesh2d import channel_mesh, o_grid
from microrans.mesh2d.geometry import Circle
from microrans.mesh3d import Mesh3D, box_mesh, extrude


def _closure(m):
    """Σ_f S_f par cellule (nul pour une cellule fermée)."""
    ni = m.n_internal
    P, N, Pb = m.owner[:ni], m.neighbour, m.owner[ni:]
    return np.column_stack([np.bincount(P, m.Sf[:ni, k], m.n_cells)
                            - np.bincount(N, m.Sf[:ni, k], m.n_cells)
                            + np.bincount(Pb, m.Sf[ni:, k], m.n_cells) for k in range(3)])


def test_box_geometry():
    m = box_mesh(0, 1, 0, 2, 0, 3, 4, 5, 6, grading=(1, 2, 0.5))
    assert m.dim == 3 and m.n_cells == 120
    assert m.cell_volumes.sum() == pytest.approx(6.0, rel=1e-13)
    assert np.abs(_closure(m)).max() < 1e-14
    assert {p.name: p.size for p in m.patches} == {"left": 30, "right": 30, "bottom": 24,
                                                   "top": 24, "back": 20, "front": 20}
    q = m.quality()
    assert q["cell_types"] == {"hexaèdres": 120}
    assert q["non_orthogonality_max_deg"] < 1e-10 and q["skewness_max"] < 1e-10


def test_distorted_hexahedra_closed_and_consistent():
    """Sommets déplacés au hasard (faces gauches) : cellules fermées, volume total exact
    (théorème de la divergence : V = Σ x_f·S_f / 3 sur la frontière)."""
    m0 = box_mesh(0, 1, 0, 1, 0, 1, 5, 5, 5)
    pts = m0.points.copy()
    inner = np.all((pts > 1e-9) & (pts < 1 - 1e-9), axis=1)
    pts[inner] += np.random.default_rng(0).uniform(-0.04, 0.04, (inner.sum(), 3))
    m = Mesh3D(pts, np.array(m0.cells_as_lists()),
               {p.name: m0.face_nodes[p.faces] for p in m0.patches})
    assert np.abs(_closure(m)).max() < 1e-14
    assert m.cell_volumes.sum() == pytest.approx(1.0, rel=1e-12)
    ni = m.n_internal
    assert np.sum(m.face_centers[ni:] * m.Sf[ni:]) / 3.0 == pytest.approx(1.0, rel=1e-12)
    assert m.quality()["non_orthogonality_max_deg"] > 1.0


def test_extrusion_keeps_2d_geometry_and_patches():
    m2 = o_grid(Circle((0.0, 0.0), 0.5), n_around=32, n_radial=16, farfield_radius=10.0,
                first_height=0.01, wall_name="cyl")
    m3 = extrude(m2, 0.0, 2.0, 4, types={"back": "symmetry", "front": "symmetry"})
    assert m3.n_cells == 4 * m2.n_cells
    assert m3.cell_volumes.sum() == pytest.approx(2.0 * m2.cell_volumes.sum(), rel=1e-12)
    assert np.abs(m3.cell_centers[:m2.n_cells, :2] - m2.cell_centers).max() < 1e-12
    types = {p.name: p.type for p in m3.patches}
    assert types == {"cyl": "wall", "farfield": m2.patch("farfield").type,
                     "back": "symmetry", "front": "symmetry"}
    # distance à la paroi : identique au 2D (exacte dans les deux cas)
    assert np.abs(m3.wall_distance[:m2.n_cells] - m2.wall_distance).max() < 1e-12


def test_extrusion_of_triangles_and_quads():
    pts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1.0]])
    from microrans.mesh2d import Mesh2D
    m2 = Mesh2D(pts, [[0, 1, 4, 3], [1, 2, 4], [2, 5, 4]],
                {"walls": [[0, 1], [1, 2], [2, 5], [5, 4], [4, 3], [3, 0]]},
                patch_types={"walls": "wall"})
    m3 = extrude(m2, 0, 1, 2)
    assert m3.quality()["cell_types"] == {"hexaèdres": 2, "prismes": 4}
    assert m3.cell_volumes.sum() == pytest.approx(2.0, rel=1e-13)
    assert np.abs(_closure(m3)).max() < 1e-14
    with pytest.raises(ValueError, match="cellule 2D à 5 sommets"):
        extrude(Mesh2D(np.array([[0, 0], [1, 0], [1.5, .5], [1, 1], [0, 1.0]]),
                       [[0, 1, 2, 3, 4]], {"w": [[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]]}),
                0, 1, 1)
    with pytest.raises(ValueError, match="nz = 0"):
        extrude(m2, 0, 1, 0)


def test_periodic_box_and_extruded_channel():
    m = box_mesh(0, 2, 0, 1, 0, 1, 8, 4, 4, names={"left": "inlet", "right": "outlet"},
                 periodic=[("inlet", "outlet")])
    assert [p.name for p in m.patches] == ["bottom", "top", "back", "front"]
    assert m.n_internal - m.n_regular_internal == 16
    assert np.allclose(m.shift[m.n_regular_internal:], [2.0, 0.0, 0.0])
    m3 = extrude(channel_mesh(nx=4, ny=16), 0, 1, 3, names={"back": "z0", "front": "z1"},
                 periodic=[("z0", "z1")])
    assert [p[:2] for p in m3.periodic_pairs] == [("inlet", "outlet"), ("z0", "z1")]
    assert [p.name for p in m3.patches] == ["bottom", "top"]


def test_wall_distance_exact_in_square_duct():
    m = box_mesh(0, 4, 0, 1, 0, 1, 16, 8, 8,
                 types={"bottom": "wall", "top": "wall", "back": "wall", "front": "wall"})
    C = m.cell_centers
    exact = np.minimum.reduce([C[:, 1], 1 - C[:, 1], C[:, 2], 1 - C[:, 2]])
    assert np.abs(m.wall_distance - exact).max() < 1e-14


def test_wall_distance_exact_with_very_long_wall_faces():
    """Faces de paroi 100 fois plus longues que hautes : le plus proche centre de face n'est
    pas la face la plus proche ; la distance reste exacte."""
    m = box_mesh(0, 100, 0, 1, 0, 1, 2, 20, 1, types={"bottom": "wall"},
                 grading=(1.0, 1.0, 1.0))
    assert np.abs(m.wall_distance - m.cell_centers[:, 1]).max() < 1e-13
