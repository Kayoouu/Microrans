"""Mailleur 2D : géométrie, générateurs, cohérence topologique, formats."""
import warnings

import numpy as np
import pytest

from microrans.cli import main
from microrans.mesh2d import (NACA4, Circle, Rectangle, Spline, backward_facing_step_mesh,
                              block_mesh, cavity_mesh, channel_mesh, grading_distribution,
                              hybrid_mesh, o_grid, read_curve, read_mesh, rectangle_mesh,
                              shape_from_dict, triangle_quality, triangulate, write_mesh)
from microrans.mesh2d.io import write_openfoam


def closed_cells(mesh):
    """Σ S_f sortant de chaque cellule = 0 (cellules fermées, faces bien orientées)."""
    s = np.zeros((mesh.n_cells, 2))
    ni = mesh.n_internal
    for k in range(2):
        s[:, k] += np.bincount(mesh.owner, mesh.Sf[:, k], mesh.n_cells)
        s[:, k] -= np.bincount(mesh.neighbour, mesh.Sf[:ni, k], mesh.n_cells)
    return np.max(np.abs(s)) / mesh.magSf.max()


# ------------------------------------------------------------------ géométrie
def test_sdf_signs_and_csg():
    c = Circle((0, 0), 1.0)
    r = Rectangle(-2, -2, 2, 2)
    pts = np.array([[0, 0], [1.5, 0], [3, 0]])
    assert np.allclose(c.sdf(pts), [-1, 0.5, 2])
    assert np.allclose(r.sdf(pts), [-2, -0.5, 1])
    d = (r - c).sdf(pts)
    assert d[0] > 0 and d[1] < 0 and d[2] > 0


def test_polygon_fast_sdf_matches_brute_force():
    poly = NACA4("2412")
    rng = np.random.default_rng(1)
    p = rng.uniform([-0.5, -0.5], [1.5, 0.5], (500, 2))
    assert np.allclose(poly._sdf_fast(p), poly._sdf_brute(p), atol=1e-12)


def test_naca_is_closed_ccw_and_thickness():
    pts = NACA4("0012").points
    assert pts[:, 1].max() == pytest.approx(0.06, rel=0.01)   # 12 % → demi-épaisseur 0.06
    x, y = pts[:, 0], pts[:, 1]
    assert 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) > 0


def test_transform_incidence_nose_up():
    s = shape_from_dict({"type": "naca", "code": "0012", "incidence": 10.0,
                         "rotation_center": [0.25, 0.0]})
    curve = s.boundary_curve(n=100)
    le = curve[np.argmin(curve[:, 0])]
    te = curve[np.argmax(curve[:, 0])]
    assert le[1] > te[1]


def test_read_curve_formats(tmp_path):
    pts = NACA4("0012", n=41).points
    (tmp_path / "selig.dat").write_text("NACA0012\n" + "\n".join(f"{x} {y}" for x, y in pts))
    up = pts[: len(pts) // 2 + 1][::-1]
    lo = pts[len(pts) // 2:]
    (tmp_path / "led.dat").write_text(f"NACA0012\n{len(up)}. {len(lo)}.\n\n"
                                      + "\n".join(f"{x} {y}" for x, y in up) + "\n\n"
                                      + "\n".join(f"{x} {y}" for x, y in lo))
    (tmp_path / "c.csv").write_text("x,y\n" + "\n".join(f"{x},{y}" for x, y in pts))
    (tmp_path / "s.svg").write_text('<svg><path d="M 0 0 L 10 0 L 10 -5 Z"/></svg>')
    (tmp_path / "d.dxf").write_text("0\nSECTION\n2\nENTITIES\n"
                                    "0\nLINE\n10\n0\n20\n0\n11\n1\n21\n0\n"
                                    "0\nLINE\n10\n1\n20\n0\n11\n1\n21\n1\n"
                                    "0\nLINE\n10\n0\n20\n1\n11\n1\n21\n1\n"
                                    "0\nLINE\n10\n0\n20\n1\n11\n0\n21\n0\n0\nENDSEC\n0\nEOF\n")
    area = lambda c: abs(0.5 * np.sum(c[:, 0] * np.roll(c[:, 1], -1) - np.roll(c[:, 0], -1) * c[:, 1]))
    a_ref = area(pts)
    assert area(read_curve(tmp_path / "selig.dat")) == pytest.approx(a_ref, rel=1e-9)
    assert area(read_curve(tmp_path / "led.dat")) == pytest.approx(a_ref, rel=1e-9)
    assert area(read_curve(tmp_path / "c.csv")) == pytest.approx(a_ref, rel=1e-9)
    assert area(read_curve(tmp_path / "s.svg")) == pytest.approx(25.0)
    assert area(read_curve(tmp_path / "d.dxf")) == pytest.approx(1.0)


# ------------------------------------------------------------------ structurés
def test_grading_distribution():
    s = grading_distribution(10, 5.0)
    d = np.diff(s)
    assert s[0] == 0 and s[-1] == 1
    assert d[-1] / d[0] == pytest.approx(5.0)
    s2 = grading_distribution(20, [(0.5, 0.5, 4.0), (0.5, 0.5, 0.25)])
    assert np.allclose(np.diff(s2), np.diff(s2)[::-1])


def test_rectangle_and_multiblock_areas():
    m = rectangle_mesh(0, 2, 0, 1, 10, 5, grading=(2.0, 3.0))
    assert m.cell_volumes.sum() == pytest.approx(2.0)
    assert closed_cells(m) < 1e-12
    bfs = backward_facing_step_mesh(step=1, upstream=4, downstream=30, channel=1)
    assert bfs.cell_volumes.sum() == pytest.approx(4 * 1 + 30 * 2)
    assert {p.name for p in bfs.patches} == {"inlet", "outlet", "walls"}
    assert closed_cells(bfs) < 1e-12


def test_curved_block_edges_quarter_annulus():
    V = [[1, 0], [2, 0], [0, 2], [0, 1]]
    edges = [{"type": "arc", "vertices": [1, 2], "point": [2 ** 0.5, 2 ** 0.5]},
             {"type": "arc", "vertices": [0, 3], "point": [0.5 ** 0.5, 0.5 ** 0.5]}]
    m = block_mesh(V, [{"vertices": [0, 1, 2, 3], "cells": [20, 40]}], edges,
                   {"inner": {"type": "wall", "faces": [[3, 0]]},
                    "other": {"type": "patch", "faces": [[0, 1], [1, 2], [2, 3]]}})
    assert m.cell_volumes.sum() == pytest.approx(np.pi * (4 - 1) / 4, rel=2e-3)


def test_periodic_channel_geometry():
    m = channel_mesh(1.0, 2.0, 4, 32, first_height=1e-2)
    assert [p.name for p in m.patches] == ["bottom", "top"]
    assert len(m.periodic_pairs) == 1
    assert closed_cells(m) < 1e-12
    # les faces périodiques relient bien des cellules voisines « à travers » la frontière
    per = np.any(m.shift != 0, axis=1)
    assert np.allclose(np.abs(m.d_PN[per, 0]), 0.25)
    assert m.first_cell_height == pytest.approx(5e-3, rel=1e-6)


def test_ogrid_circle():
    m = o_grid(Circle((0, 0), 0.5, "cyl"), 64, 32, 10.0, 1e-3)
    assert m.patch("cyl").type == "wall"
    assert m.cell_volumes.sum() == pytest.approx(np.pi * (100 - 0.25), rel=2e-3)
    assert m.quality()["non_orthogonality_max_deg"] < 1.0
    assert closed_cells(m) < 1e-10


def test_ogrid_airfoil_valid():
    m = o_grid(NACA4("2412", name="wing"), 128, 48, 15.0, 1e-5)
    assert np.all(m.cell_volumes > 0)
    assert m.first_cell_height < 1e-5


# ------------------------------------------------------------------ non structurés
def test_triangulate_cylinder_in_box():
    dom = Rectangle(-3, -3, 6, 3, names={"left": "inlet", "right": "outlet",
                                         "bottom": "sym", "top": "sym"})
    cyl = Circle((0, 0), 0.5, "cylinder").as_wall()
    m = triangulate(dom - cyl, 0.6, [{"shape": cyl, "h": 0.1, "growth": 0.2}])
    assert {p.name for p in m.patches} == {"inlet", "outlet", "sym", "cylinder"}
    assert m.patch("cylinder").type == "wall"
    assert m.cell_volumes.sum() == pytest.approx(54 - np.pi * 0.25, rel=5e-3)
    assert triangle_quality(m).min() > 0.3
    assert closed_cells(m) < 1e-10


def test_hybrid_mesh_is_conformal():
    dom = Rectangle(-3, -3, 6, 3, names={"left": "inlet", "right": "outlet",
                                         "bottom": "sym", "top": "sym"})
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # un raccord non conforme émet un warning
        m = hybrid_mesh(dom, [Circle((0, 0), 0.5, "cylinder").as_wall()], 0.6, 0.1,
                        n_layers=5, first_height=1e-2, ratio=1.2, max_iter=150)
    q = m.quality()
    assert set(q["cell_types"]) == {"triangles", "quadrilatères"}
    assert not any(p.name.startswith("_") for p in m.patches)
    assert m.first_cell_height == pytest.approx(5e-3, rel=0.05)
    assert m.cell_volumes.sum() == pytest.approx(54 - np.pi * 0.25, rel=5e-3)


def test_spline_body():
    s = Spline([[0, 0], [1, 0.3], [2, 0], [1, -0.3]], name="blob")
    assert s.sdf(np.array([[1.0, 0.0]]))[0] < 0


# ------------------------------------------------------------------ formats
@pytest.mark.parametrize("ext", ["msh", "su2"])
def test_mesh_roundtrip(tmp_path, ext):
    m = cavity_mesh(8)
    write_mesh(m, tmp_path / f"m.{ext}")
    r = read_mesh(tmp_path / f"m.{ext}", {"walls": "wall", "lid": "wall"})
    assert r.n_cells == m.n_cells
    assert np.allclose(np.sort(r.cell_volumes), np.sort(m.cell_volumes))
    assert r.quality()["patches"] == m.quality()["patches"]


def test_vtk_and_openfoam_export(tmp_path):
    m = channel_mesh(1.0, 2.0, 3, 8)
    write_mesh(m, tmp_path / "m.vtk", {"a": m.cell_volumes, "c": m.cell_centers})
    assert "CELL_TYPES 24" in (tmp_path / "m.vtk").read_text()
    d = write_openfoam(m, tmp_path / "case")
    b = (d / "boundary").read_text()
    assert "cyclic" in b and "neighbourPatch  outlet" in b and "empty" in b
    nfaces = int((d / "faces").read_text().split("\n\n", 1)[1].split("\n")[0])
    # faces internes régulières + parois + 2 patches cycliques + avant/arrière
    assert nfaces == m.n_regular_internal + 2 * 3 + 2 * 8 + 2 * m.n_cells


def test_cli_mesh_preset(tmp_path):
    assert main(["mesh", "--preset", "cavity", "-o", str(tmp_path), "-f", "msh", "su2", "-q"]) == 0
    assert (tmp_path / "cavity.msh").exists() and (tmp_path / "cavity.png").exists()
    assert (tmp_path / "quality.json").exists()


def test_shape_errors_are_explained(tmp_path):
    """Audit M1/M7/C5 : corps sans type, paramètre manquant, dimensions impossibles,
    contour sans fichier : messages clairs au lieu de KeyError / IsADirectoryError ou d'un
    maillage absurde accepté en silence."""
    import pytest

    from microrans.mesh2d.builder import build_mesh
    from microrans.mesh2d.geometry import shape_from_dict
    bad = [({"radius": 0.5}, "sans « type »"),
           ({"type": "circle"}, "manquant.*radius"),
           ({"type": "circle", "radius": -0.5}, "radius doit être > 0"),
           ({"type": "ellipse", "a": 1.0, "b": 0.0}, "b doit être > 0"),
           ({"type": "rectangle", "x0": 1, "x1": 0, "y0": 0, "y1": 1}, "x0 < x1"),
           ({"type": "file", "path": ""}, "manquant.*path"),
           ({"type": "file", "path": "absent.dat"}, "introuvable"),
           ({"type": "carre"}, "inconnu")]
    for spec, msg in bad:
        with pytest.raises(ValueError, match=msg):
            shape_from_dict(spec, tmp_path)
    # [domain] sans type : rectangle par défaut (interface)
    cfg = {"mesh": {"type": "unstructured", "h_max": 4.0, "h_surface": 0.5},
           "domain": {"x0": -4.0, "x1": 8.0, "y0": -4.0, "y1": 4.0},
           "bodies": [{"type": "circle", "radius": 0.5, "name": "c"}]}
    assert build_mesh(cfg).n_cells > 0
