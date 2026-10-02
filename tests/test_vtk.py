"""Sortie VTK (lot E3) : binaire par défaut (relecture exacte, au bit près), texte en option
(format d'avant, inchangé), choix par [output] vtk_format, lecture par le lecteur officiel
VTK quand il est installé (pas une dépendance : test ignoré sinon)."""
import warnings

import numpy as np
import pytest

from microrans.fv2d.validate import check_case
from microrans.mesh2d import Circle, Rectangle, rectangle_mesh, triangulate
from microrans.mesh2d.io import read_vtk, write_mesh, write_vtk
from microrans.mesh2d.mesh import Mesh2D
from microrans.mesh3d import box_mesh, extrude


def _meshes():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tri = triangulate(Rectangle(0, 0, 2, 1) - Circle((0.7, 0.5), 0.25), h_max=0.2)
        pts = np.array([[0, 0], [1, 0], [2, .3], [1.5, 1], [.5, 1.2], [-.5, .6], [3, 1],
                        [2.5, 2]], float)
        poly = Mesh2D(pts, [[0, 1, 2, 3, 4], [0, 4, 5], [2, 6, 7, 3]])
    return {"quads": (rectangle_mesh(0, 1, 0, 1, 5, 3), {9}),
            "triangles": (tri, {5}),
            "polygones": (poly, {5, 7, 9}),
            "hexaèdres": (box_mesh(0, 1, 0, 2, 0, 3, 3, 4, 2), {12}),
            "prismes": (extrude(tri, 0, 0.5, 2), {13})}


def _rows(m):
    return m.cells_as_lists() if getattr(m, "dim", 2) == 3 else [
        r[:k] for r, k in zip(m.cell_nodes, m.cell_nv)]


def _flat(m):
    return np.concatenate([[len(r), *r] for r in _rows(m)])


@pytest.mark.parametrize("name", ["quads", "triangles", "polygones", "hexaèdres", "prismes"])
def test_binary_roundtrip_is_exact(tmp_path, name):
    m, types = _meshes()[name]
    rng = np.random.default_rng(0)
    d = m.points.shape[1]
    data = {"p": rng.normal(size=m.n_cells), "U": rng.normal(size=(m.n_cells, d)),
            "nu t": rng.normal(size=m.n_cells) * 1e-300, "température": np.arange(m.n_cells)}
    write_vtk(m, tmp_path / "f.vtk", data)
    r = read_vtk(tmp_path / "f.vtk")
    assert r["binary"] and r["title"] == "microrans"
    P = m.points if d == 3 else np.column_stack([m.points, np.zeros(m.n_points)])
    assert np.array_equal(r["points"], P)
    assert np.array_equal(r["cells"], _flat(m))
    assert set(r["cell_types"].tolist()) == types and len(r["cell_types"]) == m.n_cells
    assert np.array_equal(r["cell_data"]["p"], data["p"])
    assert np.array_equal(r["cell_data"]["nu_t"], data["nu t"])     # nom nettoyé
    assert np.array_equal(r["cell_data"]["température"], data["température"])
    U3 = data["U"] if d == 3 else np.column_stack([data["U"], np.zeros(m.n_cells)])
    assert np.array_equal(r["cell_data"]["U"], U3)
    # le texte relu donne la même chose à 10 chiffres près
    write_vtk(m, tmp_path / "t.vtk", data, binary=False)
    t = read_vtk(tmp_path / "t.vtk")
    assert not t["binary"] and np.array_equal(t["cells"], r["cells"])
    for k in ("p", "U"):
        assert np.allclose(t["cell_data"][k], r["cell_data"][k], rtol=1e-9, atol=0)


def test_ascii_option_keeps_the_former_text_format(tmp_path):
    m = rectangle_mesh(0, 2, 0, 1, 2, 1)
    write_mesh(m, tmp_path / "m.vtk", {"p": [1.5, -2.0], "U": [[1, 0], [0.25, 1 / 3]]},
               binary=False)
    assert (tmp_path / "m.vtk").read_text().splitlines() == [
        "# vtk DataFile Version 3.0", "microrans", "ASCII", "DATASET UNSTRUCTURED_GRID",
        "POINTS 6 double", "0 0 0", "0 1 0", "1 0 0", "1 1 0", "2 0 0", "2 1 0",
        "CELLS 2 10", "4 0 2 3 1", "4 2 4 5 3", "CELL_TYPES 2", "9", "9",
        "CELL_DATA 2", "SCALARS p double 1", "LOOKUP_TABLE default", "1.5", "-2",
        "VECTORS U double", "1 0 0", "0.25 0.3333333333 0"]
    head = (tmp_path / "b.vtk")
    write_mesh(m, head)
    assert head.read_bytes().split(b"\n")[:5] == [
        b"# vtk DataFile Version 3.0", b"microrans", b"BINARY", b"DATASET UNSTRUCTURED_GRID",
        b"POINTS 6 double"]


def test_vtk_format_key_in_case_file(tmp_path):
    from microrans.fv2d.case import run_case
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0, "y1": 1, "nx": 4,
                    "ny": 4, "names": {"left": "w", "right": "w", "bottom": "w", "top": "lid"}},
           "physics": {"nu": 0.1},
           "boundary": {"lid": {"type": "wall", "U": [1, 0]}, "w": {"type": "wall"}},
           "solver": {"max_iter": 3},
           "output": {"plots": False, "checkpoint": False, "vtk_format": "texte"}}
    with pytest.raises(ValueError, match="vtk_format = 'texte' inconnu"):
        check_case(cfg)
    for fmt, binary in (("ascii", False), ("ASCII", False), ("binary", True)):
        cfg["output"]["vtk_format"] = fmt
        check_case(cfg)
        run_case(cfg, out_dir=tmp_path / fmt, verbose=False, plot=False)
        r = read_vtk(tmp_path / fmt / "fields.vtk")
        assert r["binary"] is binary and r["cell_data"]["U"].shape == (16, 3)
    del cfg["output"]["vtk_format"]
    run_case(cfg, out_dir=tmp_path / "defaut", verbose=False, plot=False)
    assert read_vtk(tmp_path / "defaut" / "fields.vtk")["binary"]


def test_official_vtk_reader(tmp_path):
    """Le lecteur de ParaView (vtkUnstructuredGridReader) relit nos fichiers binaires à
    l'identique. VTK n'est pas une dépendance : ignoré s'il n'est pas installé."""
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import vtk_to_numpy
    for name, (m, types) in _meshes().items():
        U = np.random.default_rng(1).normal(size=(m.n_cells, m.points.shape[1]))
        write_vtk(m, tmp_path / "f.vtk", {"U": U, "p": U[:, 0]})
        rd = vtk.vtkUnstructuredGridReader()
        rd.SetFileName(str(tmp_path / "f.vtk"))
        rd.ReadAllScalarsOn()
        rd.ReadAllVectorsOn()
        rd.Update()
        g = rd.GetOutput()
        assert (g.GetNumberOfPoints(), g.GetNumberOfCells()) == (m.n_points, m.n_cells), name
        assert {g.GetCellType(i) for i in range(m.n_cells)} == types, name
        r = read_vtk(tmp_path / "f.vtk")
        assert np.array_equal(vtk_to_numpy(g.GetCellData().GetArray("U")), r["cell_data"]["U"])
        assert np.array_equal(vtk_to_numpy(g.GetCellData().GetArray("p")), U[:, 0])
        ids = vtk.vtkIdList()
        for i in (0, m.n_cells - 1):
            g.GetCellPoints(i, ids)
            row = [ids.GetId(k) for k in range(ids.GetNumberOfIds())]
            assert row == list(_rows(m)[i]), name
