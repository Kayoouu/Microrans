"""Cas 3D décrits par le fichier de cas : [mesh] type = "box" et [mesh.extrude], sortie VTK
des hexaèdres et prismes, sondes et profils en 3D, reprise, vérifications avant calcul."""
import numpy as np
import pytest

from microrans.fv2d.case import build_solver, case_mesh
from microrans.fv2d.restart import load_checkpoint, save_checkpoint
from microrans.fv2d.sampling import Sampler, locate, parse_points
from microrans.fv2d.validate import check_case
from microrans.mesh2d import Circle, Rectangle
from microrans.mesh2d.io import write_mesh, write_vtk
from microrans.mesh2d.unstructured import triangulate
from microrans.mesh3d import box_mesh, extrude

SIDES = ("left", "right", "bottom", "top", "back", "front")


def _box_cfg(**mesh):
    return {"mesh": {"type": "box", "x0": 0, "x1": 1, "y0": 0, "y1": 2, "z0": 0, "z1": 3,
                     "nx": 3, "ny": 4, "nz": 5, **mesh},
            "physics": {"nu": 0.1},
            "boundary": {k: {"type": "wall"} for k in SIDES}}


def test_box_and_extrude_from_case_file():
    m = case_mesh(_box_cfg(grading=[1, 2, 0.5], names={"back": "z0"},
                           patch_types={"z0": "symmetry"}))
    assert m.dim == 3 and m.n_cells == 60
    assert m.cell_volumes.sum() == pytest.approx(6.0, rel=1e-13)
    assert {p.name: p.type for p in m.patches}["z0"] == "symmetry"
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 2, "y0": 0, "y1": 1, "nx": 4, "ny": 2,
                    "names": {"left": "in", "right": "out"}, "periodic": [["in", "out"]],
                    "extrude": {"z0": -1, "z1": 1, "nz": 3, "names": {"front": "f"},
                                "periodic": [["back", "f"]]}}}
    m = case_mesh(cfg)
    assert m.n_cells == 24 and m.cell_volumes.sum() == pytest.approx(4.0, rel=1e-13)
    assert [p[:2] for p in m.periodic_pairs] == [("in", "out"), ("back", "f")]
    m = case_mesh({"mesh": {"preset": "cavity", "extrude": {"nz": 2}}})
    assert m.dim == 3 and m.n_cells == 2 * 64 * 64


def _read_vtk_cells(path):
    lines = path.read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.startswith("CELLS"))
    n = int(lines[i].split()[1])
    cells = [[int(v) for v in ln.split()[1:]] for ln in lines[i + 1:i + 1 + n]]
    types = [int(v) for v in lines[i + 2 + n:i + 2 + 2 * n]]
    return cells, types


def test_vtk_output_of_hexahedra_and_prisms(tmp_path):
    """Numérotation VTK des sommets (volumes recalculés par la bibliothèque VTK 9.7 égaux
    aux nôtres à 1e-8, précision du fichier ASCII : vérifié hors des tests, VTK n'étant pas
    une dépendance)."""
    m2 = triangulate(Rectangle(0, 0, 2, 1) - Circle((1, .5), .2).as_wall(), 0.2, [])
    for m in (box_mesh(0, 1, 0, 1, 0, 1, 2, 3, 4), extrude(m2, 0, 1, 2)):
        U = m.cell_centers.copy()
        write_vtk(m, tmp_path / "f.vtk", {"U": U, "p": m.cell_volumes})
        cells, types = _read_vtk_cells(tmp_path / "f.vtk")
        assert [list(c) for c in m.cells_as_lists()] == cells
        assert set(types) == ({12} if m.n_cells == 24 else {13})
        txt = (tmp_path / "f.vtk").read_text()
        assert f"POINTS {m.n_points} double" in txt and "VECTORS U double" in txt
    with pytest.raises(ValueError, match="Maillage 3D : export .su2 non disponible"):
        write_mesh(m, tmp_path / "f.su2")


def test_probes_and_location_in_3d():
    m2 = triangulate(Rectangle(0, 0, 2, 1) - Circle((1, .5), .2).as_wall(), 0.1, [])
    m = extrude(m2, 0, 0.5, 3)
    X = np.random.default_rng(1).uniform([0, 0, 0], [2, 1, 0.5], (2000, 3))
    c = locate(m, X)
    r = np.hypot(X[:, 0] - 1, X[:, 1] - 0.5)
    assert np.all(c[r < 0.19] < 0) and np.all(c[r > 0.21] >= 0)    # trou : hors domaine
    assert np.all(locate(m, m.points) >= 0) and np.all(locate(m, m.face_centers) >= 0)
    # champ linéaire : reconstruction exacte dans les cellules intérieures
    mb = box_mesh(0, 1, 0, 2, 0, 3, 6, 7, 8, grading=(1, 2, 0.5))
    S = build_solver({**_box_cfg(), "mesh": {}}, mesh=mb)
    S.p = 1 + mb.cell_centers @ [1, 2, 3]
    X = np.random.default_rng(2).uniform([0, 0, 0], [1, 2, 3], (2000, 3))
    border = np.zeros(mb.n_cells, bool)
    border[mb.owner[mb.n_internal:]] = True
    X = X[~border[locate(mb, X)]]
    v = Sampler(S, X).sample(["p", "Ux", "Uy", "Uz", "U_mag"])
    assert np.abs(v["p"] - (1 + X @ [1, 2, 3])).max() < 1e-12
    assert set(v) == {"p", "Ux", "Uy", "Uz", "U_mag"}
    assert parse_points("0 1 2 ; 3 4 5", dim=3).shape == (2, 3)
    with pytest.raises(ValueError, match="x y z"):
        parse_points([[0, 1]], dim=3)


def test_restart_3d_and_dimension_mismatch(tmp_path):
    S = build_solver(_box_cfg())
    C = S.mesh.cell_centers
    S.U = np.column_stack([C[:, 0], 2 * C[:, 1], -C[:, 2]])
    save_checkpoint(S, tmp_path / "c.npz")
    S2 = build_solver(_box_cfg())
    assert load_checkpoint(S2, tmp_path / "c.npz")["mode"] == "exact"
    assert np.array_equal(S2.U, S.U)
    S3 = build_solver(_box_cfg(nx=2, ny=3, nz=4))
    assert load_checkpoint(S3, tmp_path / "c.npz")["mode"] == "interpolé"
    cfg2 = {"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0, "y1": 1, "nx": 2,
                     "ny": 2}, "physics": {"nu": 0.1},
            "boundary": {k: {"type": "wall"} for k in SIDES[:4]}}
    with pytest.raises(ValueError, match="reprise d'un calcul 3D sur un maillage 2D"):
        load_checkpoint(build_solver(cfg2), tmp_path / "c.npz")


def test_initial_velocity_formulas_2d_and_3d():
    """[initial] U en formules : documenté mais refusé par le solveur avant ce lot
    (« could not convert string to float »), en 2D comme en 3D."""
    cfg = _box_cfg()
    cfg["initial"] = {"U": [0, "x*z", 1]}
    S = build_solver(cfg)
    C = S.mesh.cell_centers
    assert np.array_equal(S.U[:, 1], C[:, 0] * C[:, 2]) and np.all(S.U[:, 2] == 1.0)
    cfg2 = {"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0, "y1": 1, "nx": 3,
                     "ny": 3}, "physics": {"nu": 0.1}, "initial": {"U": ["sin(y)", 0]},
            "boundary": {k: {"type": "wall"} for k in SIDES[:4]}}
    S = build_solver(cfg2)
    assert np.array_equal(S.U[:, 0], np.sin(S.mesh.cell_centers[:, 1]))


@pytest.mark.parametrize("change,message", [
    ({"boundary": {"left": {"type": "wall", "U": [1, 0]}}},
     r"\[boundary.left\] U = \[1, 0\] : 3 composantes attendues \[x, y, z\]"),
    ({"mesh": {"z1": -1}}, r"\[mesh\] z1 = -1 doit être > z0 = 0"),
    ({"mesh": {"grading": [1, 2]}}, "3 rapports attendus"),
    ({"mesh": {"extrude": {"nz": 2}}}, "sans objet pour un maillage « box »"),
    ({"physics": {"axisymmetric": True}, "porous": [{"region": "rectangle"}],
      "solver": {"algorithm": "coupled"}},
     "axisymétrique, zones poreuses \\[\\[porous\\]\\], solveur couplé"),
    ({"output": {"animate": "vorticity"}}, "animation"),
])
def test_case_checks_in_3d(change, message):
    cfg = _box_cfg()
    for sec, val in change.items():
        cfg[sec] = {**cfg.get(sec, {}), **val} if isinstance(val, dict) else val
    with pytest.raises(ValueError, match=message):
        check_case(cfg)


def test_case_checks_extrude():
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0, "y1": 1, "nx": 4,
                    "ny": 4, "extrude": {"z0": 1, "z1": 0, "nz": 2.0}},
           "physics": {"nu": 0.01, "body_force": [1, 0]},
           "boundary": {"left": {"type": "wall"}}}
    with pytest.raises(ValueError) as e:
        check_case(cfg)
    msg = str(e.value)
    assert "[mesh.extrude] nz = 2.0 : nombre entier attendu" in msg
    assert "[mesh.extrude] z1 = 0 doit être > z0 = 1" in msg
    assert "[physics] body_force = [1, 0] : 3 composantes" in msg
    cfg["mesh"]["extrude"] = {"z1": 2, "nz": 2000}
    cfg["mesh"]["nx"] = 500
    cfg["physics"]["body_force"] = [1, 0, 0]
    assert any("4 000 000 cellules" in w for w in check_case(cfg))
