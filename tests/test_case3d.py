"""Cas 3D décrits par le fichier de cas : [mesh] type = "box" et [mesh.extrude], sortie VTK
des hexaèdres et prismes, sondes et profils en 3D, reprise, vérifications avant calcul."""
import numpy as np
import pytest

from microrans.fv2d.case import build_solver, case_mesh
from microrans.fv2d.restart import load_checkpoint, save_checkpoint
from microrans.fv2d.sampling import Sampler, locate, parse_points
from microrans.fv2d.validate import check_case
from microrans.mesh2d import Circle, Rectangle
from microrans.mesh2d.io import read_vtk, write_mesh, write_vtk
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


def test_vtk_output_of_hexahedra_and_prisms(tmp_path):
    """Numérotation VTK des sommets (volumes recalculés par la bibliothèque VTK 9.7 égaux
    aux nôtres à 1e-8 avec le fichier texte : vérifié hors des tests ; lecture par VTK :
    tests/test_vtk.py, ignorée si VTK n'est pas installé)."""
    m2 = triangulate(Rectangle(0, 0, 2, 1) - Circle((1, .5), .2).as_wall(), 0.2, [])
    for m in (box_mesh(0, 1, 0, 1, 0, 1, 2, 3, 4), extrude(m2, 0, 1, 2)):
        U = m.cell_centers.copy()
        for binary in (True, False):
            write_vtk(m, tmp_path / "f.vtk", {"U": U, "p": m.cell_volumes}, binary=binary)
            r = read_vtk(tmp_path / "f.vtk")
            flat = np.concatenate([[len(c), *c] for c in m.cells_as_lists()])
            assert np.array_equal(r["cells"], flat)
            assert set(r["cell_types"]) == ({12} if m.n_cells == 24 else {13})
            assert r["points"].shape == (m.n_points, 3) and r["cell_data"]["U"].shape[1] == 3
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
    w = [w for w in check_case(cfg) if "4 000 000 cellules (3D)" in w]
    assert w and "~11.9 Go de mémoire, ~34 à 56 s par itération" in w[0]
    assert "~4 à 8 min" in w[0]


def _cyl_cfg(dim, alpha=0.0, nz=1, **extra):
    """Cylindre Re = 20, maillage en O grossier ; 3D : extrudé, plans de symétrie en z.
    Une couche : solution identique au 2D. Plusieurs couches : la diffusion à travers les
    faces z intérieures entre dans a_P (interpolation de Rhie-Chow) — solution différente
    à l'ordre de la discrétisation (Cd 2.17561 au lieu de 2.17575 avec 2 couches) et
    convergence SIMPLE 2 fois plus lente (276 itérations au lieu de 141)."""
    z = [0.0] * (dim - 2)
    cfg = {"mesh": {"type": "ogrid", "n_around": 32, "n_radial": 16, "farfield_radius": 15.0,
                    "first_height": 0.02},
           "bodies": [{"type": "circle", "name": "cylinder", "radius": 0.5}],
           "physics": {"reynolds": 20, "angle_of_attack": alpha},
           "boundary": {"cylinder": {"type": "wall"},
                        "farfield": {"type": "farfield", "U": [1.0, 0.0] + z}},
           "initial": {"U": [1.0, 0.0] + z},
           "solver": {"max_iter": 400, "tol": 1e-8},
           "output": {"probes": [[2.0, 0.3] + [0.25] * (dim - 2)], "moment_center": [0.25, 0.0],
                      "lines": [{"name": "wake", "start": [0.6, 0.0] + [0.25] * (dim - 2),
                                 "end": [4.0, 0.0] + [0.25] * (dim - 2), "n": 10}]}}
    if dim == 3:
        cfg["mesh"]["extrude"] = {"z1": 0.5, "nz": nz,
                                  "patch_types": {"back": "symmetry", "front": "symmetry"}}
        cfg["boundary"].update(back={"type": "symmetry"}, front={"type": "symmetry"})
    for k, v in extra.items():
        cfg[k] = {**cfg.get(k, {}), **v}
    return cfg


def test_run_case_3d_matches_2d_with_incidence(tmp_path):
    """Chaîne complète en 3D (efforts en axes écoulement, moment, sondes, profils, CSV de
    paroi, figures du plan médian, VTK) : coefficients égaux au 2D (A_ref = L × envergure)."""
    from microrans.fv2d.case import run_case
    s2 = run_case(_cyl_cfg(2, 10.0), out_dir=tmp_path / "2d", verbose=False, plot=False)
    s3 = run_case(_cyl_cfg(3, 10.0), out_dir=tmp_path / "3d", verbose=False)
    assert s3["dimension"] == 3 and s3["reference_area"] == pytest.approx(0.5)
    a, b = s2["cylinder"], s3["cylinder"]
    for k in ("Cd", "Cl", "Cm", "Cd_pressure"):
        assert b[k] == pytest.approx(a[k], rel=1e-9, abs=1e-12)
    # Cl ≈ 6e-4 (et non 0) : maillage en O de 32 mailles incliné de 10° sur l'écoulement
    assert abs(b["Cs"]) < 1e-12 and abs(a["Cl"]) < 1e-3
    p2, p3 = s2["probes"][0], s3["probes"][0]
    assert p3["z"] == 0.25 and abs(p3["Uz"]) < 1e-12
    assert p3["Ux"] == pytest.approx(p2["Ux"], rel=1e-9)
    out = tmp_path / "3d"
    for f in ("fields.vtk", "U.png", "p.png", "vorticity.png", "convergence.png",
              "line_wake.csv", "wall_cylinder.csv", "checkpoint.npz"):
        assert (out / f).is_file(), f
    assert (out / "wall_cylinder.csv").read_text().splitlines()[0] == \
        "x,y,z,tau_w,Cf,Cp,yplus"
    assert (out / "line_wake.csv").read_text().splitlines()[0].startswith("s,x,y,z,Ux,Uy,Uz")
    from microrans.fv2d.report import summary_text
    txt = summary_text(s3)
    assert ", 3D —" in txt and " Cs " in txt and "A_ref = 0.5" in txt


def test_run_case_3d_transient_vtk_and_averages(tmp_path):
    from microrans.fv2d.case import run_case
    cfg = _cyl_cfg(3, nz=2, solver={"mode": "transient", "dt": 0.2, "t_end": 1.0},
                   output={"vtk_every": 2, "average_from": 0.4})
    s = run_case(cfg, out_dir=tmp_path, verbose=False, plot=False)
    assert s["steps"] == 5 and sorted(p.name for p in tmp_path.glob("fields_*.vtk")) == [
        "fields_000002.vtk", "fields_000004.vtk"]
    names = read_vtk(tmp_path / "fields.vtk")["cell_data"]
    assert names["U"].shape[1] == 3
    head = (tmp_path / "history.csv").read_text().splitlines()[0]
    assert "Cs_cylinder" in head and "probe1_Uz" in head
    assert names["Uz_mean"].ndim == 1


def test_cli_mesh_3d_writes_vtk_only(tmp_path, capsys):
    from microrans.cli import main
    case = tmp_path / "box.toml"
    case.write_text('[mesh]\ntype = "box"\nx0 = 0\nx1 = 1\ny0 = 0\ny1 = 1\nz0 = 0\nz1 = 1\n'
                    "nx = 2\nny = 2\nnz = 2\n", encoding="utf-8")
    assert main(["mesh", str(case), "-o", str(tmp_path / "m"), "-q"]) == 0
    assert sorted(p.name for p in (tmp_path / "m").iterdir()) == ["box.vtk", "quality.json"]
    assert main(["mesh", str(case), "-o", str(tmp_path / "m2"), "-f", "su2", "vtk"]) == 0
    assert "format(s) su2 non disponible(s)" in capsys.readouterr().out


def test_example_square_duct_via_cli(tmp_path, capsys):
    """Exemple 3D fourni, maillage réduit (16²) : débit à +1.5 % de la série exacte
    (32² dans l'exemple : +0.38 %)."""
    import json

    from microrans.cli import main
    assert main(["run2d", "conduite_carree_3d", "-o", str(tmp_path), "--no-plot",
                 "--set", "mesh.ny=16", "mesh.nz=16"]) == 0
    out = capsys.readouterr().out
    assert "Cas 3D : 512 cellules" in out and "Calcul stationnaire, laminaire, 3D" in out
    s = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    assert s["converged"] and s["U_mean"][0] / 3.5144253739 - 1 == pytest.approx(0.015, abs=2e-3)
    assert (tmp_path / "line_diagonale.csv").is_file()
