"""Maillage structuré en C autour d'un profil (`[mesh] type = "cgrid"`)."""
import numpy as np
import pytest

from microrans.mesh2d.builder import build_mesh


def _cfg(**mesh):
    return {"mesh": {"type": "cgrid", "n_around": 64, "n_radial": 24, "n_wake": 12,
                     "farfield_radius": 15.0, "first_height": 1e-4, **mesh},
            "bodies": [{"type": "naca", "code": "0012", "name": "airfoil"}]}


def test_cgrid_topology_and_wake_cut():
    m = build_mesh(_cfg())
    assert m.n_cells == (64 + 2 * 12) * 24
    q = m.cell_nodes[:, :4]
    assert np.all(m.cell_nv == 4) and len(np.unique(q)) == len(m.points)
    names = {p.name: p.type for p in m.patches}
    assert names == {"airfoil": "wall", "farfield": "patch"}
    sizes = {p.name: p.size for p in m.patches}
    assert sizes["airfoil"] == 64
    assert sizes["farfield"] == (64 + 2 * 12) + 2 * 24      # bord en C + deux sorties
    # la coupure de sillage (y = 0, x > 1) est faite de faces internes, pas de frontières
    fc = m.face_centers
    cut = (np.abs(fc[:, 1]) < 1e-12) & (fc[:, 0] > 1.0)
    assert cut.sum() == 12 and np.all(np.flatnonzero(cut) < m.n_internal)
    # 1re maille : first_height à la paroi (centre à mi-hauteur)
    from microrans.mesh2d.cgrid import c_grid  # noqa: F401  (module du maillage en C)
    assert m.wall_distance[m.owner[m.patch("airfoil").faces]].min() == pytest.approx(
        0.5e-4, rel=0.05)


def test_cgrid_cells_valid_with_incidence_and_camber():
    for body in ({"type": "naca", "code": "0012", "incidence": 10.0,
                  "rotation_center": [0.25, 0.0], "name": "airfoil"},
                 {"type": "naca", "code": "4412", "trailing_edge": "sharp", "name": "airfoil"}):
        cfg = _cfg()
        cfg["bodies"] = [body]
        m = build_mesh(cfg)
        assert np.all(m.cell_volumes > 0)
        assert m.quality()["non_orthogonality_max_deg"] < 80


def test_cgrid_refuses_blunt_trailing_edge_and_two_bodies():
    cfg = _cfg()
    cfg["bodies"][0]["trailing_edge"] = "open"
    with pytest.raises(ValueError, match="bord de fuite"):
        build_mesh(cfg)
    cfg = _cfg()
    cfg["bodies"] = [{"type": "circle", "center": [0, 0], "radius": 0.5, "name": "c"}]
    with pytest.raises(ValueError, match="bord de fuite"):
        build_mesh(cfg)
    cfg = _cfg()
    cfg["bodies"].append({"type": "circle", "center": [3, 0], "radius": 0.5, "name": "c"})
    with pytest.raises(ValueError, match="exactement un"):
        build_mesh(cfg)


def test_cgrid_keys_validated_and_coarsened():
    from microrans.fv2d.fmg import coarsen_config
    from microrans.fv2d.validate import check_case
    cfg = _cfg(wake_length=20.0)
    cfg.update({"physics": {"reynolds": 1e5, "model": "laminar"},
                "boundary": {"airfoil": {"type": "wall"},
                             "farfield": {"type": "farfield", "U": [1.0, 0.0]}}})
    assert check_case(cfg) == []
    bad = {**cfg, "mesh": {**cfg["mesh"], "n_wake": 0}}
    with pytest.raises(ValueError, match="n_wake"):
        check_case(bad)
    c = coarsen_config(cfg, 1)["mesh"]
    assert (c["n_around"], c["n_radial"], c["n_wake"], c["first_height"]) == (32, 12, 6, 2e-4)


def test_cgrid_troubleshooting_quotes_real_messages():
    """docs/depannage.md cite les vrais messages du maillage en C."""
    from pathlib import Path

    from microrans.fv2d.validate import check_case
    doc = (Path(__file__).resolve().parents[1] / "docs" / "depannage.md").read_text(
        encoding="utf-8")
    cfg = _cfg()
    cfg["bodies"][0]["trailing_edge"] = "open"
    with pytest.raises(ValueError) as blunt:
        build_mesh(cfg)
    two = _cfg()
    two["bodies"].append({"type": "circle", "center": [3, 0], "radius": 0.5, "name": "c"})
    two.update({"physics": {"reynolds": 1e5}, "boundary": {"airfoil": {"type": "wall"}}})
    with pytest.raises(ValueError) as many:
        check_case(two)
    for quote, msg in (("Maillage en C : pas de bord de fuite pointu", str(blunt.value)),
                       ("le maillage en C entoure exactement un profil", str(many.value))):
        assert f"| `{quote}`" in doc and quote in msg
