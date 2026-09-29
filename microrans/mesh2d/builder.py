"""Construction d'un maillage à partir d'une configuration (dict issu d'un TOML/JSON) ou d'un
préréglage. Le type de maillage se choisit avec `mesh.type` :

  blocks        multi-blocs structuré (syntaxe proche de blockMeshDict)
  rectangle     rectangle structuré (avec grading)
  ogrid         structuré en O autour d'un corps
  unstructured  triangles (DistMesh), raffinement autour des corps
  hybrid        couches de quadrilatères aux parois + triangles
  file          lecture d'un fichier .msh (Gmsh) ou .su2
"""
from __future__ import annotations

import json
from pathlib import Path

from .blocks import (backward_facing_step_mesh, block_mesh, cavity_mesh, channel_mesh,
                     flat_plate_mesh, rectangle_mesh)
from .geometry import Circle, NACA4, Rectangle, shape_from_dict
from .io import read_mesh
from .mesh import Mesh2D
from .ogrid import o_grid
from .unstructured import hybrid_mesh, triangulate

MESH_TYPES = ("blocks", "rectangle", "ogrid", "unstructured", "hybrid", "file")


def load_config(path) -> dict:
    path = Path(path)
    if path.suffix.lower() == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    try:
        import tomllib
    except ModuleNotFoundError:          # Python 3.10
        import tomli as tomllib
    return tomllib.loads(path.read_text(encoding="utf-8"))


def _outer(cfg) -> Rectangle:
    spec = dict(cfg.get("domain", {"type": "rectangle", "x0": -10, "y0": -10, "x1": 30, "y1": 10,
                                   "names": {"left": "inlet", "right": "outlet",
                                             "bottom": "bottom", "top": "top"}}))
    return shape_from_dict(spec)


def _bodies(cfg, base_dir="."):
    out = []
    for spec in cfg.get("bodies", []):
        spec = dict(spec)
        spec.setdefault("patch_type", "wall")
        out.append(shape_from_dict(spec, base_dir))
    return out


def build_mesh(cfg: dict, base_dir=".", verbose: bool = False) -> Mesh2D:
    """Construit le maillage décrit par la section [mesh] (+ [domain], [[bodies]])."""
    m = cfg.get("mesh", cfg)
    kind = m.get("type", "unstructured").lower()
    if kind == "file":
        p = Path(m["path"])
        if not p.is_absolute():
            p = Path(base_dir) / p
        return read_mesh(p, m.get("patch_types"))
    if kind == "blocks":
        return block_mesh(m["vertices"], m["blocks"], m.get("edges"), m.get("patches"),
                          m.get("periodic"))
    if kind == "rectangle":
        return rectangle_mesh(m["x0"], m["x1"], m["y0"], m["y1"], m["nx"], m["ny"],
                              tuple(m.get("grading", (1.0, 1.0))), m.get("names"),
                              m.get("patch_types"),
                              [tuple(p) for p in m["periodic"]] if m.get("periodic") else None)
    bodies = _bodies(cfg, base_dir)
    if kind == "ogrid":
        if len(bodies) != 1:
            raise ValueError("Le maillage en O demande exactement un corps.")
        return o_grid(bodies[0], m.get("n_around", 128), m.get("n_radial", 64),
                      m.get("farfield_radius", 20.0), m.get("first_height", 1e-3),
                      m.get("center"))
    outer = _outer(cfg)
    refs = []
    for b in bodies:
        refs.append({"shape": b, "h": m.get("h_surface", 0.05), "growth": m.get("growth", 0.2)})
    for r in m.get("refinements", []):
        r = dict(r)
        shape = shape_from_dict(r.pop("shape"))
        refs.append({"shape": shape, **r})
    if kind == "unstructured":
        domain = outer
        for b in bodies:
            domain = domain - b
        return triangulate(domain, m.get("h_max", 1.0), refs, max_iter=m.get("max_iter", 300),
                           verbose=verbose)
    if kind == "hybrid":
        lay = m.get("layers", {})
        return hybrid_mesh(outer, bodies, m.get("h_max", 1.0), m.get("h_surface", 0.05),
                           lay.get("n", 10), lay.get("first_height", 1e-3), lay.get("ratio", 1.2),
                           m.get("growth", 0.2), refinements=[r for r in refs if r["shape"]
                                                              not in bodies],
                           max_iter=m.get("max_iter", 300), verbose=verbose)
    raise ValueError(f"Type de maillage inconnu '{kind}'. Choix : {', '.join(MESH_TYPES)}")


def _cyl_domain():
    return Rectangle(-10, -10, 30, 10, names={"left": "inlet", "right": "outlet",
                                              "bottom": "bottom", "top": "top"})


PRESETS = {
    "cavity": lambda: cavity_mesh(64, grading=2.0),
    "channel": lambda: channel_mesh(1.0, 2.0, 4, 96, first_height=1e-3),
    "backstep": lambda: backward_facing_step_mesh(),
    "flatplate": lambda: flat_plate_mesh(),
    "cylinder-ogrid": lambda: o_grid(Circle((0, 0), 0.5, "cylinder"), 128, 96, 30.0, 5e-3),
    "cylinder-tri": lambda: triangulate(
        _cyl_domain() - Circle((0, 0), 0.5, "cylinder").as_wall(), 1.0,
        [{"shape": Circle((0, 0), 0.5), "h": 0.04, "growth": 0.12}]),
    "cylinder-hybrid": lambda: hybrid_mesh(
        _cyl_domain(), [Circle((0, 0), 0.5, "cylinder").as_wall()], 1.0, 0.04, 12, 2e-3, 1.2, 0.12),
    "naca0012-ogrid": lambda: o_grid(NACA4("0012", name="airfoil"), 192, 96, 30.0, 1e-5),
    "naca0012-hybrid": lambda: hybrid_mesh(
        Rectangle(-5, -5, 10, 5, names={"left": "inlet", "right": "outlet", "bottom": "bottom",
                                        "top": "top"}),
        [NACA4("0012", name="airfoil").as_wall()], 0.6, 0.01, 20, 1e-5, 1.3, 0.15),
}
