"""Construction d'un maillage à partir d'une configuration (dict issu d'un TOML/JSON) ou d'un
préréglage. Le type de maillage se choisit avec `mesh.type` :

  blocks        multi-blocs structuré (syntaxe proche de blockMeshDict)
  rectangle     rectangle structuré (avec grading)
  ogrid         structuré en O autour d'un corps
  unstructured  triangles (DistMesh), raffinement autour des corps
  hybrid        couches de quadrilatères aux parois + triangles
  file          lecture d'un fichier .msh (Gmsh) ou .su2
  box           pavé 3D d'hexaèdres (x0…z1, nx, ny, nz, grading, names, periodic)

3D : [mesh.extrude] z0, z1, nz (, grading, names = {back, front}, patch_types, periodic)
extrude n'importe quel maillage 2D ci-dessus selon z (quadrilatères → hexaèdres,
triangles → prismes ; voir mesh3d.extrude).
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

from .blocks import (backward_facing_step_mesh, block_mesh, cavity_mesh, channel_mesh,
                     flat_plate_mesh, rectangle_mesh)
from .geometry import Circle, NACA4, Rectangle, shape_from_dict
from .io import read_mesh
from .mesh import Mesh2D
from .ogrid import o_grid
from .unstructured import hybrid_mesh, triangulate

MESH_TYPES = ("blocks", "rectangle", "ogrid", "unstructured", "hybrid", "file", "box")


def load_config(path) -> dict:
    """Fichier de cas TOML (ou JSON) : BOM et encodage Windows acceptés (remarque affichée
    en avertissement), erreurs de syntaxe en français."""
    from ..tomlio import loads, read_text
    path = Path(path)
    text, note = read_text(path)
    if note:
        warnings.warn(note, stacklevel=2)
    if path.suffix.lower() == ".json":
        try:
            return json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path.name} : JSON incorrect (ligne {exc.lineno}, colonne "
                             f"{exc.colno}) : {exc.msg}.") from None
    try:
        return loads(text)
    except ValueError as exc:
        raise ValueError(f"{path.name} : {exc}") from None


def _outer(cfg) -> Rectangle:
    spec = dict(cfg.get("domain", {"type": "rectangle", "x0": -10, "y0": -10, "x1": 30, "y1": 10,
                                   "names": {"left": "inlet", "right": "outlet",
                                             "bottom": "bottom", "top": "top"}}))
    spec.setdefault("type", "rectangle")      # [domain] sans type (interface) : rectangle
    return shape_from_dict(spec)


def _bodies(cfg, base_dir="."):
    out = []
    for spec in cfg.get("bodies", []):
        spec = dict(spec)
        spec.setdefault("patch_type", "wall")
        out.append(shape_from_dict(spec, base_dir))
    return out


def build_mesh(cfg: dict, base_dir=".", verbose: bool = False) -> Mesh2D:
    """Maillage décrit par la section [mesh] ; [mesh] cut_axis = true garde la moitié y > 0
    (calcul axisymétrique d'un corps de révolution, patch « axis ») ; [mesh.extrude] :
    maillage 3D extrudé selon z."""
    mesh = _build_mesh(cfg, base_dir, verbose)
    if cfg.get("mesh", cfg).get("cut_axis"):
        mesh = mesh.cut_at_axis()
    return extrude_from(mesh, cfg.get("mesh", cfg))


def extrude_from(mesh, m: dict):
    """Applique [mesh.extrude] (si présente) à un maillage 2D."""
    e = m.get("extrude")
    if not e:
        return mesh
    if getattr(mesh, "dim", 2) == 3:
        raise ValueError("[mesh.extrude] : le maillage est déjà en 3D (type « box »).")
    from ..mesh3d import extrude
    per = e.get("periodic")
    return extrude(mesh, e.get("z0", 0.0), e.get("z1", 1.0), e.get("nz", 1),
                   e.get("grading", 1.0), e.get("names"), e.get("patch_types"),
                   [tuple(p) for p in per] if per else None)


def _build_mesh(cfg: dict, base_dir=".", verbose: bool = False) -> Mesh2D:
    """Construit le maillage décrit par la section [mesh] (+ [domain], [[bodies]])."""
    m = cfg.get("mesh", cfg)
    kind = m.get("type", "unstructured").lower()
    if kind == "file":
        if not str(m.get("path") or "").strip():
            raise ValueError("[mesh] path manquant : fichier .msh (Gmsh) ou .su2 à importer.")
        p = Path(m["path"])
        if not p.is_absolute():
            p = Path(base_dir) / p
        if not p.is_file():
            raise ValueError(f"[mesh] path : fichier de maillage introuvable : {p}")
        return read_mesh(p, m.get("patch_types"))
    if kind == "blocks":
        return block_mesh(m["vertices"], m["blocks"], m.get("edges"), m.get("patches"),
                          m.get("periodic"))
    if kind == "box":
        from ..mesh3d import box_mesh
        return box_mesh(m["x0"], m["x1"], m["y0"], m["y1"], m["z0"], m["z1"], m["nx"],
                        m["ny"], m["nz"], tuple(m.get("grading", (1.0, 1.0, 1.0))),
                        m.get("names"), m.get("patch_types"),
                        [tuple(p) for p in m["periodic"]] if m.get("periodic") else None)
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
