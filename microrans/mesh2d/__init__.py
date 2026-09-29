"""Mailleur 2D : géométrie CSG, maillages structurés (blocs, O-grid), non structurés
(triangles DistMesh), hybrides (couches limites + triangles), import/export."""
from .blocks import (backward_facing_step_mesh, block_mesh, cavity_mesh, channel_mesh,
                     flat_plate_mesh, grading_distribution, rectangle_mesh)
from .geometry import (NACA4, Circle, Ellipse, Polygon, Rectangle, Shape, Spline,
                       shape_from_dict)
from .io import read_curve, read_mesh, write_mesh
from .mesh import Mesh2D, Patch
from .ogrid import o_grid
from .unstructured import hybrid_mesh, size_function, triangle_quality, triangulate

__all__ = ["Mesh2D", "Patch", "Shape", "Circle", "Rectangle", "Polygon", "Ellipse", "NACA4",
           "Spline", "shape_from_dict", "block_mesh", "rectangle_mesh", "channel_mesh",
           "cavity_mesh", "backward_facing_step_mesh", "flat_plate_mesh", "grading_distribution",
           "o_grid", "triangulate", "hybrid_mesh", "size_function", "triangle_quality",
           "read_mesh", "write_mesh", "read_curve"]
