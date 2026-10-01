"""Maillages 3D volumes finis : hexaèdres (et prismes, tétraèdres, pyramides), pavé
structuré et extrusion d'un maillage 2D."""
from .generators import box_mesh, extrude
from .mesh import Mesh3D

__all__ = ["Mesh3D", "box_mesh", "extrude"]
