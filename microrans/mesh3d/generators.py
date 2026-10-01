"""Maillages 3D hexaédriques : boîte structurée et extrusion d'un maillage 2D.

L'extrusion reprend n'importe quel maillage 2D du mailleur (rectangle, blocs, O, triangles,
hybride) : quadrilatères → hexaèdres, triangles → prismes ; les frontières 2D deviennent des
frontières 3D de même nom et même type, avec deux frontières de plus aux extrémités
(« back » en z0, « front » en z1 par défaut).
"""
from __future__ import annotations

import numpy as np

from ..mesh2d.blocks import grading_distribution
from .mesh import Mesh3D

_SIDES = ("left", "right", "bottom", "top", "back", "front")


def box_mesh(x0, x1, y0, y1, z0, z1, nx, ny, nz, grading=(1.0, 1.0, 1.0), names=None,
             types=None, periodic=None) -> Mesh3D:
    """Pavé [x0, x1] × [y0, y1] × [z0, z1] de nx × ny × nz hexaèdres.

    grading : rapport dernière / première maille par direction (ou multi-grading, comme
    rectangle_mesh) ; names : {'left', 'right', 'bottom', 'top', 'back', 'front'} → nom de
    frontière (x0, x1, y0, y1, z0, z1) ; types : {nom: wall | patch | symmetry} ;
    periodic : [(nom_A, nom_B), …]."""
    for n, lab in ((nx, "nx"), (ny, "ny"), (nz, "nz")):
        if int(n) != n or n < 1:
            raise ValueError(f"{lab} = {n} : entier ≥ 1 attendu.")
    nx, ny, nz = int(nx), int(ny), int(nz)
    if not (x1 > x0 and y1 > y0 and z1 > z0):
        raise ValueError("Boîte vide : x1 > x0, y1 > y0 et z1 > z0 attendus.")
    g = list(grading) if not np.isscalar(grading) else [grading] * 3
    xs = x0 + (x1 - x0) * grading_distribution(nx, g[0])
    ys = y0 + (y1 - y0) * grading_distribution(ny, g[1])
    zs = z0 + (z1 - z0) * grading_distribution(nz, g[2])
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    pts = np.column_stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")])

    def vid(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    i, j, k = i.ravel(order="F"), j.ravel(order="F"), k.ravel(order="F")
    cells = np.column_stack([vid(i, j, k), vid(i + 1, j, k), vid(i + 1, j + 1, k),
                             vid(i, j + 1, k), vid(i, j, k + 1), vid(i + 1, j, k + 1),
                             vid(i + 1, j + 1, k + 1), vid(i, j + 1, k + 1)])

    def quads(a, b, fixed, axis):
        """Faces d'une face du pavé : grille (a, b) des deux autres indices."""
        A, B = np.meshgrid(np.arange(a), np.arange(b), indexing="ij")
        A, B = A.ravel(), B.ravel()
        if axis == 0:
            return np.column_stack([vid(fixed, A, B), vid(fixed, A + 1, B),
                                    vid(fixed, A + 1, B + 1), vid(fixed, A, B + 1)])
        if axis == 1:
            return np.column_stack([vid(A, fixed, B), vid(A + 1, fixed, B),
                                    vid(A + 1, fixed, B + 1), vid(A, fixed, B + 1)])
        return np.column_stack([vid(A, B, fixed), vid(A + 1, B, fixed),
                                vid(A + 1, B + 1, fixed), vid(A, B + 1, fixed)])

    faces = {"left": quads(ny, nz, 0, 0), "right": quads(ny, nz, nx, 0),
             "bottom": quads(nx, nz, 0, 1), "top": quads(nx, nz, ny, 1),
             "back": quads(nx, ny, 0, 2), "front": quads(nx, ny, nz, 2)}
    return _assemble(pts, cells, faces, names, types, periodic)


def extrude(mesh2d, z0: float = 0.0, z1: float = 1.0, nz: int = 1, grading=1.0,
            names=None, types=None, periodic=None) -> Mesh3D:
    """Extrusion d'un maillage 2D selon z, de z0 à z1, en nz couches.

    names : {'back': …, 'front': …} (frontières en z0 et z1) ; types : types des frontières
    (défaut : celui de la frontière 2D, « patch » pour back / front) ; periodic : paires
    (les paires 2D sont reprises)."""
    if int(nz) != nz or nz < 1:
        raise ValueError(f"nz = {nz} : entier ≥ 1 attendu.")
    if not z1 > z0:
        raise ValueError("Extrusion : z1 > z0 attendu.")
    nz = int(nz)
    p2 = mesh2d.points
    npt = len(p2)
    zs = z0 + (z1 - z0) * grading_distribution(nz, grading)
    pts = np.column_stack([np.tile(p2, (nz + 1, 1)), np.repeat(zs, npt)])
    cells = []
    for c in mesh2d.cells_as_lists():
        if len(c) not in (3, 4):
            raise ValueError(f"Extrusion : cellule 2D à {len(c)} sommets (triangles et "
                             "quadrilatères seulement).")
        cells.append(c)
    cn = [c for c in cells]
    out = []
    for k in range(nz):
        lo, hi = k * npt, (k + 1) * npt
        out.extend(np.concatenate([np.asarray(c) + lo, np.asarray(c) + hi]) for c in cn)
    faces = {}
    for name, ptype, edges in mesh2d.all_boundary_patches():
        e = np.asarray(edges)
        q = [np.column_stack([e[:, 0] + k * npt, e[:, 1] + k * npt,
                              e[:, 1] + (k + 1) * npt, e[:, 0] + (k + 1) * npt])
             for k in range(nz)]
        faces[name] = np.vstack(q)
    base = {"back": "back", "front": "front", **(names or {})}

    def caps(off):
        tri = [c for c in cn if len(c) == 3]
        qua = [c for c in cn if len(c) == 4]
        rows = [np.concatenate([np.asarray(c), [-1]]) + np.where(np.arange(4) < 3, off, 0)
                for c in tri] + [np.asarray(c) + off for c in qua]
        return np.array(rows, dtype=np.int64)

    for side, off in (("back", 0), ("front", nz * npt)):
        nm = base[side]
        f = caps(off)
        faces[nm] = np.vstack([faces[nm], f]) if nm in faces else f
    t2 = {name: ptype for name, ptype, _ in mesh2d.all_boundary_patches() if ptype != "cyclic"}
    t2.update(types or {})
    per = [tuple(pp[:2]) for pp in mesh2d.periodic_pairs] + list(periodic or [])
    return _assemble(pts, out, faces, None, t2, per)


def _assemble(pts, cells, faces, names, types, periodic) -> Mesh3D:
    """Regroupe les faces par nom de frontière (plusieurs côtés peuvent porter le même nom)."""
    names = {**{s: s for s in _SIDES}, **(names or {})}
    boundary = {}
    for side, f in faces.items():
        nm = names.get(side, side)
        boundary[nm] = np.vstack([boundary[nm], f]) if nm in boundary else f
    return Mesh3D(pts, cells, boundary, patch_types=types or {}, periodic=periodic)
