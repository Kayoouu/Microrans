"""Mailleur structuré multi-blocs, inspiré de blockMesh (OpenFOAM).

Chaque bloc est un quadrilatère (4 sommets, sens trigonométrique : v0 bas-gauche, v1 bas-droit,
v2 haut-droit, v3 haut-gauche) découpé en nx × ny mailles, avec progression géométrique
(« simpleGrading ») ou multi-progression, et arêtes éventuellement courbes (arc, polyligne,
spline). L'intérieur est obtenu par interpolation transfinie. Les blocs sont recollés par
fusion des sommets coïncidents (les nombres de mailles doivent correspondre).
"""
from __future__ import annotations

import numpy as np

from .mesh import Mesh2D, merge_points


def grading_distribution(n: int, grading=1.0) -> np.ndarray:
    """Abscisses normalisées (n+1 valeurs de 0 à 1).

    grading : rapport taille dernière / première maille (comme blockMesh), ou liste de
    segments [(fraction_longueur, fraction_mailles, rapport), ...] (multi-grading).
    """
    if np.isscalar(grading):
        segs = [(1.0, 1.0, float(grading))]
    else:
        segs = [tuple(map(float, s)) for s in grading]
    if any(len(sg) != 3 or sg[2] <= 0 or sg[0] <= 0 or sg[1] <= 0 for sg in segs):
        raise ValueError(f"Progression (grading) invalide : {grading} — rapports > 0 attendus, "
                         "ou segments (fraction de longueur, fraction de mailles, rapport).")
    lf = np.array([s[0] for s in segs])
    nf = np.array([s[1] for s in segs])
    lf, nf = lf / lf.sum(), nf / nf.sum()
    counts = np.maximum(np.round(nf * n).astype(int), 1)
    counts[-1] = n - counts[:-1].sum()
    if counts[-1] < 1:
        raise ValueError("Trop de segments de progression pour le nombre de mailles.")
    s = [0.0]
    pos = 0.0
    for L, m, (_, _, r) in zip(lf, counts, segs):
        if m > 1 and abs(r - 1.0) > 1e-12:
            q = r ** (1.0 / (m - 1))
            d = q ** np.arange(m)
        else:
            d = np.ones(m)
        d = d / d.sum() * L
        s.extend(pos + np.cumsum(d))
        pos += L
    s = np.array(s)
    s[-1] = 1.0
    return s


def grading_for_first_cell(n: int, first_fraction: float) -> float:
    """Rapport d'expansion (dernière/première) donnant une 1re maille = first_fraction × longueur."""
    if first_fraction * n >= 1.0:
        return 1.0
    lo, hi = 1.0 + 1e-12, 1e6
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if grading_distribution(n, mid)[1] > first_fraction:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


class _Edge:
    """Courbe paramétrée par l'abscisse curviligne normalisée s ∈ [0, 1]."""

    def __init__(self, p0, p1, spec=None, n_sample=400):
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        kind = (spec or {}).get("type", "line").lower()
        if kind == "line":
            pts = np.vstack([p0, p1])
        elif kind == "arc":
            pts = _arc_points(p0, np.asarray(spec["point"], float), p1, n_sample)
        elif kind in ("polyline", "polyLine"):
            pts = np.vstack([p0, np.asarray(spec["points"], float), p1])
        elif kind == "spline":
            from scipy.interpolate import CubicSpline
            ctrl = np.vstack([p0, np.asarray(spec["points"], float), p1])
            t = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(ctrl, axis=0), axis=1))])
            pts = CubicSpline(t, ctrl)(np.linspace(0.0, t[-1], n_sample))
        else:
            raise ValueError(f"Type d'arête inconnu : {kind}")
        s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        self.s = s / s[-1]
        self.pts = pts

    def __call__(self, u):
        u = np.asarray(u, float)
        return np.column_stack([np.interp(u, self.s, self.pts[:, 0]),
                                np.interp(u, self.s, self.pts[:, 1])])


def _arc_points(a, m, b, n):
    """Arc de cercle passant par a, m, b."""
    ax, ay = a
    bx, by = m
    cx, cy = b
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-300:
        return np.vstack([a, b])
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    c = np.array([ux, uy])
    r = np.linalg.norm(a - c)
    t0 = np.arctan2(*(a - c)[::-1])
    tm = np.arctan2(*(m - c)[::-1])
    t1 = np.arctan2(*(b - c)[::-1])
    # choisir le sens qui passe par m
    def unwrap(t, ref):
        return ref + np.mod(t - ref, 2 * np.pi)
    tm_ccw, t1_ccw = unwrap(tm, t0), unwrap(t1, t0)
    if tm_ccw <= t1_ccw:
        t = np.linspace(t0, t1_ccw, n)
    else:
        t = np.linspace(t0, t1_ccw - 2 * np.pi, n)
    return c + r * np.column_stack([np.cos(t), np.sin(t)])


SIDES = {"bottom": (0, 1), "right": (1, 2), "top": (2, 3), "left": (3, 0)}


def block_mesh(vertices, blocks, edges=None, patches=None, periodic=None,
               default_patch: str = "defaultFaces") -> Mesh2D:
    """Maillage multi-blocs.

    vertices : [[x, y], ...]
    blocks   : [{"vertices": [i0, i1, i2, i3], "cells": [nx, ny], "grading": [gx, gy]}, ...]
               gx/gy : rapport d'expansion ou liste de segments (multi-grading).
    edges    : [{"type": "arc", "vertices": [i, j], "point": [x, y]},
                {"type": "polyline" | "spline", "vertices": [i, j], "points": [[x, y], ...]}]
    patches  : {nom: {"type": "wall" | "patch" | "symmetry" | "empty",
                      "faces": [[i, j], ...]}}   (faces = paires de sommets de blocs)
    periodic : [(patch_A, patch_B), ...]
    """
    V = np.asarray(vertices, float)
    edge_specs = {}
    for e in edges or []:
        i, j = e["vertices"]
        edge_specs[(i, j)] = (e, False)
        edge_specs[(j, i)] = (e, True)

    def make_edge(i, j):
        if (i, j) in edge_specs:
            spec, rev = edge_specs[(i, j)]
            if rev:
                ed = _Edge(V[j], V[i], spec)
                return lambda u: ed(1.0 - np.asarray(u))
            return _Edge(V[i], V[j], spec)
        return _Edge(V[i], V[j])

    all_pts, all_cells, sides_pts = [], [], {}
    offset = 0
    for b, blk in enumerate(blocks):
        v0, v1, v2, v3 = blk["vertices"]
        nx, ny = blk["cells"]
        gx, gy = blk.get("grading", [1.0, 1.0])
        sx, sy = grading_distribution(nx, gx), grading_distribution(ny, gy)
        B, T = make_edge(v0, v1)(sx), make_edge(v3, v2)(sx)
        L, R = make_edge(v0, v3)(sy), make_edge(v1, v2)(sy)
        xi, eta = sx[:, None, None], sy[None, :, None]
        P = ((1 - eta) * B[:, None] + eta * T[:, None] + (1 - xi) * L[None] + xi * R[None]
             - ((1 - xi) * (1 - eta) * V[v0] + xi * (1 - eta) * V[v1]
                + xi * eta * V[v2] + (1 - xi) * eta * V[v3]))
        idx = offset + np.arange((nx + 1) * (ny + 1)).reshape(nx + 1, ny + 1)
        all_pts.append(P.reshape(-1, 2))
        q = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], axis=-1)
        all_cells.append(q.reshape(-1, 4))
        sides_pts[b] = {"bottom": idx[:, 0], "right": idx[-1, :], "top": idx[::-1, -1],
                        "left": idx[0, ::-1]}
        offset += idx.size
    P, new = merge_points(np.vstack(all_pts))
    cells = new[np.vstack(all_cells)]

    boundary, types = {}, {}
    for name, spec in (patches or {}).items():
        types[name] = spec.get("type", "patch")
        edges_list = []
        for i, j in spec["faces"]:
            found = False
            for b, blk in enumerate(blocks):
                bv = list(blk["vertices"])
                for side, (a, c) in SIDES.items():
                    if {bv[a], bv[c]} == {i, j}:
                        chain = new[sides_pts[b][side]]
                        edges_list.append(np.column_stack([chain[:-1], chain[1:]]))
                        found = True
            if not found:
                raise ValueError(f"Patch '{name}' : la face ({i}, {j}) n'est pas un côté de bloc.")
        boundary[name] = np.vstack(edges_list)
    return Mesh2D(P, cells, boundary, types, periodic=periodic, default_patch=default_patch)


# ------------------------------------------------------------ maillages prêts à l'emploi
def rectangle_mesh(x0, x1, y0, y1, nx, ny, grading=(1.0, 1.0), names=None, types=None,
                   periodic=None) -> Mesh2D:
    """Rectangle structuré ; names = {'left':.., 'right':.., 'bottom':.., 'top':..}."""
    names = {**{"left": "left", "right": "right", "bottom": "bottom", "top": "top"},
             **(names or {})}
    types = types or {}
    faces = {"bottom": [0, 1], "right": [1, 2], "top": [2, 3], "left": [3, 0]}
    patches = {}
    for side, f in faces.items():
        patches.setdefault(names[side], {"type": types.get(names[side], "patch"), "faces": []})
        patches[names[side]]["faces"].append(f)
    return block_mesh([[x0, y0], [x1, y0], [x1, y1], [x0, y1]],
                      [{"vertices": [0, 1, 2, 3], "cells": [nx, ny], "grading": list(grading)}],
                      patches=patches, periodic=periodic)


def channel_mesh(length: float = 1.0, height: float = 2.0, nx: int = 4, ny: int = 64,
                 first_height: float | None = None, periodic: bool = True) -> Mesh2D:
    """Canal plan [0, L] × [0, H], resserré aux deux parois ('bottom', 'top' = wall)."""
    half = ny // 2
    if first_height is None:
        gy = 1.0
    else:
        r = grading_for_first_cell(half, first_height / (0.5 * height))
        gy = [(0.5, 0.5, r), (0.5, 0.5, 1.0 / r)]
    names = {"left": "inlet", "right": "outlet", "bottom": "bottom", "top": "top"}
    return rectangle_mesh(0.0, length, 0.0, height, nx, ny, grading=(1.0, gy), names=names,
                          types={"bottom": "wall", "top": "wall"},
                          periodic=[("inlet", "outlet")] if periodic else None)


def cavity_mesh(n: int = 64, grading: float = 1.0) -> Mesh2D:
    """Cavité carrée [0,1]² : 'lid' (haut, mobile) et 'walls' ; grading > 1 resserre aux parois."""
    g = [(0.5, 0.5, grading), (0.5, 0.5, 1.0 / grading)] if grading != 1.0 else 1.0
    return rectangle_mesh(0, 1, 0, 1, n, n, grading=(g, g),
                          names={"left": "walls", "right": "walls", "bottom": "walls", "top": "lid"},
                          types={"walls": "wall", "lid": "wall"})


def backward_facing_step_mesh(step: float = 1.0, upstream: float = 4.0, downstream: float = 30.0,
                              channel: float = 1.0, n_step: int = 20, n_channel: int = 20,
                              n_up: int = 30, n_down: int = 120,
                              wall_grading: float = 4.0) -> Mesh2D:
    """Marche descendante (3 blocs) : entrée à gauche, marche de hauteur `step` en x = 0."""
    h, H = step, channel
    V = [[-upstream, h], [0, h], [0, h + H], [-upstream, h + H],      # 0..3 amont
         [0, 0], [downstream, 0], [downstream, h], [downstream, h + H]]  # 4..7
    gy_ch = [(0.5, 0.5, wall_grading), (0.5, 0.5, 1.0 / wall_grading)]
    gy_st = [(0.5, 0.5, wall_grading), (0.5, 0.5, 1.0 / wall_grading)]
    blocks = [
        {"vertices": [0, 1, 2, 3], "cells": [n_up, n_channel], "grading": [1.0 / wall_grading, gy_ch]},
        {"vertices": [1, 6, 7, 2], "cells": [n_down, n_channel], "grading": [8.0, gy_ch]},
        {"vertices": [4, 5, 6, 1], "cells": [n_down, n_step], "grading": [8.0, gy_st]},
    ]
    patches = {
        "inlet": {"type": "patch", "faces": [[3, 0]]},
        "outlet": {"type": "patch", "faces": [[6, 7], [5, 6]]},
        "walls": {"type": "wall", "faces": [[0, 1], [1, 4], [4, 5], [2, 3], [7, 2]]},
    }
    return block_mesh(V, blocks, patches=patches)


def flat_plate_mesh(length: float = 2.0, upstream: float = 0.333, height: float = 1.0,
                    nx_up: int = 16, nx_plate: int = 96, ny: int = 64,
                    first_height: float = 5e-6, le_fraction: float = 0.002) -> Mesh2D:
    """Plaque plane (géométrie du NASA TMR) : symétrie amont, plaque 'wall' de x=0 à L."""
    ry = grading_for_first_cell(ny, first_height / height)
    rx = grading_for_first_cell(nx_plate, le_fraction * length / length)
    V = [[-upstream, 0], [0, 0], [length, 0], [length, height], [0, height], [-upstream, height]]
    blocks = [
        {"vertices": [0, 1, 4, 5], "cells": [nx_up, ny], "grading": [0.1, ry]},
        {"vertices": [1, 2, 3, 4], "cells": [nx_plate, ny], "grading": [rx, ry]},
    ]
    patches = {
        "inlet": {"type": "patch", "faces": [[5, 0]]},
        "outlet": {"type": "patch", "faces": [[2, 3]]},
        "top": {"type": "patch", "faces": [[3, 4], [4, 5]]},
        "symmetry": {"type": "symmetry", "faces": [[0, 1]]},
        "plate": {"type": "wall", "faces": [[1, 2]]},
    }
    return block_mesh(V, blocks, patches=patches)
