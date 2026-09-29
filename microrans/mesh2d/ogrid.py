"""Maillage structuré en O autour d'un corps (cylindre, ellipse, profil, contour importé).

Principe : les points du corps sont extrudés le long de la normale (couches orthogonales près
de la paroi, 1re maille imposée pour contrôler y⁺), puis raccordés progressivement à un cercle
de champ lointain. Progression radiale géométrique.
"""
from __future__ import annotations

import numpy as np

from .geometry import Shape, signed_area
from .mesh import Mesh2D


def _normals(curve: np.ndarray, smooth: int = 5) -> np.ndarray:
    t = np.roll(curve, -1, axis=0) - np.roll(curve, 1, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    n = np.column_stack([t[:, 1], -t[:, 0]])      # normale extérieure (contour CCW)
    for _ in range(smooth):
        n = 0.25 * np.roll(n, 1, axis=0) + 0.5 * n + 0.25 * np.roll(n, -1, axis=0)
        n /= np.linalg.norm(n, axis=1)[:, None]
    return n


def geometric_layers(n: int, first: float, total: float) -> np.ndarray:
    """n+1 abscisses de 0 à total, 1re épaisseur `first`, progression géométrique."""
    if first * n >= total:
        return np.linspace(0.0, total, n + 1)
    lo, hi = 1.0 + 1e-12, 10.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        if first * (q ** n - 1) / (q - 1) < total:
            lo = q
        else:
            hi = q
    q = 0.5 * (lo + hi)
    r = first * (q ** np.arange(n + 1) - 1) / (q - 1)
    r[-1] = total
    return r


def o_grid(body: Shape, n_around: int = 128, n_radial: int = 64, farfield_radius: float = 20.0,
           first_height: float = 1e-3, center=None, wall_name: str | None = None,
           farfield_name: str = "farfield", index_blend: float = 0.5,
           normal_smoothing: int = 5, curve: np.ndarray | None = None) -> Mesh2D:
    """Maillage en O.

    farfield_radius : rayon absolu du cercle extérieur (centré sur le centre du corps).
    first_height : épaisseur de la 1re couche à la paroi.
    index_blend : 0 = points extérieurs à l'angle des points du corps, 1 = répartis
        uniformément sur le cercle ; 0.5 évite que le resserrement des points du corps
        (bords d'attaque/de fuite d'un profil) se propage jusqu'au champ lointain et croise
        les normales extrudées.
    """
    xb = np.asarray(curve, float) if curve is not None else body.boundary_curve(n=n_around)
    if signed_area(xb) < 0:
        xb = xb[::-1]
    n = len(xb)
    c = np.asarray(center, float) if center is not None else xb.mean(axis=0)
    rel = xb - c
    ang = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0]))
    if ang[-1] < ang[0]:
        raise ValueError("Contour non étoilé par rapport au centre : O-grid impossible "
                         "(utiliser le maillage non structuré ou hybride).")
    theta = ang[0] + (1 - index_blend) * (ang - ang[0]) + index_blend * 2 * np.pi * np.arange(n) / n
    outer = c + farfield_radius * np.column_stack([np.cos(theta), np.sin(theta)])
    nrm = _normals(xb, normal_smoothing)
    L = np.linalg.norm(outer - xb, axis=1)
    eta = geometric_layers(n_radial, first_height / L.mean(), 1.0)
    blend = 3 * eta ** 2 - 2 * eta ** 3
    P = xb[:, None] + eta[None, :, None] * L[:, None, None] * nrm[:, None]
    Q = xb[:, None] + eta[None, :, None] * (outer - xb)[:, None]
    X = (1 - blend)[None, :, None] * P + blend[None, :, None] * Q
    idx = np.arange(n * (n_radial + 1)).reshape(n, n_radial + 1)
    ip = np.roll(np.arange(n), -1)
    cells = np.stack([idx[:, :-1], idx[ip, :-1], idx[ip, 1:], idx[:, 1:]], axis=-1).reshape(-1, 4)
    # contrôle de validité : aucune cellule retournée
    pts = X.reshape(-1, 2)
    q = pts[cells]
    x, y = q[..., 0], q[..., 1]
    area = 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)
    if np.any(np.sign(area) != np.sign(np.median(area))):
        raise ValueError("O-grid : cellules retournées (corps trop concave ou first_height trop "
                         "grand). Essayer normal_smoothing plus grand ou le maillage hybride.")
    wall = wall_name or body.name
    boundary = {
        wall: np.column_stack([idx[:, 0], idx[ip, 0]]),
        farfield_name: np.column_stack([idx[:, -1], idx[ip, -1]]),
    }
    return Mesh2D(pts, cells, boundary, {wall: "wall", farfield_name: "patch"})
