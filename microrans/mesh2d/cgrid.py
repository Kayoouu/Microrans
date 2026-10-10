"""Maillage structuré en C autour d'un profil (bord de fuite fermé ou pointu).

La ligne de maillage j = 0 part de la sortie, longe le sillage (côté intrados) jusqu'au bord
de fuite, fait le tour du profil (intrados, bord d'attaque, extrados), revient au bord de
fuite et longe le sillage (côté extrados) jusqu'à la sortie. Les deux bords de la coupure de
sillage sont fusionnés : faces internes, le sillage est maillé finement sur toute sa
longueur, contrairement au maillage en O où les lignes radiales s'écartent derrière le profil.

Comme pour le maillage en O : couches orthogonales extrudées le long de la normale près de la
paroi (1re maille imposée), raccordées progressivement au bord extérieur ; ce bord est un
demi-cercle de rayon R centré au bord de fuite, prolongé par deux droites y = ±R jusqu'à la
sortie x = x_BF + wake_length.
"""
from __future__ import annotations

import numpy as np

from .geometry import Shape, signed_area
from .mesh import Mesh2D
from .ogrid import geometric_layers


def _open_normals(curve: np.ndarray, smooth: int = 5) -> np.ndarray:
    """Normales à gauche d'une courbe ouverte parcourue dans l'ordre des points."""
    t = np.empty_like(curve)
    t[1:-1] = curve[2:] - curve[:-2]
    t[0], t[-1] = curve[1] - curve[0], curve[-1] - curve[-2]
    t /= np.linalg.norm(t, axis=1)[:, None]
    n = np.column_stack([-t[:, 1], t[:, 0]])
    for _ in range(smooth):
        n[1:-1] = 0.25 * n[:-2] + 0.5 * n[1:-1] + 0.25 * n[2:]
        n /= np.linalg.norm(n, axis=1)[:, None]
    return n


def _layers(n: int, first: np.ndarray) -> np.ndarray:
    """geometric_layers(n, first[i], 1) pour chaque ligne i, vectorisé (bissection sur la
    raison q de la progression)."""
    first = np.asarray(first, float)
    lo, hi = np.full_like(first, 1.0 + 1e-12), np.full_like(first, 10.0)
    for _ in range(200):
        q = 0.5 * (lo + hi)
        below = np.log(first) + n * np.log(q) - np.log(q - 1) < np.log(1.0 + first / (q - 1))
        lo, hi = np.where(below, q, lo), np.where(below, hi, q)
    q = 0.5 * (lo + hi)[:, None]
    r = first[:, None] * (q ** np.arange(n + 1) - 1) / (q - 1)
    r[:, -1] = 1.0
    uniform = first * n >= 1.0
    r[uniform] = np.linspace(0.0, 1.0, n + 1)
    return r


def trailing_edge_index(xb: np.ndarray, max_angle: float = 60.0) -> int:
    """Indice du bord de fuite (point d'abscisse maximale) d'un contour fermé ; refuse un
    bord de fuite épais ou arrondi (angle entre les deux côtés > max_angle degrés)."""
    i = int(np.argmax(xb[:, 0]))
    a, b = xb[i - 1] - xb[i], xb[(i + 1) % len(xb)] - xb[i]
    ang = np.degrees(np.arccos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)),
                                       -1.0, 1.0)))
    if ang > max_angle:
        raise ValueError(
            f"Maillage en C : pas de bord de fuite pointu (angle {ang:.0f}° au point "
            f"d'abscisse maximale). Profils à bord de fuite fermé ou pointu seulement "
            "(NACA : trailing_edge = \"closed\" ou \"sharp\") ; sinon maillage en O ou hybride.")
    return i


def c_grid(body: Shape, n_around: int = 128, n_radial: int = 64, n_wake: int = 32,
           farfield_radius: float = 20.0, wake_length: float | None = None,
           first_height: float = 1e-5, wall_name: str | None = None,
           farfield_name: str = "farfield", index_blend: float = 0.5,
           normal_smoothing: int = 5, curve: np.ndarray | None = None) -> Mesh2D:
    """Maillage en C : (n_around + 2 n_wake) × n_radial cellules.

    n_around : mailles sur le profil ; n_wake : mailles le long du sillage (de chaque côté
    de la coupure) ; farfield_radius : rayon R du demi-cercle extérieur, centré au bord de
    fuite ; wake_length : longueur du sillage maillé (défaut R) ; first_height : épaisseur
    de la 1re couche (paroi et coupure de sillage) ; index_blend : comme pour le maillage
    en O (0 : points extérieurs répartis comme ceux du profil, 1 : uniformément)."""
    xb = np.asarray(curve, float) if curve is not None else body.boundary_curve(n=n_around)
    if signed_area(xb) < 0:
        xb = xb[::-1]
    xb = np.roll(xb, -trailing_edge_index(xb), axis=0)       # contour direct, BF en tête
    te = xb[0]
    # profil en sens horaire : BF → intrados → bord d'attaque → extrados → BF
    air = np.vstack([xb[:1], xb[:0:-1], xb[:1]])
    ds_te = 0.5 * (np.linalg.norm(air[1] - air[0]) + np.linalg.norm(air[-1] - air[-2]))
    R = float(farfield_radius)
    Lw = float(wake_length) if wake_length is not None else R
    s = geometric_layers(n_wake, ds_te, Lw)[1:]               # abscisses du sillage, BF exclu
    wake = te + np.column_stack([s, np.zeros_like(s)])
    C = np.vstack([wake[::-1], air, wake])                    # ligne j = 0
    na = len(air)
    # bord extérieur correspondant : droite y = y_BF − R, demi-cercle, droite y = y_BF + R
    seg = np.linalg.norm(np.diff(air, axis=0), axis=1)
    frac = np.concatenate([[0.0], np.cumsum(seg)]) / seg.sum()
    frac = (1 - index_blend) * frac + index_blend * np.linspace(0.0, 1.0, na)
    th = -0.5 * np.pi - np.pi * frac
    # sillage : même mélange, sinon les lignes resserrées au bord de fuite restent serrées
    # jusqu'au bord extérieur (cellules de rapport d'aspect ~1000, non-orthogonalité 86°)
    xo = te[0] + (1 - index_blend) * s + index_blend * Lw * np.arange(1, n_wake + 1) / n_wake
    outer = np.vstack([np.column_stack([xo[::-1], np.full(n_wake, te[1] - R)]),
                       te + R * np.column_stack([np.cos(th), np.sin(th)]),
                       np.column_stack([xo, np.full(n_wake, te[1] + R)])])
    nrm = _open_normals(C, normal_smoothing)
    nrm[:n_wake] = (0.0, -1.0)                                # sillage : normales exactes
    nrm[-n_wake:] = (0.0, 1.0)
    # Près du bord de fuite, profil et sillage forment un coin rentrant (172° côté fluide
    # pour le NACA 0012) : les normales y convergent et les lignes se croisent. On impose
    # une rotation monotone des normales le long de la ligne j = 0 (lignes parallèles ou
    # divergentes) : jusqu'au bord d'attaque, l'angle ne fait que décroître ; ensuite, en
    # partant de la sortie côté extrados, il ne fait que croître.
    phi = np.unwrap(np.arctan2(nrm[:, 1], nrm[:, 0]))
    mid = n_wake + na // 2
    phi[:mid] = np.minimum.accumulate(phi[:mid])
    phi[mid:] = np.maximum.accumulate(phi[mid:][::-1])[::-1]
    nrm = np.column_stack([np.cos(phi), np.sin(phi)])
    L = np.linalg.norm(outer - C, axis=1)
    # progression propre à chaque ligne : 1re maille = first_height partout (une épaisseur
    # commune en fraction de L donnait des mailles voisines d'épaisseurs différentes, d'où
    # une asymétrie de 300 sur la coupure de sillage)
    eta = _layers(n_radial, first_height / L)                 # (lignes, n_radial + 1)
    blend = 3 * eta ** 2 - 2 * eta ** 3
    P = C[:, None] + (eta * L[:, None])[..., None] * nrm[:, None]
    Q = C[:, None] + eta[..., None] * (outer - C)[:, None]
    X = (1 - blend)[..., None] * P + blend[..., None] * Q
    n = len(C)
    node = np.arange(n * (n_radial + 1)).reshape(n, n_radial + 1)
    # coupure de sillage : les points j = 0 du côté extrados (et le BF) sont ceux de l'intrados
    k = np.arange(n_wake + 1)
    node[n - 1 - k, 0] = node[k, 0]
    used, inv = np.unique(node, return_inverse=True)
    node = inv.reshape(node.shape)
    pts = X.reshape(-1, 2)[used]
    cells = np.stack([node[:-1, :-1], node[1:, :-1], node[1:, 1:], node[:-1, 1:]],
                     axis=-1).reshape(-1, 4)
    q = pts[cells]
    x, y = q[..., 0], q[..., 1]
    area = 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)
    if np.any(np.sign(area) != np.sign(np.median(area))) or np.any(area == 0.0):
        raise ValueError("Maillage en C : cellules retournées (first_height trop grand, "
                         "profil trop cambré ou trop peu de points). Essayer le maillage en O "
                         "ou hybride.")
    i1, i2 = n_wake, n_wake + na - 1                          # bords de fuite sur la ligne j = 0
    wall = wall_name or body.name
    boundary = {
        wall: np.column_stack([node[i1:i2, 0], node[i1 + 1:i2 + 1, 0]]),
        farfield_name: np.vstack([
            np.column_stack([node[:-1, -1], node[1:, -1]]),          # bord en C
            np.column_stack([node[0, :-1], node[0, 1:]]),            # sortie, côté intrados
            np.column_stack([node[-1, :-1], node[-1, 1:]])]),        # sortie, côté extrados
    }
    return Mesh2D(pts, cells, boundary, {wall: "wall", farfield_name: "patch"})
