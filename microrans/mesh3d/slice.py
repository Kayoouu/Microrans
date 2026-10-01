"""Coupe d'un maillage 3D en couches selon z (pavé, extrusion) par un plan z = constante,
pour les figures 2D (champs dans le plan médian) : les cellules de la couche coupée, vues
par leur face inférieure (polygone en x, y)."""
from __future__ import annotations

import numpy as np


class ZSlice:
    """Maillage 2D « vu » par les figures (points, cell_nodes, cell_nv, cell_centers, bbox,
    n_cells) ; `cells` : indices des cellules 3D coupées, `z` : cote du plan."""

    dim = 2

    def __init__(self, mesh, z=None):
        cn, nv = mesh.cell_nodes, mesh.cell_nv
        if not np.all(np.isin(nv, (6, 8))):
            raise ValueError("Coupe en z : prismes et hexaèdres empilés selon z seulement.")
        half = nv // 2
        rows = np.arange(len(cn))[:, None]
        k = np.arange(4)[None, :]
        lo = np.where(k < half[:, None], cn[rows, np.minimum(k, 3)], -1)
        hi = np.where(k < half[:, None], cn[rows, np.minimum(k + half[:, None], 7)], -1)
        P = mesh.points
        ok = lo >= 0
        same_xy = np.all(~ok[..., None] | (np.abs(P[np.maximum(lo, 0), :2]
                                                 - P[np.maximum(hi, 0), :2]) < 1e-9 * (
            np.ptp(P[:, :2], axis=0).max())), axis=(1, 2))
        if not same_xy.all():
            raise ValueError("Coupe en z : cellules non empilées selon z (maillage ni pavé ni "
                             "extrudé).")
        zlo = np.where(ok, P[np.maximum(lo, 0), 2], np.inf).min(axis=1)
        zhi = np.where(ok, P[np.maximum(hi, 0), 2], -np.inf).max(axis=1)
        if z is None:
            z = 0.5 * (P[:, 2].min() + P[:, 2].max())
        sel = np.nonzero((zlo <= z) & (z < zhi))[0]
        if not len(sel):                                  # plan sur la face supérieure
            sel = np.nonzero((zlo < z) & (z <= zhi))[0]
        if not len(sel):
            raise ValueError(f"plan z = {z:g} hors du domaine (z de {P[:, 2].min():g} à "
                             f"{P[:, 2].max():g}).")
        self.z = float(z)
        self.cells = sel
        used, inv = np.unique(lo[sel][ok[sel]], return_inverse=True)
        self.points = P[used, :2]
        nodes = -np.ones((len(sel), 4), dtype=np.int64)
        nodes[ok[sel]] = inv
        self.cell_nodes = nodes
        self.cell_nv = half[sel]
        self.n_cells = len(sel)
        self.cell_centers = mesh.cell_centers[sel, :2]

    def bbox(self):
        (x0, y0), (x1, y1) = self.points.min(axis=0), self.points.max(axis=0)
        return x0, y0, x1, y1
