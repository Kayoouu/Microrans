"""Coupes d'un maillage 3D par un plan, pour les figures 2D :

- `ZSlice` : plan z = constante d'un maillage en couches selon z (pavé, extrusion) ; les
  cellules de la couche coupée, vues par leur face inférieure (polygone en x, y) ;
- `PlaneSlice` : plan x = cte ou y = cte (ou z = cte) d'un maillage quelconque
  (hexaèdres, prismes, pyramides, tétraèdres) ; polygone d'intersection de chaque cellule
  coupée (arêtes coupées par le plan, sommets rangés par angle : cellules convexes).
- `slice_mesh(mesh, axis, value)` : ZSlice pour z, PlaneSlice pour x et y.

Les deux exposent ce que lisent les tracés (points, cell_nodes, cell_nv, cell_centers,
n_cells, bbox()) et `cells` (cellules 3D coupées), `axis`, `value`, `labels` (axes du plan)."""
from __future__ import annotations

import numpy as np

# arêtes des cellules selon leur nombre de sommets (numérotation VTK)
EDGES = {4: [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
         5: [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
         6: [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
         8: [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5),
             (2, 6), (3, 7)]}


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
        self.z = self.value = float(z)
        self.axis, self.labels = "z", ("x", "y")
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


class PlaneSlice:
    """Coupe par le plan {axis} = value (défaut : plan médian) d'un maillage 3D quelconque.
    Plan sur une face entre deux cellules : la cellule du côté des cotes croissantes (comme
    ZSlice) ; sur la face extrême du domaine : la cellule du bord. Axes de la figure :
    x -> (y, z), y -> (x, z), z -> (x, y)."""

    dim = 2

    def __init__(self, mesh, axis, value=None):
        a = "xyz".index(axis)
        u, v = [k for k in range(3) if k != a]
        P = mesh.points
        lo, hi = P[:, a].min(), P[:, a].max()
        if value is None:
            value = 0.5 * (lo + hi)
        if not lo <= value <= hi:
            raise ValueError(f"plan {axis} = {value:g} hors du domaine ({axis} de {lo:g} à "
                             f"{hi:g}).")
        cn, nv = np.asarray(mesh.cell_nodes), np.asarray(mesh.cell_nv)
        bad = sorted(set(np.unique(nv).tolist()) - set(EDGES))
        if bad:
            raise ValueError(f"Coupe : cellules à {bad} sommets non prises en charge.")
        s = P[:, a] - value
        sn = np.where(cn >= 0, s[np.maximum(cn, 0)], np.nan)
        smin, smax = np.nanmin(sn, axis=1), np.nanmax(sn, axis=1)
        above_strict = True                              # « dessus » : s > 0
        cand = np.nonzero((smin <= 0) & (smax > 0))[0]
        if not len(cand):                                # plan sur la face extrême haute
            above_strict = False                         # « dessus » : s ≥ 0
            cand = np.nonzero((smin < 0) & (smax >= 0))[0]
        tol = 1e-12 * float(np.ptp(P, axis=0).max())
        cells, polys = [], []
        for k in np.unique(nv[cand]):
            idx = cand[nv[cand] == k]
            E = np.array(EDGES[int(k)])
            A, B = cn[idx][:, E[:, 0]], cn[idx][:, E[:, 1]]          # (n, arêtes)
            sa, sb = s[A], s[B]
            up_a, up_b = (sa > 0, sb > 0) if above_strict else (sa >= 0, sb >= 0)
            cut = up_a != up_b
            t = np.where(cut, sa / np.where(cut, sa - sb, 1.0), 0.0)
            X = P[A][..., [u, v]] + t[..., None] * (P[B][..., [u, v]] - P[A][..., [u, v]])
            # sommet sur le plan relié à plusieurs arêtes coupées : un seul point
            d = np.max(np.abs(X[:, :, None, :] - X[:, None, :, :]), axis=3)
            ne = len(E)
            earlier = np.tril(np.ones((ne, ne), bool), -1)[None]
            dup = np.any((d <= tol) & cut[:, None, :] & earlier, axis=2)
            ok = cut & ~dup
            m = ok.sum(axis=1)
            keep = m >= 3                                # contact par une arête : pas d'aire
            idx, X, ok, m = idx[keep], X[keep], ok[keep], m[keep]
            c = np.sum(np.where(ok[..., None], X, 0.0), axis=1) / m[:, None]
            ang = np.where(ok, np.arctan2(X[..., 1] - c[:, None, 1], X[..., 0] - c[:, None, 0]),
                           np.inf)
            order = np.argsort(ang, axis=1)
            X = np.take_along_axis(X, order[..., None], axis=1)
            cells.append(idx)
            polys.append((X, m))
        if not cells:
            raise ValueError(f"plan {axis} = {value:g} : aucune cellule coupée.")
        cells_all = np.concatenate(cells)
        order = np.argsort(cells_all, kind="stable")
        width = max(int(m.max()) for _, m in polys)     # hexaèdre : 6 au plus
        pts, nodes, nvs, cent = [], [], [], []
        base = 0
        for X, m in polys:
            n = len(m)
            j = np.arange(X.shape[1])[None, :]
            valid = j < m[:, None]
            pts.append(X[valid])
            ids = base + (np.cumsum(valid.ravel()) - 1).reshape(valid.shape)
            row = -np.ones((n, width), dtype=np.int64)
            row[:, :X.shape[1]][valid[:, :width]] = ids[valid]
            nodes.append(row)
            nvs.append(m)
            cent.append(np.sum(np.where(valid[..., None], X, 0.0), axis=1) / m[:, None])
            base += int(valid.sum())
        self.axis, self.value = axis, float(value)
        self.labels = ("xyz"[u], "xyz"[v])
        self.cells = cells_all[order]
        self.points = np.concatenate(pts)
        self.cell_nodes = np.concatenate(nodes)[order]
        self.cell_nv = np.concatenate(nvs)[order]
        self.cell_centers = np.concatenate(cent)[order]
        self.n_cells = len(self.cells)

    def bbox(self):
        (x0, y0), (x1, y1) = self.points.min(axis=0), self.points.max(axis=0)
        return x0, y0, x1, y1


def slice_mesh(mesh, axis="z", value=None):
    """Coupe pour les figures : z -> ZSlice (maillages en couches selon z), x / y ->
    PlaneSlice."""
    if axis == "z":
        return ZSlice(mesh, value)
    return PlaneSlice(mesh, axis, value)
