"""Maillage 2D non structuré pour les volumes finis (cellules polygonales quelconques).

Organisation « à la OpenFOAM » :
- faces internes d'abord (owner < neighbour, triées), puis faces frontières groupées par patch ;
- le vecteur surface S_f d'une face pointe de l'owner vers le neighbour (ou vers l'extérieur) ;
- les patches ont un type géométrique : wall, patch, symmetry, empty (+ périodicité par paires,
  convertie en faces internes avec un vecteur de translation).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

PATCH_TYPES = ("wall", "patch", "symmetry", "empty")


@dataclass
class Patch:
    name: str
    type: str
    start: int      # indice global de la première face
    size: int

    @property
    def faces(self) -> slice:
        return slice(self.start, self.start + self.size)


def _pad_cells(cells, n_points):
    if isinstance(cells, np.ndarray) and cells.ndim == 2:
        cn = cells.astype(np.int64)
        nv = np.sum(cn >= 0, axis=1)
    else:
        cells = [np.asarray(c, dtype=np.int64) for c in cells]
        m = max(len(c) for c in cells)
        cn = -np.ones((len(cells), m), dtype=np.int64)
        nv = np.array([len(c) for c in cells])
        for i, c in enumerate(cells):
            cn[i, :len(c)] = c
    if np.any(nv < 3):
        raise ValueError("Toutes les cellules doivent avoir au moins 3 sommets.")
    if cn.max() >= n_points:
        raise ValueError("Indice de sommet hors limites dans les cellules.")
    return cn, nv


def _edge_key(a, b, n):
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    return lo * n + hi


class Mesh2D:
    """Maillage 2D volumes finis.

    Paramètres
    ----------
    points : (Np, 2) coordonnées des sommets.
    cells : liste de polygones (indices de sommets) ou tableau (Nc, k) complété par -1.
    boundary : {nom_patch: arêtes (Ne, 2)} ; les arêtes frontières non listées vont dans
        `default_patch` (type wall), avec un avertissement.
    patch_types : {nom_patch: type} parmi wall, patch, symmetry, empty (défaut : patch).
    periodic : liste de paires (patch_A, patch_B) ou (patch_A, patch_B, (tx, ty)) ; la
        translation (A -> B) est déduite des centres si elle n'est pas donnée.
    """

    dim = 2

    def __init__(self, points, cells, boundary: dict | None = None,
                 patch_types: dict | None = None, periodic=None,
                 default_patch: str = "defaultFaces"):
        self.points = np.asarray(points, dtype=float).reshape(-1, 2)
        n_pts = len(self.points)
        cn, nv = _pad_cells(cells, n_pts)
        self.cell_nodes, self.cell_nv = self._orient(cn, nv)
        self.n_cells = len(nv)
        boundary = {k: np.asarray(v, dtype=np.int64).reshape(-1, 2)
                    for k, v in (boundary or {}).items()}
        patch_types = dict(patch_types or {})
        self._build_faces(boundary, patch_types, periodic or [], default_patch)
        self._build_geometry()
        self._wall_distance = None

    # ------------------------------------------------------------------ topologie
    def _orient(self, cn, nv):
        idx = np.arange(cn.shape[1])[None, :]
        valid = idx < nv[:, None]
        nxt = np.where(idx + 1 < nv[:, None], idx + 1, 0)
        rows = np.arange(len(nv))[:, None]
        p = self.points
        a = np.where(valid, cn, 0)
        b = np.where(valid, cn[rows, nxt], 0)
        cross = np.where(valid, p[a, 0] * p[b, 1] - p[b, 0] * p[a, 1], 0.0)
        area = 0.5 * cross.sum(axis=1)
        if np.any(np.abs(area) < 1e-300):
            raise ValueError("Cellule(s) d'aire nulle.")
        flip = area < 0
        if np.any(flip):
            rev = np.where(valid, nv[:, None] - 1 - idx, idx)
            cn = np.where(flip[:, None], cn[rows, rev], cn)
        return cn, nv

    def _half_edges(self):
        cn, nv = self.cell_nodes, self.cell_nv
        idx = np.arange(cn.shape[1])[None, :]
        valid = idx < nv[:, None]
        nxt = np.where(idx + 1 < nv[:, None], idx + 1, 0)
        rows = np.arange(len(nv))[:, None]
        a = cn[valid]
        b = cn[rows, nxt][valid]
        cell = np.broadcast_to(rows, cn.shape)[valid]
        return a, b, cell

    def _build_faces(self, boundary, patch_types, periodic, default_patch):
        n_pts = len(self.points)
        a, b, cell = self._half_edges()
        keys = _edge_key(a, b, n_pts)
        order = np.argsort(keys, kind="stable")
        ks = keys[order]
        uniq, start, counts = np.unique(ks, return_index=True, return_counts=True)
        if np.any(counts > 2):
            raise ValueError("Maillage non manifold : une arête est partagée par plus de 2 cellules.")
        # faces internes
        pair = start[counts == 2]
        h0, h1 = order[pair], order[pair + 1]
        c0, c1 = cell[h0], cell[h1]
        own_h = np.where(c0 <= c1, h0, h1)
        owner_i = cell[own_h]
        neigh_i = np.where(c0 <= c1, c1, c0)
        fn_i = np.column_stack([a[own_h], b[own_h]])
        srt = np.lexsort((neigh_i, owner_i))
        owner_i, neigh_i, fn_i = owner_i[srt], neigh_i[srt], fn_i[srt]
        # faces frontières
        hb = order[start[counts == 1]]
        fn_b = np.column_stack([a[hb], b[hb]])
        owner_b = cell[hb]
        kb = keys[hb]
        names = list(boundary)
        tag = np.full(len(hb), -1)
        for i, name in enumerate(names):
            e = boundary[name]
            if len(e) == 0:
                continue
            ke = np.sort(_edge_key(e[:, 0], e[:, 1], n_pts))
            pos = np.searchsorted(ke, kb)
            pos = np.minimum(pos, len(ke) - 1)
            hit = (ke[pos] == kb) & (tag < 0)
            tag[hit] = i
        if np.any(tag < 0):
            n_def = int(np.sum(tag < 0))
            if default_patch not in names:
                names.append(default_patch)
                patch_types.setdefault(default_patch, "wall")
            tag[tag < 0] = names.index(default_patch)
            warnings.warn(f"{n_def} arête(s) frontière sans patch -> '{default_patch}' (paroi).",
                          stacklevel=3)
        # périodicité : paires de patches converties en faces internes
        self.periodic_pairs = []
        shift_list, own_p, nei_p, fn_p = [], [], [], []
        remove = np.zeros(len(hb), dtype=bool)
        self._periodic_boundary = {}
        self._periodic_owner = {}
        for spec in periodic:
            pa, pb = spec[0], spec[1]
            ia, ib = names.index(pa), names.index(pb)
            fa, fb = np.nonzero(tag == ia)[0], np.nonzero(tag == ib)[0]
            if len(fa) != len(fb):
                raise ValueError(f"Périodicité {pa}/{pb} : nombres de faces différents.")
            ca = 0.5 * (self.points[fn_b[fa, 0]] + self.points[fn_b[fa, 1]])
            cb = 0.5 * (self.points[fn_b[fb, 0]] + self.points[fn_b[fb, 1]])
            t = np.asarray(spec[2], float) if len(spec) > 2 else cb.mean(axis=0) - ca.mean(axis=0)
            from scipy.spatial import cKDTree
            dist, j = cKDTree(cb).query(ca + t)
            lmin = np.min(np.linalg.norm(self.points[fn_b[fa, 0]] - self.points[fn_b[fa, 1]], axis=1))
            if np.any(dist > 1e-6 * max(lmin, 1e-300) + 1e-12) or len(np.unique(j)) != len(j):
                raise ValueError(f"Périodicité {pa}/{pb} : faces non appariées (translation {t}).")
            own_p.append(owner_b[fa])
            nei_p.append(owner_b[fb[j]])
            fn_p.append(fn_b[fa])
            shift_list.append(np.broadcast_to(t, (len(fa), 2)))
            remove[fa] = True
            remove[fb] = True
            self.periodic_pairs.append((pa, pb, t))
            self._periodic_boundary[pa] = fn_b[fa]
            self._periodic_boundary[pb] = fn_b[fb[j]]        # même ordre que A
            self._periodic_owner[pa] = owner_b[fa]
            self._periodic_owner[pb] = owner_b[fb[j]]
        n_reg = len(owner_i)
        self.n_regular_internal = n_reg
        if own_p:
            owner_i = np.concatenate([owner_i] + own_p)
            neigh_i = np.concatenate([neigh_i] + nei_p)
            fn_i = np.vstack([fn_i] + fn_p)
        shift = np.zeros((len(owner_i), 2))
        if shift_list:
            shift[n_reg:] = np.vstack(shift_list)
        keep = ~remove
        fn_b, owner_b, tag = fn_b[keep], owner_b[keep], tag[keep]
        # regroupement par patch
        self.patches: list[Patch] = []
        n_int = len(owner_i)
        fn_blocks, own_blocks = [], []
        pos = n_int
        periodic_names = {p for pair in self.periodic_pairs for p in pair[:2]}
        for i, name in enumerate(names):
            if name in periodic_names:
                continue
            sel = np.nonzero(tag == i)[0]
            if len(sel) == 0 and name.startswith("_"):
                continue        # patch interne temporaire (raccord) devenu vide
            ptype = patch_types.get(name, "patch")
            if ptype not in PATCH_TYPES:
                raise ValueError(f"Type de patch inconnu '{ptype}' ({name}).")
            self.patches.append(Patch(name, ptype, pos, len(sel)))
            fn_blocks.append(fn_b[sel])
            own_blocks.append(owner_b[sel])
            pos += len(sel)
        self.n_internal = n_int
        self.n_faces = pos
        self.face_nodes = np.vstack([fn_i] + fn_blocks) if fn_blocks else fn_i
        self.owner = np.concatenate([owner_i] + own_blocks) if own_blocks else owner_i
        self.neighbour = neigh_i
        self.shift = shift
        self.patch_types = {p.name: p.type for p in self.patches}

    # ------------------------------------------------------------------ géométrie
    def _build_geometry(self):
        p = self.points
        cn, nv = self.cell_nodes, self.cell_nv
        idx = np.arange(cn.shape[1])[None, :]
        valid = idx < nv[:, None]
        nxt = np.where(idx + 1 < nv[:, None], idx + 1, 0)
        rows = np.arange(len(nv))[:, None]
        ia = np.where(valid, cn, 0)
        ib = np.where(valid, cn[rows, nxt], 0)
        xa, ya, xb, yb = p[ia, 0], p[ia, 1], p[ib, 0], p[ib, 1]
        cross = np.where(valid, xa * yb - xb * ya, 0.0)
        area = 0.5 * cross.sum(axis=1)
        cx = np.sum((xa + xb) * cross, axis=1) / (6.0 * area)
        cy = np.sum((ya + yb) * cross, axis=1) / (6.0 * area)
        self.cell_volumes = area
        self.cell_centers = np.column_stack([cx, cy])

        fa, fb = p[self.face_nodes[:, 0]], p[self.face_nodes[:, 1]]
        self.face_centers = 0.5 * (fa + fb)
        e = fb - fa
        self.Sf = np.column_stack([e[:, 1], -e[:, 0]])
        self.magSf = np.linalg.norm(self.Sf, axis=1)
        self.nf = self.Sf / self.magSf[:, None]

        ni = self.n_internal
        P, N = self.owner[:ni], self.neighbour
        C = self.cell_centers
        cN = C[N] - self.shift
        self.d_PN = cN - C[P]
        n_i = self.nf[:ni]
        dP = np.sum((self.face_centers[:ni] - C[P]) * n_i, axis=1)
        dN = np.sum((cN - self.face_centers[:ni]) * n_i, axis=1)
        if np.any(dP + dN <= 0):
            raise ValueError("Maillage invalide : centres de cellules du mauvais côté d'une face.")
        self.weights = dN / (dP + dN)                       # poids de l'owner
        dS = np.sum(self.d_PN * self.Sf[:ni], axis=1)
        self.orth_coeff = self.magSf[:ni] ** 2 / dS         # |S|²/(d·S)
        self.nonorth_vec = self.Sf[:ni] - self.orth_coeff[:, None] * self.d_PN
        # frontières
        Pb = self.owner[ni:]
        self.d_Pb = self.face_centers[ni:] - C[Pb]
        self.d_perp_b = np.sum(self.d_Pb * self.nf[ni:], axis=1)
        if np.any(self.d_perp_b <= 0):
            raise ValueError("Maillage invalide : centre de cellule hors du domaine (frontière).")

    # ------------------------------------------------------------------ accès
    def patch(self, name: str) -> Patch:
        for p in self.patches:
            if p.name == name:
                return p
        raise KeyError(f"Patch inconnu '{name}'. Disponibles : {[p.name for p in self.patches]}")

    def patch_face_nodes(self, name: str) -> np.ndarray:
        if name in self._periodic_boundary:
            return self._periodic_boundary[name]
        return self.face_nodes[self.patch(name).faces]

    def all_boundary_patches(self) -> list[tuple[str, str, np.ndarray]]:
        """(nom, type, arêtes) de tous les patches, y compris les périodiques (pour l'export)."""
        out = [(p.name, p.type, self.face_nodes[p.faces]) for p in self.patches]
        for pa, pb, _ in self.periodic_pairs:
            out.append((pa, "cyclic", self._periodic_boundary[pa]))
            out.append((pb, "cyclic", self._periodic_boundary[pb]))
        return out

    @property
    def n_points(self) -> int:
        return len(self.points)

    def cells_as_lists(self) -> list[np.ndarray]:
        return [row[:k] for row, k in zip(self.cell_nodes, self.cell_nv)]

    @property
    def wall_patches(self) -> list[str]:
        return [p.name for p in self.patches if p.type == "wall"]

    @property
    def wall_distance(self) -> np.ndarray:
        """Distance de chaque centre de cellule à la paroi la plus proche (patches 'wall')."""
        if self._wall_distance is None:
            self._wall_distance = self.distance_to_patches(self.wall_patches)
        return self._wall_distance

    def distance_to_patches(self, names) -> np.ndarray:
        return self._nearest_on_patches(names)[0]

    def wall_normal(self) -> np.ndarray:
        """Vecteur unitaire allant du point de paroi le plus proche vers chaque centre de
        cellule (= ∇d, normale locale à la paroi ; sert au modèle de transition γ)."""
        d, vec = self._nearest_on_patches(self.wall_patches)
        return vec / np.maximum(d, 1e-300)[:, None]

    def _nearest_on_patches(self, names):
        """(distance, vecteur point le plus proche → centre) aux segments des patches."""
        segs = [self.face_nodes[self.patch(n).faces] for n in names]
        if not segs or sum(len(s) for s in segs) == 0:
            vec = np.zeros((self.n_cells, 2))
            vec[:, 1] = 1e30
            return np.full(self.n_cells, 1e30), vec
        seg = np.vstack(segs)
        a, b = self.points[seg[:, 0]], self.points[seg[:, 1]]
        ab = b - a
        ab2 = np.maximum(np.sum(ab * ab, axis=1), 1e-300)
        out = np.empty(self.n_cells)
        vec = np.empty((self.n_cells, 2))
        C = self.cell_centers
        chunk = max(1, 2_000_000 // max(len(seg), 1))
        for s in range(0, self.n_cells, chunk):
            q = C[s:s + chunk, None, :]
            ap = q - a[None]
            t = np.clip(np.sum(ap * ab[None], axis=2) / ab2[None], 0.0, 1.0)
            r = ap - t[..., None] * ab[None]
            d2 = np.sum(r ** 2, axis=2)
            j = d2.argmin(axis=1)
            rows = np.arange(len(j))
            out[s:s + chunk] = np.sqrt(d2[rows, j])
            vec[s:s + chunk] = r[rows, j]
        return out, vec

    @property
    def first_cell_height(self) -> float:
        walls = self.wall_patches
        if not walls:
            return float(np.sqrt(self.cell_volumes.min()))
        sel = np.concatenate([np.arange(self.n_faces)[self.patch(n).faces] for n in walls])
        return float(self.d_perp_b[sel - self.n_internal].min())

    def bbox(self):
        mn, mx = self.points.min(axis=0), self.points.max(axis=0)
        return (mn[0], mn[1], mx[0], mx[1])

    # ------------------------------------------------------------------ qualité
    def quality(self) -> dict:
        """Indicateurs de qualité inspirés de `checkMesh` (OpenFOAM)."""
        ni = self.n_internal
        d = self.d_PN
        S = self.Sf[:ni]
        cosang = np.sum(d * S, axis=1) / (np.linalg.norm(d, axis=1) * self.magSf[:ni])
        nonorth = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
        # asymétrie : distance entre le centre de face et l'intersection de la droite PN
        C = self.cell_centers
        lam = np.sum((self.face_centers[:ni] - C[self.owner[:ni]]) * self.nf[:ni], axis=1) / \
            np.sum(d * self.nf[:ni], axis=1)
        xi = C[self.owner[:ni]] + lam[:, None] * d
        skew = np.linalg.norm(self.face_centers[:ni] - xi, axis=1) / np.linalg.norm(d, axis=1)
        # rapport d'aspect : plus longue / plus courte arête de chaque cellule
        lf = self.magSf
        allf = np.concatenate([self.owner, self.neighbour])
        allL = np.concatenate([lf, lf[:ni]])
        lmax = np.full(self.n_cells, 0.0)
        lmin = np.full(self.n_cells, np.inf)
        np.maximum.at(lmax, allf, allL)
        np.minimum.at(lmin, allf, allL)
        aspect = lmax / lmin
        types = {3: "triangles", 4: "quadrilatères"}
        counts = {}
        for k, c in zip(*np.unique(self.cell_nv, return_counts=True)):
            counts[types.get(int(k), f"polygones à {k} côtés")] = int(c)
        return {
            "n_points": int(self.n_points), "n_cells": int(self.n_cells),
            "n_faces": int(self.n_faces), "n_internal_faces": int(ni),
            "cell_types": counts,
            "patches": {p.name: {"type": p.type, "faces": p.size} for p in self.patches}
            | {pa: {"type": "cyclic", "faces": len(self._periodic_boundary[pa])}
               for pair in self.periodic_pairs for pa in pair[:2]},
            "area_min": float(self.cell_volumes.min()), "area_max": float(self.cell_volumes.max()),
            "non_orthogonality_max_deg": float(nonorth.max()) if ni else 0.0,
            "non_orthogonality_mean_deg": float(nonorth.mean()) if ni else 0.0,
            "skewness_max": float(skew.max()) if ni else 0.0,
            "aspect_ratio_max": float(aspect.max()),
            "first_center_wall_distance": self.first_cell_height,
        }

    def check(self) -> list[str]:
        """Liste d'avertissements de qualité (seuils usuels d'OpenFOAM)."""
        q = self.quality()
        msgs = []
        if q["non_orthogonality_max_deg"] > 70:
            msgs.append(f"non-orthogonalité max {q['non_orthogonality_max_deg']:.1f}° > 70°")
        if q["skewness_max"] > 4:
            msgs.append(f"asymétrie max {q['skewness_max']:.2f} > 4")
        return msgs

    def __repr__(self):
        return (f"Mesh2D({self.n_cells} cellules, {self.n_points} sommets, "
                f"patches={[p.name for p in self.patches]})")

    # ------------------------------------------------------------------ assemblage
    def set_patch_types(self, types: dict):
        """Change le type de patches existants (ex. 'wall' imposé par le cas de calcul)."""
        for p in self.patches:
            if p.name in types:
                if types[p.name] not in PATCH_TYPES:
                    raise ValueError(f"Type de patch inconnu '{types[p.name]}'.")
                p.type = types[p.name]
        self.patch_types = {p.name: p.type for p in self.patches}
        self._wall_distance = None

    def subset(self, keep, new_patch: str = "cut", new_type: str = "patch") -> "Mesh2D":
        """Sous-maillage des cellules `keep` (masque booléen) ; les faces coupées forment le
        patch `new_patch` (comme subsetMesh d'OpenFOAM)."""
        keep = np.asarray(keep, dtype=bool)
        if self.periodic_pairs:
            raise ValueError("subset : maillages périodiques non pris en charge.")
        if not keep.any():
            raise ValueError("subset : aucune cellule conservée.")
        ni = self.n_internal
        bnd = {}
        for p in self.patches:
            f = np.arange(p.start, p.start + p.size)
            bnd[p.name] = self.face_nodes[f[keep[self.owner[f]]]]
        P, N = self.owner[:ni], self.neighbour
        cut = keep[P] != keep[N]
        edges = self.face_nodes[:ni][cut]
        # arête orientée comme vue depuis la cellule conservée
        flip = ~keep[P][cut]
        edges[flip] = edges[flip][:, ::-1]
        bnd[new_patch] = np.vstack([bnd.get(new_patch, np.zeros((0, 2), int)), edges])
        cells = [c for c, k in zip(self.cells_as_lists(), keep) if k]
        used = np.unique(np.concatenate(cells))
        remap = -np.ones(self.n_points, dtype=np.int64)
        remap[used] = np.arange(len(used))
        types = {p.name: p.type for p in self.patches}
        types[new_patch] = new_type
        return Mesh2D(self.points[used], [remap[c] for c in cells],
                      {k: remap[v] for k, v in bnd.items() if len(v)}, types)

    def cut_at_axis(self, name: str = "axis") -> "Mesh2D":
        """Moitié y > 0 d'un maillage symétrique par rapport à l'axe x (calcul
        axisymétrique) ; les sommets à |y| < 1e-9·taille sont placés exactement sur l'axe et
        la coupe forme le patch `name` (type symmetry)."""
        tol = 1e-9 * max(float(np.ptp(self.points, axis=0).max()), 1e-300)
        m = self.subset(self.cell_centers[:, 1] > 0.0, name, "symmetry")
        on_axis = np.abs(m.points[:, 1]) < tol
        if not np.all(on_axis[np.unique(m.patch_face_nodes(name))]):
            raise ValueError("cut_at_axis : la coupe ne suit pas des arêtes du maillage sur "
                             "y = 0 (utiliser un nombre pair de points autour du corps, "
                             "symétrique par rapport à l'axe).")
        pts = m.points.copy()
        pts[on_axis, 1] = 0.0
        return Mesh2D(pts, m.cells_as_lists(), {p.name: m.face_nodes[p.faces]
                                               for p in m.patches}, m.patch_types)

    @staticmethod
    def merge(meshes: list["Mesh2D"], tol: float = 1e-9, patch_types=None) -> "Mesh2D":
        """Fusionne des maillages en confondant les sommets coïncidents (arêtes communes -> internes)."""
        pts, cells, bnd, types = [], [], {}, {}
        off = 0
        for m in meshes:
            pts.append(m.points)
            cells += [c + off for c in m.cells_as_lists()]
            for name, ptype, edges in m.all_boundary_patches():
                bnd.setdefault(name, []).append(edges + off)
                types.setdefault(name, "patch" if ptype == "cyclic" else ptype)
            off += m.n_points
        P, new = merge_points(np.vstack(pts), tol)
        cells = [new[c] for c in cells]
        bnd = {k: new[np.vstack(v)] for k, v in bnd.items()}
        types.update(patch_types or {})
        return Mesh2D(P, cells, bnd, types)


def merge_points(P: np.ndarray, tol: float = 1e-9):
    """Confond les points distants de moins de tol × taille du domaine.

    Retourne (points_uniques, nouvel_indice_de_chaque_point).
    """
    from scipy.spatial import cKDTree
    scale = max(np.ptp(P[:, 0]), np.ptp(P[:, 1]), 1e-300)
    pairs = cKDTree(P).query_pairs(tol * scale, output_type="ndarray")
    parent = np.arange(len(P))
    if len(pairs):
        # union-find itératif vectorisé : propage le plus petit indice
        a, b = pairs[:, 0], pairs[:, 1]
        while True:
            m = np.minimum(parent[a], parent[b])
            changed = (parent[a] != m) | (parent[b] != m)
            parent[a] = np.minimum(parent[a], m)
            parent[b] = np.minimum(parent[b], m)
            parent = parent[parent]
            if not np.any(changed):
                break
    uniq, new = np.unique(parent, return_inverse=True)
    return P[uniq], new
