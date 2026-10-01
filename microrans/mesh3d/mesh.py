"""Maillage 3D volumes finis (hexaèdres, prismes, tétraèdres, pyramides).

Même interface que Mesh2D (owner / neighbour, S_f, centres, volumes, poids, coefficients
d'orthogonalité, patches, périodicité, distance à la paroi, qualité) : les opérateurs
volumes finis et le solveur, écrits face par face, s'en servent tels quels, avec des
vecteurs à 3 composantes. Géométrie des faces et des cellules calculée comme OpenFOAM
(primitiveMesh) : face découpée en triangles autour de la moyenne de ses sommets, cellule
en pyramides de sommet la moyenne des centres de ses faces — exact pour les faces planes,
cohérent (Σ S_f = 0) pour les faces gauches.
"""
from __future__ import annotations

import warnings

import numpy as np

from ..mesh2d.mesh import PATCH_TYPES, Patch

# faces locales des cellules (numérotation VTK) ; l'orientation est corrigée ensuite
# (normale sortante) à partir de la géométrie, l'ordre de ces listes est indifférent
CELL_FACES = {
    4: [(0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1)],                          # tétraèdre
    5: [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],            # pyramide
    6: [(0, 1, 2), (3, 5, 4), (0, 3, 4, 1), (1, 4, 5, 2), (2, 5, 3, 0)],      # prisme
    8: [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1),
        (4, 5, 6, 7)],                                                         # hexaèdre
}
VTK_TYPES = {4: 10, 5: 14, 6: 13, 8: 12}


def face_geometry(points, fn, chunk: int = 250_000):
    """(centres, vecteurs surface) de faces polygonales (nf, 4) complétées par −1 (triangles).

    Triangles autour de la moyenne des sommets (OpenFOAM) : S = Σ ½ (x_i − x̄) × (x_i+1 − x̄),
    centre = moyenne des centres des triangles pondérée par leur aire (projetée sur S).
    Par paquets de `chunk` faces : les tableaux intermédiaires (nf, 4, 3) des 6 millions de
    faces locales d'un maillage de 10⁶ hexaèdres faisaient monter la mémoire à 5.2 Go."""
    if len(fn) > chunk:
        parts = [_face_geometry(points, fn[i:i + chunk]) for i in range(0, len(fn), chunk)]
        return np.vstack([a for a, _ in parts]), np.vstack([b for _, b in parts])
    return _face_geometry(points, fn)


def _face_geometry(points, fn):
    k = np.sum(fn >= 0, axis=1)
    fn0 = np.where(fn >= 0, fn, fn[:, :1])
    X = points[fn0]                                              # (nf, 4, 3)
    valid = (np.arange(fn.shape[1])[None, :] < k[:, None])
    xbar = np.sum(np.where(valid[..., None], X, 0.0), axis=1) / k[:, None]
    nxt = np.where(np.arange(fn.shape[1])[None, :] + 1 < k[:, None],
                   np.arange(fn.shape[1])[None, :] + 1, 0)
    Xn = np.take_along_axis(X, nxt[..., None], axis=1)
    tri = 0.5 * np.cross(X - xbar[:, None], Xn - xbar[:, None])     # (nf, 4, 3)
    tri = np.where(valid[..., None], tri, 0.0)
    S = tri.sum(axis=1)
    magS = np.linalg.norm(S, axis=1)
    nS = S / np.maximum(magS, 1e-300)[:, None]
    w = np.sum(tri * nS[:, None, :], axis=2)                     # aires projetées
    ctr = (X + Xn + xbar[:, None]) / 3.0
    fc = np.sum(w[..., None] * ctr, axis=1) / np.maximum(w.sum(axis=1), 1e-300)[:, None]
    return fc, S


def _flip(fn):
    """Inverse l'orientation (a, b, c, d) -> (a, d, c, b) ; triangles (a, b, c) -> (a, c, b)."""
    out = fn.copy()
    quad = fn[:, 3] >= 0
    out[quad, 1], out[quad, 3] = fn[quad, 3], fn[quad, 1]
    out[~quad, 1], out[~quad, 2] = fn[~quad, 2], fn[~quad, 1]
    return out


def _face_keys(fn):
    """Clé indépendante de l'ordre des sommets : sommets triés (−1 en tête pour les
    triangles)."""
    return np.sort(fn, axis=1)


class Mesh3D:
    """Maillage 3D volumes finis.

    points : (Np, 3) ; cells : tableau (Nc, 8) d'hexaèdres ou liste de cellules (4, 5, 6 ou
    8 sommets, numérotation VTK) ; boundary : {patch: faces (Nf, 3 ou 4)} (sommets de chaque
    face frontière, ordre quelconque) ; patch_types : {patch: wall | patch | symmetry |
    empty} ; periodic : [(patch_A, patch_B[, (tx, ty, tz)])].
    """

    dim = 3

    def __init__(self, points, cells, boundary: dict | None = None,
                 patch_types: dict | None = None, periodic=None,
                 default_patch: str = "defaultFaces"):
        self.points = np.asarray(points, dtype=float).reshape(-1, 3)
        if isinstance(cells, np.ndarray) and cells.ndim == 2:
            cl = [np.asarray(cells, dtype=np.int64)]
        else:
            groups = {}
            for c in cells:
                c = np.asarray(c, dtype=np.int64)
                groups.setdefault(len(c), []).append(c)
            cl = [np.array(v) for v in groups.values()]
        self._cells = []                     # [(indices globaux, sommets (n, k))]
        n = 0
        for arr in cl:
            if arr.shape[1] not in CELL_FACES:
                raise ValueError(f"Cellule à {arr.shape[1]} sommets : 4 (tétraèdre), 5 "
                                 "(pyramide), 6 (prisme) ou 8 (hexaèdre) attendus.")
            self._cells.append((np.arange(n, n + len(arr)), arr))
            n += len(arr)
        self.n_cells = n
        boundary = {k: np.asarray(v, dtype=np.int64).reshape(len(v), -1) if len(v) else
                    np.zeros((0, 4), dtype=np.int64) for k, v in (boundary or {}).items()}
        self._build_faces(boundary, dict(patch_types or {}), periodic or [], default_patch)
        self._build_geometry()
        self._wall_distance = None
        self._wall_vec = None

    # ------------------------------------------------------------------ topologie
    def _local_faces(self):
        fn, cell = [], []
        for idx, arr in self._cells:
            for f in CELL_FACES[arr.shape[1]]:
                a = -np.ones((len(arr), 4), dtype=np.int64)
                a[:, :len(f)] = arr[:, list(f)]
                fn.append(a)
                cell.append(idx)
        fn, cell = np.vstack(fn), np.concatenate(cell)
        # orientation sortante : normale dirigée du centre (moyenne des sommets) vers la face
        cc = np.zeros((self.n_cells, 3))
        for idx, arr in self._cells:
            cc[idx] = self.points[arr].mean(axis=1)
        fc, S = face_geometry(self.points, fn)
        inward = np.sum(S * (fc - cc[cell]), axis=1) < 0.0
        fn[inward] = _flip(fn[inward])
        return fn, cell

    def _build_faces(self, boundary, patch_types, periodic, default_patch):
        fn, cell = self._local_faces()
        keys = _face_keys(fn)
        uniq, inv, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
        inv = inv.ravel()
        if np.any(counts > 2):
            raise ValueError("Maillage non manifold : une face est partagée par plus de 2 "
                             "cellules.")
        order = np.argsort(inv, kind="stable")
        start = np.concatenate([[0], np.cumsum(counts)[:-1]])
        first = order[start]
        # faces internes
        two = counts == 2
        h0, h1 = first[two], order[start[two] + 1]
        c0, c1 = cell[h0], cell[h1]
        own_h = np.where(c0 <= c1, h0, h1)
        owner_i = cell[own_h]
        neigh_i = np.where(c0 <= c1, c1, c0)
        fn_i = fn[own_h]
        srt = np.lexsort((neigh_i, owner_i))
        owner_i, neigh_i, fn_i = owner_i[srt], neigh_i[srt], fn_i[srt]
        # faces frontières
        hb = first[~two]
        fn_b, owner_b = fn[hb], cell[hb]
        kb = keys[hb]
        names = list(boundary)
        tag = np.full(len(hb), -1)
        for i, name in enumerate(names):
            e = boundary[name]
            if len(e) == 0:
                continue
            if e.shape[1] == 3:
                e = np.column_stack([e, -np.ones(len(e), dtype=np.int64)])
            ke = _face_keys(e)
            both = np.vstack([ke, kb])
            _, iv = np.unique(both, axis=0, return_inverse=True)
            iv = iv.ravel()
            hit = np.isin(iv[len(ke):], iv[:len(ke)]) & (tag < 0)
            tag[hit] = i
        if np.any(tag < 0):
            n_def = int(np.sum(tag < 0))
            if default_patch not in names:
                names.append(default_patch)
                patch_types.setdefault(default_patch, "wall")
            tag[tag < 0] = names.index(default_patch)
            warnings.warn(f"{n_def} face(s) frontière sans patch -> '{default_patch}' "
                          "(paroi).", stacklevel=3)
        fcb, _ = face_geometry(self.points, fn_b)
        # périodicité
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
            t = (np.asarray(spec[2], float) if len(spec) > 2
                 else fcb[fb].mean(axis=0) - fcb[fa].mean(axis=0))
            from scipy.spatial import cKDTree
            dist, j = cKDTree(fcb[fb]).query(fcb[fa] + t)
            h = np.sqrt(np.min(np.linalg.norm(face_geometry(self.points, fn_b[fa])[1],
                                              axis=1)))          # longueur de maille
            if np.any(dist > 1e-6 * max(h, 1e-300) + 1e-12) or len(np.unique(j)) != len(j):
                raise ValueError(f"Périodicité {pa}/{pb} : faces non appariées (translation "
                                 f"{t}).")
            own_p.append(owner_b[fa])
            nei_p.append(owner_b[fb[j]])
            fn_p.append(fn_b[fa])
            shift_list.append(np.broadcast_to(t, (len(fa), 3)))
            remove[fa] = True
            remove[fb] = True
            self.periodic_pairs.append((pa, pb, t))
            self._periodic_boundary[pa] = fn_b[fa]
            self._periodic_boundary[pb] = fn_b[fb[j]]          # même ordre que A
            self._periodic_owner[pa] = owner_b[fa]
            self._periodic_owner[pb] = owner_b[fb[j]]
        n_reg = len(owner_i)
        self.n_regular_internal = n_reg
        if own_p:
            owner_i = np.concatenate([owner_i] + own_p)
            neigh_i = np.concatenate([neigh_i] + nei_p)
            fn_i = np.vstack([fn_i] + fn_p)
        shift = np.zeros((len(owner_i), 3))
        if shift_list:
            shift[n_reg:] = np.vstack(shift_list)
        keep = ~remove
        fn_b, owner_b, tag = fn_b[keep], owner_b[keep], tag[keep]
        self.patches: list[Patch] = []
        n_int = len(owner_i)
        fn_blocks, own_blocks = [], []
        pos = n_int
        periodic_names = {p for pair in self.periodic_pairs for p in pair[:2]}
        for i, name in enumerate(names):
            if name in periodic_names:
                continue
            sel = np.nonzero(tag == i)[0]
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
        fc, S = face_geometry(self.points, self.face_nodes)
        self.face_centers, self.Sf = fc, S
        self.magSf = np.linalg.norm(S, axis=1)
        if np.any(self.magSf < 1e-300):
            raise ValueError("Face(s) d'aire nulle.")
        self.nf = S / self.magSf[:, None]
        ni, nc = self.n_internal, self.n_cells
        P, N, Pb = self.owner[:ni], self.neighbour, self.owner[ni:]
        # centre approché : moyenne des centres de faces (côté neighbour : face translatée
        # pour les faces périodiques)
        allc = np.concatenate([P, N, Pb])
        fcs = np.vstack([fc[:ni], fc[:ni] + self.shift, fc[ni:]])
        cnt = np.bincount(allc, minlength=nc)
        cest = np.column_stack([np.bincount(allc, fcs[:, k], nc) for k in range(3)])
        cest /= cnt[:, None]
        # pyramides (sommet cest, base = face) : V = Σ S·(x_f − cest)/3 (normale sortante)
        sgn = np.concatenate([np.ones(ni), -np.ones(ni), np.ones(len(Pb))])
        Sall = np.vstack([S[:ni], S[:ni], S[ni:]]) * sgn[:, None]
        vp = np.sum(Sall * (fcs - cest[allc]), axis=1) / 3.0
        if np.any(vp <= 0.0):
            raise ValueError("Maillage invalide : pyramide de volume négatif (cellule "
                             "retournée ou trop déformée).")
        V = np.bincount(allc, vp, nc)
        cpy = 0.75 * fcs + 0.25 * cest[allc]
        C = np.column_stack([np.bincount(allc, vp * cpy[:, k], nc) for k in range(3)])
        self.cell_volumes = V
        self.cell_centers = C / V[:, None]
        C = self.cell_centers
        cN = C[N] - self.shift
        self.d_PN = cN - C[P]
        n_i = self.nf[:ni]
        dP = np.sum((fc[:ni] - C[P]) * n_i, axis=1)
        dN = np.sum((cN - fc[:ni]) * n_i, axis=1)
        if np.any(dP + dN <= 0):
            raise ValueError("Maillage invalide : centres de cellules du mauvais côté d'une "
                             "face.")
        self.weights = dN / (dP + dN)
        dS = np.sum(self.d_PN * S[:ni], axis=1)
        self.orth_coeff = self.magSf[:ni] ** 2 / dS
        self.nonorth_vec = S[:ni] - self.orth_coeff[:, None] * self.d_PN
        self.d_Pb = fc[ni:] - C[Pb]
        self.d_perp_b = np.sum(self.d_Pb * self.nf[ni:], axis=1)
        if np.any(self.d_perp_b <= 0):
            raise ValueError("Maillage invalide : centre de cellule hors du domaine "
                             "(frontière).")

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
        """(nom, type, faces) de tous les patches, y compris les périodiques (export)."""
        out = [(p.name, p.type, self.face_nodes[p.faces]) for p in self.patches]
        for pa, pb, _ in self.periodic_pairs:
            out.append((pa, "cyclic", self._periodic_boundary[pa]))
            out.append((pb, "cyclic", self._periodic_boundary[pb]))
        return out

    @property
    def n_points(self) -> int:
        return len(self.points)

    @property
    def wall_patches(self) -> list[str]:
        return [p.name for p in self.patches if p.type == "wall"]

    def cells_as_lists(self) -> list[np.ndarray]:
        out = [None] * self.n_cells
        for idx, arr in self._cells:
            for i, row in zip(idx, arr):
                out[i] = row
        return out

    @property
    def cell_nodes(self) -> np.ndarray:
        """Sommets des cellules (Nc, 8), complétés par −1 (comme Mesh2D.cell_nodes)."""
        if getattr(self, "_cell_nodes", None) is None:
            cn = -np.ones((self.n_cells, 8), dtype=np.int64)
            for idx, arr in self._cells:
                cn[idx, :arr.shape[1]] = arr
            self._cell_nodes = cn
        return self._cell_nodes

    @property
    def cell_nv(self) -> np.ndarray:
        nv = np.empty(self.n_cells, dtype=np.int64)
        for idx, arr in self._cells:
            nv[idx] = arr.shape[1]
        return nv

    def bbox(self):
        mn, mx = self.points.min(axis=0), self.points.max(axis=0)
        return tuple(mn) + tuple(mx)

    @property
    def wall_distance(self) -> np.ndarray:
        if self._wall_distance is None:
            self._wall_distance, self._wall_vec = self._nearest_on_patches(self.wall_patches)
        return self._wall_distance

    def distance_to_patches(self, names) -> np.ndarray:
        return self._nearest_on_patches(names)[0]

    def wall_normal(self) -> np.ndarray:
        d = self.wall_distance
        return self._wall_vec / np.maximum(d, 1e-300)[:, None]

    def _nearest_on_patches(self, names, k: int = 8):
        """(distance, vecteur point le plus proche → centre) aux faces des patches, exacte.

        Distance aux triangles (découpage des faces autour de la moyenne des sommets) des k
        faces de centres les plus proches ; une face plus lointaine est à au moins
        d_k − R (d_k : distance au k-ième centre, R : plus grand rayon de face), donc le
        résultat est exact si la meilleure distance trouvée est ≤ d_k − R. Sinon (faces très
        allongées), les faces de centre à moins de d + R sont toutes examinées."""
        fns = [self.patch_face_nodes(n) for n in names]
        nc = self.n_cells
        if not fns or sum(len(f) for f in fns) == 0:
            vec = np.zeros((nc, 3))
            vec[:, 2] = 1e30
            return np.full(nc, 1e30), vec
        fn = np.vstack(fns)
        fc, _ = face_geometry(self.points, fn)
        kk = np.sum(fn >= 0, axis=1)
        fn0 = np.where(fn >= 0, fn, fn[:, :1])
        X = self.points[fn0]                                     # (nf, 4, 3)
        xbar = np.sum(np.where((np.arange(4)[None, :] < kk[:, None])[..., None], X, 0.0),
                      axis=1) / kk[:, None]
        R = float(np.max(np.linalg.norm(X - fc[:, None, :], axis=2)))
        from scipy.spatial import cKDTree
        tree = cKDTree(fc)
        k = min(k, len(fn))
        dk, cand = tree.query(self.cell_centers, k=k)
        cand, dk = cand.reshape(nc, k), dk.reshape(nc, k)
        best, vec = self._dist_to_faces(self.cell_centers, cand, X, xbar, kk)
        unsure = np.nonzero(best > dk[:, -1] - R)[0] if k < len(fn) else np.zeros(0, int)
        # toutes les faces de centre à moins de d + R, par paires (cellule, face) ; par paquets
        # de cellules : les listes d'indices de query_ball_point (~150 entiers Python par
        # cellule loin des parois) prenaient ~5 Go à 1 million de cellules
        for u0 in range(0, len(unsure), 20_000):
            uc = unsure[u0:u0 + 20_000]
            balls = tree.query_ball_point(self.cell_centers[uc], best[uc] + R + 1e-12)
            cnt = np.array([len(b) for b in balls])
            cells = np.repeat(uc, cnt)
            faces = np.fromiter((f for b in balls for f in b), dtype=np.int64,
                                count=int(cnt.sum()))
            del balls
            step = 500_000
            for s0 in range(0, len(cells), step):
                c, f = cells[s0:s0 + step], faces[s0:s0 + step]
                b, v = self._dist_to_faces(self.cell_centers[c], f[:, None], X, xbar, kk)
                order = np.lexsort((b, c))                   # meilleure face par cellule
                c, b, v = c[order], b[order], v[order]
                first = np.concatenate([[True], c[1:] != c[:-1]])
                c, b, v = c[first], b[first], v[first]
                upd = b < best[c]
                best[c[upd]], vec[c[upd]] = b[upd], v[upd]
        return best, vec

    @staticmethod
    def _dist_to_faces(C, cand, X, xbar, kk):
        """Distance de chaque point C[i] aux faces cand[i, :] (triangles autour de xbar)."""
        n = len(C)
        best = np.full(n, np.inf)
        vec = np.zeros((n, 3))
        for c in range(cand.shape[1]):
            f = cand[:, c]
            for j in range(4):
                ok = j < kk[f]
                jn = np.where(j + 1 < kk[f], j + 1, 0)
                q = _closest_on_triangle(C, xbar[f], X[f, j], X[f, jn])
                r = C - q
                d = np.linalg.norm(r, axis=1)
                upd = ok & (d < best)
                best[upd] = d[upd]
                vec[upd] = r[upd]
        return best, vec

    @property
    def first_cell_height(self) -> float:
        walls = self.wall_patches
        if not walls:
            return float(np.cbrt(self.cell_volumes.min()))
        sel = np.concatenate([np.arange(self.n_faces)[self.patch(n).faces] for n in walls])
        return float(self.d_perp_b[sel - self.n_internal].min())

    def quality(self) -> dict:
        ni = self.n_internal
        d, S = self.d_PN, self.Sf[:ni]
        cosang = np.sum(d * S, axis=1) / (np.linalg.norm(d, axis=1) * self.magSf[:ni])
        nonorth = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
        C = self.cell_centers
        lam = np.sum((self.face_centers[:ni] - C[self.owner[:ni]]) * self.nf[:ni], axis=1) / \
            np.sum(d * self.nf[:ni], axis=1)
        xi = C[self.owner[:ni]] + lam[:, None] * d
        skew = np.linalg.norm(self.face_centers[:ni] - xi, axis=1) / np.linalg.norm(d, axis=1)
        lf = np.sqrt(self.magSf)
        allf = np.concatenate([self.owner, self.neighbour])
        allL = np.concatenate([lf, lf[:ni]])
        lmax = np.full(self.n_cells, 0.0)
        lmin = np.full(self.n_cells, np.inf)
        np.maximum.at(lmax, allf, allL)
        np.minimum.at(lmin, allf, allL)
        names = {4: "tétraèdres", 5: "pyramides", 6: "prismes", 8: "hexaèdres"}
        counts = {names[arr.shape[1]]: len(arr) for _, arr in self._cells}
        return {
            "n_points": int(self.n_points), "n_cells": int(self.n_cells),
            "n_faces": int(self.n_faces), "n_internal_faces": int(ni),
            "cell_types": counts,
            "patches": {p.name: {"type": p.type, "faces": p.size} for p in self.patches}
            | {pa: {"type": "cyclic", "faces": len(self._periodic_boundary[pa])}
               for pair in self.periodic_pairs for pa in pair[:2]},
            "volume_min": float(self.cell_volumes.min()),
            "volume_max": float(self.cell_volumes.max()),
            "non_orthogonality_max_deg": float(nonorth.max()) if ni else 0.0,
            "non_orthogonality_mean_deg": float(nonorth.mean()) if ni else 0.0,
            "skewness_max": float(skew.max()) if ni else 0.0,
            "aspect_ratio_max": float(np.max(lmax / lmin)),
            "first_center_wall_distance": self.first_cell_height,
        }

    def check(self) -> list[str]:
        q = self.quality()
        msgs = []
        if q["non_orthogonality_max_deg"] > 70:
            msgs.append(f"non-orthogonalité max {q['non_orthogonality_max_deg']:.1f}° > 70°")
        if q["skewness_max"] > 4:
            msgs.append(f"asymétrie max {q['skewness_max']:.2f} > 4")
        return msgs

    def set_patch_types(self, types: dict):
        walls_before = self.wall_patches
        for p in self.patches:
            if p.name in types:
                if types[p.name] not in PATCH_TYPES:
                    raise ValueError(f"Type de patch inconnu '{types[p.name]}'.")
                p.type = types[p.name]
        self.patch_types = {p.name: p.type for p in self.patches}
        if self.wall_patches != walls_before:      # distance à recalculer
            self._wall_distance = None
            self._wall_vec = None

    def __repr__(self):
        return (f"Mesh3D({self.n_cells} cellules, {self.n_points} sommets, "
                f"patches={[p.name for p in self.patches]})")


def _closest_on_triangle(p, a, b, c):
    """Point le plus proche de p sur le triangle (a, b, c), vectorisé (Ericson, « Real-Time
    Collision Detection », § 5.1.5). Tableaux (n, 3)."""
    ab, ac, ap = b - a, c - a, p - a
    d1 = np.sum(ab * ap, axis=1)
    d2 = np.sum(ac * ap, axis=1)
    bp = p - b
    d3 = np.sum(ab * bp, axis=1)
    d4 = np.sum(ac * bp, axis=1)
    cp = p - c
    d5 = np.sum(ab * cp, axis=1)
    d6 = np.sum(ac * cp, axis=1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    den = np.where(np.abs(va + vb + vc) > 0, va + vb + vc, 1e-300)
    v = vb / den
    w = vc / den
    out = a + ab * v[:, None] + ac * w[:, None]                  # intérieur
    # arêtes
    e_bc = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
    t = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) != 0, (d4 - d3) + (d5 - d6), 1e-300)
    out = np.where(e_bc[:, None], b + (c - b) * t[:, None], out)
    e_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    t = d2 / np.where(d2 - d6 != 0, d2 - d6, 1e-300)
    out = np.where(e_ac[:, None], a + ac * t[:, None], out)
    e_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    t = d1 / np.where(d1 - d3 != 0, d1 - d3, 1e-300)
    out = np.where(e_ab[:, None], a + ab * t[:, None], out)
    # sommets
    out = np.where(((d6 >= 0) & (d5 <= d6))[:, None], c, out)
    out = np.where(((d3 >= 0) & (d4 <= d3))[:, None], b, out)
    out = np.where(((d1 <= 0) & (d2 <= 0))[:, None], a, out)
    return out
