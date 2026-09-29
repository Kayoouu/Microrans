"""Géométrie 2D : primitives, opérations booléennes (CSG), transformations.

Chaque objet définit une fonction distance signée `sdf(p)` (négative à l'intérieur), ce qui
permet de combiner les objets (union `|`, différence `-`, intersection `&`) comme dans les
modeleurs CSG, et une courbe frontière fermée (sens trigonométrique) pour les mailleurs
structurés. Le nom d'un objet devient le nom du « patch » (condition aux limites) de sa frontière.
"""
from __future__ import annotations

import numpy as np

_CHUNK = 4096


def _as_points(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    return p.reshape(-1, 2)


def resample_closed(points: np.ndarray, n: int | None = None, h: float | None = None,
                    keep_corners: bool = True, corner_angle: float = 30.0) -> np.ndarray:
    """Rééchantillonne une polyligne fermée à pas d'abscisse curviligne ~constant.

    Les sommets anguleux (déviation > corner_angle degrés) sont conservés.
    """
    pts = _as_points(points)
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    closed = np.vstack([pts, pts[:1]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    if n is None:
        if h is None:
            raise ValueError("Il faut n ou h.")
        n = max(int(np.ceil(total / h)), 3)
    corners = []
    if keep_corners:
        t_prev = pts - np.roll(pts, 1, axis=0)
        t_next = np.roll(pts, -1, axis=0) - pts
        cosang = np.sum(t_prev * t_next, axis=1) / (
            np.linalg.norm(t_prev, axis=1) * np.linalg.norm(t_next, axis=1) + 1e-300)
        corners = list(np.nonzero(cosang < np.cos(np.radians(corner_angle)))[0])
    if not corners:
        u = np.linspace(0.0, total, n, endpoint=False)
        return np.column_stack([np.interp(u, s, closed[:, 0]), np.interp(u, s, closed[:, 1])])
    # Répartition par tronçons entre coins, proportionnelle à la longueur
    out = []
    sc = s[corners]
    sc_next = np.append(sc[1:], sc[0] + total)
    for a, b in zip(sc, sc_next):
        m = max(int(round(n * (b - a) / total)), 1)
        u = np.mod(np.linspace(a, b, m, endpoint=False), total)
        out.append(np.column_stack([np.interp(u, s, closed[:, 0]), np.interp(u, s, closed[:, 1])]))
    return np.vstack(out)


def signed_area(points: np.ndarray) -> float:
    x, y = points[:, 0], points[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


class Shape:
    """Objet géométrique 2D (classe de base)."""

    name: str = "body"
    patch_type: str = "patch"     # type du patch créé sur la frontière (wall pour un corps)

    # -- interface ---------------------------------------------------------------
    def sdf(self, p) -> np.ndarray:
        raise NotImplementedError

    def bbox(self) -> tuple[float, float, float, float]:
        raise NotImplementedError

    def boundary_curve(self, n: int | None = None, h: float | None = None) -> np.ndarray:
        """Polyligne fermée (sens trigonométrique, sans point dupliqué)."""
        raise NotImplementedError(f"{type(self).__name__} n'a pas de frontière paramétrée unique.")

    def corners(self) -> np.ndarray:
        """Points anguleux à conserver comme points fixes lors du maillage."""
        return np.empty((0, 2))

    def leaves(self) -> list["Shape"]:
        return [self]

    def patch_of(self, p) -> np.ndarray:
        """Nom de patch pour des points situés sur la frontière de cet objet."""
        return np.full(len(_as_points(p)), self.name, dtype=object)

    def patch_names(self) -> list[str]:
        return [self.name]

    # -- CSG -----------------------------------------------------------------------
    def __or__(self, other):
        return Union(self, other)

    def __sub__(self, other):
        return Difference(self, other)

    def __and__(self, other):
        return Intersection(self, other)

    # -- transformations --------------------------------------------------------------
    def translate(self, dx: float, dy: float) -> "Shape":
        return Transformed(self, offset=(dx, dy))

    def rotate(self, angle_deg: float, center=(0.0, 0.0)) -> "Shape":
        return Transformed(self, angle=angle_deg, center=center)

    def scale(self, factor: float, center=(0.0, 0.0)) -> "Shape":
        return Transformed(self, factor=factor, center=center)

    def renamed(self, name: str) -> "Shape":
        self.name = name
        return self

    def as_wall(self) -> "Shape":
        """Marque la frontière de l'objet comme paroi (patch de type wall)."""
        self.patch_type = "wall"
        for leaf in self.leaves():
            if leaf is not self:
                leaf.as_wall()
        return self


class Circle(Shape):
    def __init__(self, center=(0.0, 0.0), radius: float = 0.5, name: str = "circle"):
        self.center = np.asarray(center, dtype=float)
        self.radius = float(radius)
        self.name = name

    def sdf(self, p):
        return np.linalg.norm(_as_points(p) - self.center, axis=1) - self.radius

    def bbox(self):
        c, r = self.center, self.radius
        return (c[0] - r, c[1] - r, c[0] + r, c[1] + r)

    def boundary_curve(self, n=None, h=None):
        if n is None:
            n = max(int(np.ceil(2 * np.pi * self.radius / h)), 8)
        t = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
        return self.center + self.radius * np.column_stack([np.cos(t), np.sin(t)])


class Rectangle(Shape):
    """Rectangle ; chaque côté peut porter son propre nom de patch."""

    SIDES = ("bottom", "right", "top", "left")

    def __init__(self, x0, y0, x1, y1, name: str = "rectangle", names: dict | None = None):
        self.x0, self.y0 = float(min(x0, x1)), float(min(y0, y1))
        self.x1, self.y1 = float(max(x0, x1)), float(max(y0, y1))
        self.name = name
        self.names = {s: (names or {}).get(s, name) for s in self.SIDES}

    def sdf(self, p):
        p = _as_points(p)
        c = np.array([0.5 * (self.x0 + self.x1), 0.5 * (self.y0 + self.y1)])
        half = np.array([0.5 * (self.x1 - self.x0), 0.5 * (self.y1 - self.y0)])
        q = np.abs(p - c) - half
        outside = np.linalg.norm(np.maximum(q, 0.0), axis=1)
        inside = np.minimum(np.max(q, axis=1), 0.0)
        return outside + inside

    def bbox(self):
        return (self.x0, self.y0, self.x1, self.y1)

    def corners(self):
        return np.array([[self.x0, self.y0], [self.x1, self.y0], [self.x1, self.y1],
                         [self.x0, self.y1]])

    def boundary_curve(self, n=None, h=None):
        pts = self.corners()
        if h is None:
            h = 2 * ((self.x1 - self.x0) + (self.y1 - self.y0)) / n
        return resample_closed(pts, h=h)

    def patch_of(self, p):
        p = _as_points(p)
        dist = np.column_stack([np.abs(p[:, 1] - self.y0), np.abs(p[:, 0] - self.x1),
                                np.abs(p[:, 1] - self.y1), np.abs(p[:, 0] - self.x0)])
        side = np.argmin(dist, axis=1)
        names = np.array([self.names[s] for s in self.SIDES], dtype=object)
        return names[side]

    def patch_names(self):
        return list(dict.fromkeys(self.names[s] for s in self.SIDES))


class Polygon(Shape):
    """Polygone fermé quelconque (distance signée exacte, signe par parité des croisements)."""

    def __init__(self, points, name: str = "polygon", sharp_angle: float = 30.0):
        pts = _as_points(points)
        if np.allclose(pts[0], pts[-1]):
            pts = pts[:-1]
        if len(pts) < 3:
            raise ValueError("Un polygone demande au moins 3 points.")
        if signed_area(pts) < 0:
            pts = pts[::-1]
        self.points = pts
        self.name = name
        self.sharp_angle = sharp_angle

    def sdf(self, p):
        p = _as_points(p)
        if len(self.points) > 48 and len(p) > 64:
            return self._sdf_fast(p)
        return self._sdf_brute(p)

    def _prepare_fast(self):
        """Arêtes subdivisées (longueur bornée) + KD-tree des milieux, chemin matplotlib."""
        from matplotlib.path import Path as MplPath
        from scipy.spatial import cKDTree
        a = self.points
        b = np.roll(a, -1, axis=0)
        L = np.linalg.norm(b - a, axis=1)
        lmax = max(np.median(L), 1e-300)
        k = np.maximum(np.ceil(L / lmax).astype(int), 1)
        t0 = np.concatenate([np.arange(m) / m for m in k])
        t1 = np.concatenate([np.arange(1, m + 1) / m for m in k])
        seg = np.repeat(np.arange(len(a)), k)
        A = a[seg] + t0[:, None] * (b - a)[seg]
        B = a[seg] + t1[:, None] * (b - a)[seg]
        self._fast = (A, B, cKDTree(0.5 * (A + B)), MplPath(np.vstack([a, a[:1]]), closed=True))

    def _sdf_fast(self, p):
        if getattr(self, "_fast", None) is None:
            self._prepare_fast()
        A, B, tree, path = self._fast
        kk = min(8, len(A))
        _, j = tree.query(p, k=kk)
        a, b = A[j], B[j]
        ab = b - a
        ap = p[:, None, :] - a
        t = np.clip(np.sum(ap * ab, axis=2) / np.maximum(np.sum(ab * ab, axis=2), 1e-300), 0.0, 1.0)
        dist = np.sqrt(np.min(np.sum((ap - t[..., None] * ab) ** 2, axis=2), axis=1))
        inside = path.contains_points(p)
        return np.where(inside, -dist, dist)

    def _sdf_brute(self, p):
        a = self.points
        b = np.roll(a, -1, axis=0)
        ab = b - a
        ab2 = np.maximum(np.sum(ab * ab, axis=1), 1e-300)
        out = np.empty(len(p))
        for s in range(0, len(p), _CHUNK):
            q = p[s:s + _CHUNK, None, :]
            ap = q - a[None]
            t = np.clip(np.sum(ap * ab[None], axis=2) / ab2[None], 0.0, 1.0)
            d2 = np.sum((ap - t[..., None] * ab[None]) ** 2, axis=2)
            dist = np.sqrt(np.min(d2, axis=1))
            # parité des croisements d'une demi-droite horizontale
            y = q[..., 1]
            cond = (a[None, :, 1] > y) != (b[None, :, 1] > y)
            with np.errstate(divide="ignore", invalid="ignore"):
                xint = a[None, :, 0] + (y - a[None, :, 1]) * ab[None, :, 0] / ab[None, :, 1]
            inside = np.sum(cond & (q[..., 0] < xint), axis=1) % 2 == 1
            out[s:s + _CHUNK] = np.where(inside, -dist, dist)
        return out

    def bbox(self):
        mn, mx = self.points.min(axis=0), self.points.max(axis=0)
        return (mn[0], mn[1], mx[0], mx[1])

    def corners(self):
        p = self.points
        t_prev = p - np.roll(p, 1, axis=0)
        t_next = np.roll(p, -1, axis=0) - p
        cosang = np.sum(t_prev * t_next, axis=1) / (
            np.linalg.norm(t_prev, axis=1) * np.linalg.norm(t_next, axis=1) + 1e-300)
        return p[cosang < np.cos(np.radians(self.sharp_angle))]

    def boundary_curve(self, n=None, h=None):
        if n is None and h is None:
            return self.points.copy()
        return resample_closed(self.points, n=n, h=h, corner_angle=self.sharp_angle)


class Ellipse(Polygon):
    def __init__(self, center=(0.0, 0.0), a: float = 1.0, b: float = 0.5, n: int = 256,
                 name: str = "ellipse"):
        t = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
        pts = np.asarray(center, float) + np.column_stack([a * np.cos(t), b * np.sin(t)])
        super().__init__(pts, name=name)


def naca4_points(code: str = "0012", chord: float = 1.0, n: int = 201,
                 closed_te: bool = True) -> np.ndarray:
    """Profil NACA 4 chiffres (répartition en cosinus, bord d'attaque raffiné).

    Retourne un contour fermé : extrados du bord de fuite vers le bord d'attaque, puis intrados.
    """
    code = str(code).strip()
    if len(code) != 4 or not code.isdigit():
        raise ValueError("Code NACA 4 chiffres attendu, ex. '0012' ou '2412'.")
    m, p, t = int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
    npts = max(n // 2, 20)
    beta = np.linspace(0.0, np.pi, npts)
    x = 0.5 * (1.0 - np.cos(beta))
    a4 = -0.1036 if closed_te else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 + a4 * x ** 4)
    if m > 0 and p > 0:
        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2),
                      m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
        dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    upper = np.column_stack([xu, yu])[::-1]          # BF -> BA
    lower = np.column_stack([xl, yl])[1:]            # BA -> BF (sans doubler le BA)
    pts = np.vstack([upper, lower]) * chord
    if closed_te:
        pts = pts[:-1]                               # le BF de l'intrados = BF de l'extrados
    return pts


class NACA4(Polygon):
    def __init__(self, code: str = "0012", chord: float = 1.0, n: int = 201,
                 closed_te: bool = True, name: str = "airfoil"):
        self.code, self.chord, self.closed_te = code, chord, closed_te
        super().__init__(naca4_points(code, chord, n, closed_te), name=name, sharp_angle=60.0)

    def boundary_curve(self, n=None, h=None):
        # Répartition native en cosinus (raffinée aux bords d'attaque et de fuite)
        if n is not None:
            return Polygon(naca4_points(self.code, self.chord, n + 2, self.closed_te)).points
        return super().boundary_curve(n, h)


class Spline(Polygon):
    """Courbe fermée lisse (spline cubique périodique) passant par des points de contrôle."""

    def __init__(self, control_points, n: int = 256, name: str = "spline"):
        from scipy.interpolate import CubicSpline
        c = _as_points(control_points)
        if np.allclose(c[0], c[-1]):
            c = c[:-1]
        closed = np.vstack([c, c[:1]])
        s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))])
        cs = CubicSpline(s, closed, bc_type="periodic")
        pts = cs(np.linspace(0.0, s[-1], n, endpoint=False))
        super().__init__(pts, name=name, sharp_angle=179.0)


# ------------------------------------------------------------------------- CSG
class _Boolean(Shape):
    def __init__(self, a: Shape, b: Shape):
        self.a, self.b = a, b
        self.name = a.name

    def leaves(self):
        return self.a.leaves() + self.b.leaves()

    def patch_names(self):
        out = []
        for leaf in self.leaves():
            out += leaf.patch_names()
        return list(dict.fromkeys(out))

    def patch_of(self, p):
        p = _as_points(p)
        leaves = self.leaves()
        d = np.column_stack([np.abs(leaf.sdf(p)) for leaf in leaves])
        which = np.argmin(d, axis=1)
        out = np.empty(len(p), dtype=object)
        for i, leaf in enumerate(leaves):
            m = which == i
            if np.any(m):
                out[m] = leaf.patch_of(p[m])
        return out

    def corners(self):
        cands = [leaf.corners() for leaf in self.leaves()]
        cands = np.vstack(cands) if cands else np.empty((0, 2))
        if len(cands) == 0:
            return cands
        bb = self.bbox()
        tol = 1e-9 * max(bb[2] - bb[0], bb[3] - bb[1])
        return cands[np.abs(self.sdf(cands)) < tol]


class Union(_Boolean):
    def sdf(self, p):
        return np.minimum(self.a.sdf(p), self.b.sdf(p))

    def bbox(self):
        a, b = self.a.bbox(), self.b.bbox()
        return (min(a[0], b[0]), min(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3]))


class Difference(_Boolean):
    def sdf(self, p):
        return np.maximum(self.a.sdf(p), -self.b.sdf(p))

    def bbox(self):
        return self.a.bbox()


class Intersection(_Boolean):
    def sdf(self, p):
        return np.maximum(self.a.sdf(p), self.b.sdf(p))

    def bbox(self):
        a, b = self.a.bbox(), self.b.bbox()
        return (max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3]))


class Transformed(Shape):
    """Similitude (rotation, homothétie, translation) appliquée à un objet."""

    def __init__(self, shape: Shape, angle: float = 0.0, factor: float = 1.0,
                 offset=(0.0, 0.0), center=(0.0, 0.0)):
        self.shape = shape
        self.name = shape.name
        self.patch_type = shape.patch_type
        th = np.radians(angle)
        self.R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        self.factor = float(factor)
        self.center = np.asarray(center, float)
        self.offset = np.asarray(offset, float)

    def forward(self, p):
        p = _as_points(p)
        return (p - self.center) @ self.R.T * self.factor + self.center + self.offset

    def inverse(self, p):
        p = _as_points(p)
        return ((p - self.offset - self.center) / self.factor) @ self.R + self.center

    def sdf(self, p):
        return self.shape.sdf(self.inverse(p)) * self.factor

    def bbox(self):
        try:
            pts = self.forward(self.shape.boundary_curve(n=720))
        except NotImplementedError:
            b = self.shape.bbox()
            pts = self.forward(np.array([[b[0], b[1]], [b[2], b[1]], [b[2], b[3]], [b[0], b[3]]]))
        return (pts[:, 0].min(), pts[:, 1].min(), pts[:, 0].max(), pts[:, 1].max())

    def boundary_curve(self, n=None, h=None):
        return self.forward(self.shape.boundary_curve(n=n, h=None if h is None else h / self.factor))

    def corners(self):
        c = self.shape.corners()
        return self.forward(c) if len(c) else c

    def leaves(self):
        if isinstance(self.shape, _Boolean):
            return [_TransformedLeaf(leaf, self) for leaf in self.shape.leaves()]
        return [self]

    def patch_of(self, p):
        return self.shape.patch_of(self.inverse(p))

    def patch_names(self):
        return self.shape.patch_names()

    def renamed(self, name):
        self.name = name
        self.shape.renamed(name)
        return self


class _TransformedLeaf(Shape):
    def __init__(self, leaf: Shape, tr: Transformed):
        self.leaf, self.tr = leaf, tr
        self.name = leaf.name
        self.patch_type = leaf.patch_type

    def sdf(self, p):
        return self.leaf.sdf(self.tr.inverse(p)) * self.tr.factor

    def patch_of(self, p):
        return self.leaf.patch_of(self.tr.inverse(p))

    def corners(self):
        c = self.leaf.corners()
        return self.tr.forward(c) if len(c) else c

    def patch_names(self):
        return self.leaf.patch_names()


# ------------------------------------------------------------ construction depuis un dict
def shape_from_dict(spec: dict, base_dir=".") -> Shape:
    """Construit un objet depuis une description (fichier de configuration TOML/JSON).

    Types : circle, rectangle, ellipse, polygon, naca, spline, file.
    Transformations optionnelles (autour de `rotation_center`) :
      - angle     : rotation en degrés, sens trigonométrique ;
      - incidence : angle d'attaque en degrés (convention aéro : nez vers le haut = horaire) ;
      - scale, translate.
    """
    from .io import read_curve

    spec = dict(spec)
    kind = spec.pop("type").lower()
    name = spec.pop("name", kind)
    patch_type = spec.pop("patch_type", None)
    angle = spec.pop("angle", 0.0) - spec.pop("incidence", 0.0)
    rot_center = spec.pop("rotation_center", (0.0, 0.0))
    scale = spec.pop("scale", 1.0)
    translate = spec.pop("translate", (0.0, 0.0))
    if kind == "circle":
        s = Circle(spec.get("center", (0.0, 0.0)), spec["radius"], name=name)
    elif kind == "rectangle":
        s = Rectangle(spec["x0"], spec["y0"], spec["x1"], spec["y1"], name=name,
                      names=spec.get("names"))
    elif kind == "ellipse":
        s = Ellipse(spec.get("center", (0.0, 0.0)), spec["a"], spec["b"], name=name)
    elif kind == "polygon":
        s = Polygon(spec["points"], name=name)
    elif kind == "naca":
        s = NACA4(str(spec.get("code", "0012")), spec.get("chord", 1.0), spec.get("n", 201),
                  spec.get("closed_te", True), name=name)
    elif kind == "spline":
        s = Spline(spec["points"], spec.get("n", 256), name=name)
    elif kind == "file":
        from pathlib import Path
        path = Path(spec["path"])
        if not path.is_absolute():
            path = Path(base_dir) / path
        s = Polygon(read_curve(path), name=name,
                    sharp_angle=spec.get("sharp_angle", 60.0))
    else:
        raise ValueError(f"Type de géométrie inconnu : {kind}")
    if patch_type:
        s.patch_type = patch_type
    if angle or scale != 1.0 or any(translate):
        s = Transformed(s, angle=angle, factor=scale, offset=translate, center=rot_center)
    return s
