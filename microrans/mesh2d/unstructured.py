"""Maillages non structurés : triangles (DistMesh) et hybride couches limites + triangles.

- `triangulate` : algorithme DistMesh (Persson & Strang, SIAM Review 46, 2004) : les points
  sont des nœuds d'un réseau de ressorts dont la longueur au repos suit une fonction taille
  h(x) ; le domaine est défini par une distance signée (objets CSG). Triangles de bonne qualité,
  raffinement local par zones.
- `hybrid_mesh` : couches de quadrilatères extrudées depuis les parois (« inflation layers »
  de Fluent, « addLayers » de snappyHexMesh) pour résoudre la couche limite, puis triangles.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.spatial import Delaunay

from ..stop import check_stop
from .geometry import Polygon, Shape, signed_area
from .mesh import Mesh2D
from .ogrid import _normals


def size_function(h_max: float, refinements=()):
    """Fonction taille h(p) = min(h_max, min_i(h_i + g_i · distance_à_l_objet_i)).

    refinements : liste de dicts {"shape": Shape, "h": taille près de l'objet,
                  "growth": croissance (≈ 0.1-0.3)}.
    """
    refs = list(refinements)

    def h(p):
        out = np.full(len(p), float(h_max))
        for r in refs:
            d = np.abs(r["shape"].sdf(p))
            out = np.minimum(out, r["h"] + r.get("growth", 0.2) * d)
        return out
    return h


def _sample_curve_by_size(curve: np.ndarray, h, closed: bool = True) -> np.ndarray:
    """Échantillonne une polyligne avec un pas local ≈ h(x)."""
    pts = np.vstack([curve, curve[:1]]) if closed else curve
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    # suréchantillonnage fin puis intégration de ds/h
    m = max(int(s[-1] / (0.25 * np.min(h(pts)))) + 2, 8 * len(pts))
    u = np.linspace(0.0, s[-1], min(m, 200_000))
    fine = np.column_stack([np.interp(u, s, pts[:, 0]), np.interp(u, s, pts[:, 1])])
    w = 1.0 / h(fine)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(u))])
    n = max(int(np.round(cum[-1])), 3)
    targets = np.linspace(0.0, cum[-1], n, endpoint=not closed)
    uu = np.interp(targets, cum, u)
    return np.column_stack([np.interp(uu, s, pts[:, 0]), np.interp(uu, s, pts[:, 1])])


def _initial_points(domain: Shape, h, bbox, h_min, h_max, rng):
    """Points intérieurs : réseaux hexagonaux imbriqués (densité ∝ 1/h²)."""
    x0, y0, x1, y1 = bbox
    pts = []
    level_h = h_max
    while True:
        check_stop()                               # U20 : « Arrêter » de l'interface
        dx = level_h
        dy = dx * np.sqrt(3) / 2
        xs = np.arange(x0, x1 + dx, dx)
        ys = np.arange(y0, y1 + dy, dy)
        X, Y = np.meshgrid(xs, ys)
        X = X + (np.arange(len(ys))[:, None] % 2) * dx / 2
        P = np.column_stack([X.ravel(), Y.ravel()])
        hp = h(P)
        upper = np.inf if level_h >= h_max * 0.999 else 2 * level_h
        keep = (hp >= level_h) & (hp < upper) if level_h > h_min * 1.001 else (hp < upper)
        P = P[keep]
        if len(P):
            pts.append(P + rng.uniform(-0.05, 0.05, P.shape) * level_h)
        if level_h <= h_min * 1.001:
            break
        level_h = max(level_h / 2, h_min)
    P = np.vstack(pts) if pts else np.empty((0, 2))
    d = domain.sdf(P)
    return P[d < -0.4 * h(P)] if len(P) else P


def triangulate(domain: Shape, h_max: float, refinements=(), h=None, fixed_points=None,
                sample_leaves=None, max_iter: int = 300, dptol: float = 2e-3,
                seed: int = 0, verbose: bool = False, constrained=()) -> Mesh2D:
    """Maillage triangulaire d'un domaine CSG (DistMesh).

    domain : objet géométrique (ex. Rectangle(...) - Circle(...)).
    h_max  : taille maximale ; refinements : zones de raffinement (voir size_function).
    fixed_points : points imposés (ex. sommets d'une couche limite).
    sample_leaves : objets dont la frontière est échantillonnée (défaut : tous).
    constrained : objets dont la frontière ne doit contenir QUE des points fixes (raccord
        conforme avec des couches limites).
    """
    rng = np.random.default_rng(seed)
    hfun = h or size_function(h_max, refinements)
    bb = domain.bbox()
    pad = 1e-9 * max(bb[2] - bb[0], bb[3] - bb[1])
    # taille minimale effective
    probe = np.column_stack([rng.uniform(bb[0], bb[2], 20000), rng.uniform(bb[1], bb[3], 20000)])
    leaves = sample_leaves if sample_leaves is not None else domain.leaves()
    bpts = []
    for leaf in leaves:
        try:
            curve = leaf.boundary_curve(h=0.25 * min(hfun(probe).min(), h_max))
        except NotImplementedError:
            continue
        bpts.append(_sample_curve_by_size(curve, hfun))
    bpts = np.vstack(bpts) if bpts else np.empty((0, 2))
    if len(bpts):
        bpts = bpts[np.abs(domain.sdf(bpts)) < 1e-6 * (bb[2] - bb[0]) + pad]
    h_min = min(hfun(probe).min(), hfun(bpts).min() if len(bpts) else np.inf)
    pfix = domain.corners()
    if fixed_points is not None and len(fixed_points):
        pfix = np.vstack([pfix, np.asarray(fixed_points, float)]) if len(pfix) else \
            np.asarray(fixed_points, float)
    pint = _initial_points(domain, hfun, bb, h_min, h_max, rng)
    # retirer les points trop proches des points fixes / de frontière
    from scipy.spatial import cKDTree
    anchors = np.vstack([a for a in (pfix, bpts) if len(a)]) if (len(pfix) + len(bpts)) else None
    if anchors is not None and len(bpts):
        dist, _ = cKDTree(pfix).query(bpts) if len(pfix) else (np.full(len(bpts), np.inf), None)
        bpts = bpts[dist > 0.5 * hfun(bpts)]
        anchors = np.vstack([a for a in (pfix, bpts) if len(a)])
    if anchors is not None and len(pint):
        dist, _ = cKDTree(anchors).query(pint)
        pint = pint[dist > 0.6 * hfun(pint)]
    nfix = len(pfix)
    p = np.vstack([a for a in (pfix, bpts, pint) if len(a)])

    geps = 1e-3 * h_min
    deps = np.sqrt(np.finfo(float).eps) * h_min
    Fscale, deltat, ttol = 1.2, 0.2, 0.1
    pold = np.full_like(p, np.inf)
    bars = None
    for it in range(max_iter):
        check_stop()                               # U20 : avant, maillage mené à son terme
        hp = hfun(p)
        # re-triangulation si les points ont notablement bougé (centile : quelques nœuds
        # oscillant près d'une frontière ne doivent pas tout déclencher)
        if bars is None or np.percentile(np.linalg.norm(p - pold, axis=1) / hp, 99) > ttol:
            pold = p.copy()
            tri = Delaunay(p).simplices
            tri = tri[domain.sdf(p[tri].mean(axis=1)) < -geps]
            e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
            bars = np.unique(np.sort(e, axis=1), axis=0)
        bv = p[bars[:, 0]] - p[bars[:, 1]]
        L = np.linalg.norm(bv, axis=1)
        hb = 0.5 * (hp[bars[:, 0]] + hp[bars[:, 1]])
        L0 = hb * Fscale * np.sqrt(np.sum(L ** 2) / np.sum(hb ** 2))
        F = np.maximum(L0 - L, 0.0)
        Fv = (F / L)[:, None] * bv
        Ft = np.zeros_like(p)
        for k in range(2):
            Ft[:, k] += np.bincount(bars[:, 0], Fv[:, k], len(p))
            Ft[:, k] -= np.bincount(bars[:, 1], Fv[:, k], len(p))
        Ft[:nfix] = 0.0
        p = p + deltat * Ft
        d = domain.sdf(p)
        out = d > 0
        if np.any(out):
            q = p[out]
            gx = (domain.sdf(q + [deps, 0]) - d[out]) / deps
            gy = (domain.sdf(q + [0, deps]) - d[out]) / deps
            g2 = np.maximum(gx ** 2 + gy ** 2, 1e-300)
            p[out] = q - (d[out] / g2)[:, None] * np.column_stack([gx, gy])
        inner = d < -geps
        move = np.percentile(deltat * np.linalg.norm(Ft[inner], axis=1) / hp[inner], 99) \
            if np.any(inner) else 0.0
        if verbose and it % 50 == 0:
            print(f"  distmesh it {it}: {len(p)} points, déplacement max {move:.2e}")
        if move < dptol:
            break
    if constrained:
        # retirer les nœuds non fixes reprojetés sur une frontière contrainte
        hp = hfun(p)
        drop = np.zeros(len(p), dtype=bool)
        for c in constrained:
            drop |= np.abs(c.sdf(p)) < 0.3 * hp
        drop[:nfix] = False
        p = p[~drop]
    tri = Delaunay(p).simplices
    tri = tri[domain.sdf(p[tri].mean(axis=1)) < -geps]
    # triangles dégénérés
    a, b, c = p[tri[:, 0]], p[tri[:, 1]], p[tri[:, 2]]
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))
    tri = tri[area > 1e-12 * h_min ** 2]
    used = np.unique(tri)
    remap = -np.ones(len(p), dtype=int)
    remap[used] = np.arange(len(used))
    p, tri = p[used], remap[tri]
    return _mesh_from_triangles(p, tri, domain, hfun)


def _mesh_from_triangles(p, tri, domain: Shape, hfun, extra_boundary=None) -> Mesh2D:
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    es = np.sort(e, axis=1)
    uniq, cnt = np.unique(es, axis=0, return_counts=True)
    bedges = uniq[cnt == 1]
    mid = 0.5 * (p[bedges[:, 0]] + p[bedges[:, 1]])
    names = domain.patch_of(mid)
    inner = domain.sdf(mid) < -0.25 * hfun(mid)
    if np.any(inner):
        warnings.warn(f"{int(inner.sum())} arête(s) frontière à l'intérieur du domaine "
                      "(trou dans la triangulation).", stacklevel=2)
        names[inner] = "defaultFaces"
    boundary = {}
    for nm in dict.fromkeys(names):
        boundary[nm] = bedges[names == nm]
    types = {}
    for leaf in domain.leaves():
        for nm in leaf.patch_names():
            types.setdefault(nm, getattr(leaf, "patch_type", "patch"))
    return Mesh2D(p, tri, boundary, types)


def triangle_quality(mesh: Mesh2D) -> np.ndarray:
    """Qualité q = 2 r_inscrit / r_circonscrit des triangles (1 = équilatéral)."""
    tris = mesh.cell_nodes[mesh.cell_nv == 3][:, :3]
    p = mesh.points
    a = np.linalg.norm(p[tris[:, 1]] - p[tris[:, 2]], axis=1)
    b = np.linalg.norm(p[tris[:, 2]] - p[tris[:, 0]], axis=1)
    c = np.linalg.norm(p[tris[:, 0]] - p[tris[:, 1]], axis=1)
    return (b + c - a) * (c + a - b) * (a + b - c) / (a * b * c)


def inflation_layers(body: Shape, n_layers: int, first_height: float, ratio: float = 1.2,
                     h_surface: float | None = None, n_surface: int | None = None,
                     normal_smoothing: int = 5):
    """Couches de quadrilatères extrudées depuis la frontière d'un corps.

    Retourne (points (n_s, n_layers+1, 2), épaisseur totale).
    """
    curve = body.boundary_curve(n=n_surface) if n_surface else body.boundary_curve(h=h_surface)
    if signed_area(curve) < 0:
        curve = curve[::-1]
    nrm = _normals(curve, normal_smoothing)
    t = first_height * (ratio ** np.arange(n_layers + 1) - 1) / (ratio - 1) if ratio != 1 else \
        first_height * np.arange(n_layers + 1)
    return curve[:, None] + t[None, :, None] * nrm[:, None], float(t[-1])


def hybrid_mesh(outer: Shape, bodies: list[Shape], h_max: float, h_surface: float,
                n_layers: int = 10, first_height: float = 1e-3, ratio: float = 1.2,
                growth: float = 0.2, n_surface: dict | None = None, refinements=(),
                max_iter: int = 300, verbose: bool = False) -> Mesh2D:
    """Maillage hybride : couches de quadrilatères autour des corps + triangles ailleurs.

    outer : domaine extérieur (ex. Rectangle avec noms inlet/outlet/...) ;
    bodies : corps (leurs noms deviennent des patches 'wall').
    """
    quad_meshes, holes, fixed = [], [], []
    refs = list(refinements)
    for body in bodies:
        ns = (n_surface or {}).get(body.name)
        X, thick = inflation_layers(body, n_layers, first_height, ratio, h_surface, ns)
        n = X.shape[0]
        idx = np.arange(n * (n_layers + 1)).reshape(n, n_layers + 1)
        ip = np.roll(np.arange(n), -1)
        cells = np.stack([idx[:, :-1], idx[ip, :-1], idx[ip, 1:], idx[:, 1:]], axis=-1).reshape(-1, 4)
        pts = X.reshape(-1, 2)
        q = pts[cells]
        area = 0.5 * np.sum(q[..., 0] * np.roll(q[..., 1], -1, axis=1)
                            - np.roll(q[..., 0], -1, axis=1) * q[..., 1], axis=1)
        if np.any(np.sign(area) != np.sign(np.median(area))):
            raise ValueError(f"Couches limites retournées autour de '{body.name}' : réduire "
                             "n_layers/first_height/ratio.")
        top = X[:, -1]
        quad_meshes.append(Mesh2D(pts, cells, {body.name: np.column_stack([idx[:, 0], idx[ip, 0]]),
                                               f"_{body.name}_top": np.column_stack([idx[:, -1], idx[ip, -1]])},
                                  {body.name: "wall"}))
        hole = Polygon(top, name=f"_{body.name}_top", sharp_angle=179.9)
        holes.append(hole)
        fixed.append(top)
        h_top = float(np.mean(np.linalg.norm(np.diff(np.vstack([top, top[:1]]), axis=0), axis=1)))
        refs.append({"shape": hole, "h": h_top, "growth": growth})
    domain = outer
    for hole in holes:
        domain = domain - hole
    tri_mesh = triangulate(domain, h_max, refinements=refs, fixed_points=np.vstack(fixed),
                           sample_leaves=outer.leaves(), max_iter=max_iter, verbose=verbose,
                           constrained=holes)
    merged = Mesh2D.merge(quad_meshes + [tri_mesh], patch_types={b.name: "wall" for b in bodies})
    left = [p for p in merged.patches if p.name.startswith("_") and p.size > 0]
    if left:
        warnings.warn("Raccord couches/triangles non conforme sur "
                      f"{sum(p.size for p in left)} face(s).", stacklevel=2)
    return merged
