"""Échantillonnage des champs : sondes (points de suivi), profils le long d'une ligne,
moyennes temporelles (moyenne et écart-type des fluctuations) en instationnaire.

Valeur en un point : cellule qui le contient, puis reconstruction linéaire
φ_P + ∇φ_P·(x − x_P) avec le gradient de Green-Gauss du solveur (ordre 2, cohérent avec
les volumes finis ; valeurs de paroi respectées à l'ordre 2). Points hors du domaine (dans
un corps, hors maillage) : NaN.
"""
from __future__ import annotations

import numpy as np


def parse_points(spec) -> np.ndarray:
    """[[x, y], ...] ou « x y ; x y » / « x, y ; x, y »."""
    if spec is None:
        return np.zeros((0, 2))
    if isinstance(spec, str):
        pts = [[float(v) for v in p.replace(",", " ").split()] for p in spec.split(";")
               if p.strip()]
    else:
        pts = [[float(v) for v in p] for p in spec]
    if any(len(p) != 2 for p in pts):
        raise ValueError(f"Points invalides : {spec!r} (attendu : x y ; x y …)")
    return np.asarray(pts, float).reshape(-1, 2)


def locate(mesh, pts, k: int = 12) -> np.ndarray:
    """Indice de la cellule contenant chaque point (−1 si aucune)."""
    from scipy.spatial import cKDTree
    pts = np.asarray(pts, float).reshape(-1, 2)
    if not len(pts):
        return np.zeros(0, dtype=int)
    k = min(k, mesh.n_cells)
    _, cand = cKDTree(mesh.cell_centers).query(pts, k=k)
    cand = cand.reshape(len(pts), k)
    out = -np.ones(len(pts), dtype=int)
    P = mesh.points
    cn, nv = mesh.cell_nodes, mesh.cell_nv
    for j in range(k):
        todo = out < 0
        if not todo.any():
            break
        c = cand[todo, j]
        inside = _inside(P, cn[c], nv[c], pts[todo])
        idx = np.nonzero(todo)[0][inside]
        out[idx] = c[inside]
    return out


def _inside(P, cn, nv, x):
    """Point dans le polygone (orienté trigonométrique, convexe ou non) : nombre
    d'enroulement non nul, tolérance aux arêtes."""
    m = cn.shape[1]
    idx = np.arange(m)[None, :]
    valid = idx < nv[:, None]
    nxt = np.where(idx + 1 < nv[:, None], idx + 1, 0)
    rows = np.arange(len(cn))[:, None]
    a = P[np.where(valid, cn, 0)]
    b = P[np.where(valid, cn[rows, nxt], 0)]
    ax, ay = a[..., 0] - x[:, None, 0], a[..., 1] - x[:, None, 1]
    bx, by = b[..., 0] - x[:, None, 0], b[..., 1] - x[:, None, 1]
    cross = ax * by - ay * bx
    scale = np.hypot(b[..., 0] - a[..., 0], b[..., 1] - a[..., 1])
    up = (ay <= 0) & (by > 0) & (cross > 0)
    down = (ay > 0) & (by <= 0) & (cross < 0)
    wind = np.sum(np.where(valid, up.astype(int) - down.astype(int), 0), axis=1)
    on_edge = np.any(valid & (np.abs(cross) <= 1e-12 * scale ** 2 + 1e-300)
                     & (ax * bx + ay * by <= 0), axis=1)
    return (wind != 0) | on_edge


class Sampler:
    """Valeurs des champs du solveur en des points fixes (cellule + gradient)."""

    def __init__(self, solver, pts):
        self.solver = solver
        self.pts = np.asarray(pts, float).reshape(-1, 2)
        self.cell = locate(solver.mesh, self.pts)
        self.ok = self.cell >= 0
        c = np.where(self.ok, self.cell, 0)
        self._c = c
        self._d = self.pts - solver.mesh.cell_centers[c]

    def _rec(self, phi, grad):
        h = self.solver.backend.to_host
        phi, grad = np.asarray(h(phi)), np.asarray(h(grad))
        v = phi[self._c] + np.sum(grad[self._c] * self._d, axis=1)
        return np.where(self.ok, v, np.nan)

    def sample(self, fields=None) -> dict:
        """{nom: valeurs aux points} : Ux, Uy, U_mag, p, variables de turbulence, T,
        moyennes temporelles si elles existent (fields : sous-ensemble à calculer)."""
        s = self.solver
        fvm = s.fvm
        want = (lambda k: True) if fields is None else (lambda k: k in fields)
        out = {}
        if want("Ux") or want("Uy") or want("U_mag"):
            gU = s.grad_U(s.U)
            out["Ux"] = self._rec(s.U[:, 0], gU[:, 0, :])
            out["Uy"] = self._rec(s.U[:, 1], gU[:, 1, :])
            out["U_mag"] = np.hypot(out["Ux"], out["Uy"])
        if want("p"):
            out["p"] = self._rec(s.p, fvm.grad(s.p, s.boundary_p(s.p)))
        for k, v in s.state.items():
            if want(k):
                a, b, _, _ = s.scalar_bc(k)
                out[k] = self._rec(v, fvm.grad(v, a * v[fvm.Pb] + b))
        if s.energy is not None and want("T"):
            out["T"] = self._rec(s.T, fvm.grad(s.T, s.boundary_T(s.T)))
        for k, v in s.scalars.items():
            if want(k):
                out[k] = self._rec(v, fvm.grad(v, s.boundary_scalar(k)))
        if s.rheology is not None and want("viscosity"):
            vv = np.asarray(s.backend.to_host(s.nu_lam))
            out["viscosity"] = np.where(self.ok, vv[self._c], np.nan)
        for k, v in s.mean_fields().items():
            if want(k):
                # moyennes : valeur de la cellule (pas de gradient stocké)
                vv = np.asarray(s.backend.to_host(v))
                out[k] = np.where(self.ok, vv[self._c], np.nan)
        if fields is not None:
            out = {k: out[k] for k in fields if k in out}
        return out


def line_points(start, end, n: int = 200):
    start, end = np.asarray(start, float), np.asarray(end, float)
    t = np.linspace(0.0, 1.0, int(n))
    pts = start[None, :] + t[:, None] * (end - start)[None, :]
    return t * np.linalg.norm(end - start), pts


def write_lines(solver, lines, out_dir, plot=True):
    """[[output.lines]] : name, start, end, n (défaut 200) → line_<name>.csv (+ .png)."""
    from pathlib import Path
    out_dir = Path(out_dir)
    written = []
    for i, ln in enumerate(lines or []):
        import re
        name = str(ln.get("name", f"ligne{i + 1}"))
        if not re.fullmatch(r"[\w-]{1,64}", name):
            raise ValueError(f"[[output.lines]] name = {name!r} : lettres, chiffres, _ ou - "
                             f"seulement (nom de fichier).")
        s, pts = line_points(ln["start"], ln["end"], ln.get("n", 200))
        vals = Sampler(solver, pts).sample()
        cols = {"s": s, "x": pts[:, 0], "y": pts[:, 1], **vals}
        path = out_dir / f"line_{name}.csv"
        np.savetxt(path, np.column_stack(list(cols.values())), delimiter=",",
                   header=",".join(cols), comments="", encoding="utf-8")
        written.append(path)
        if plot:
            _plot_line(cols, out_dir / f"line_{name}.png", name)
    return written


def _plot_line(cols, path, name):
    from matplotlib.figure import Figure
    keys = [k for k in cols if k not in ("s", "x", "y", "U_mag") and "_mean" not in k
            and "_rms" not in k and k not in ("k", "omega", "epsilon", "nu_tilde")][:6]
    fig = Figure(figsize=(4.0 * len(keys), 3.4), layout="constrained")
    axes = np.atleast_1d(fig.subplots(1, len(keys)))
    for ax, k in zip(axes, keys):
        ax.plot(cols["s"], cols[k], color="#2a78d6", lw=1.6)
        ax.set(xlabel="abscisse s le long de la ligne", ylabel=k)
        ax.grid(True, alpha=0.3)
    fig.suptitle(f"Profil « {name} »")
    fig.savefig(path, dpi=120)


class TimeAverage:
    """Moyenne temporelle pondérée par Δt et écart-type des fluctuations (RMS) de U, p (, T),
    à partir de t_start — comme fieldAverage d'OpenFOAM."""


    def __init__(self, solver, t_start: float):
        self.solver, self.t_start = solver, float(t_start)
        self.weight = 0.0
        self.sum, self.sumsq = {}, {}

    def _values(self):
        s = self.solver
        v = {"Ux": s.U[:, 0], "Uy": s.U[:, 1], "p": s.p}
        if s.energy is not None:
            v["T"] = s.T
        v.update(s.scalars)
        return v

    def update(self):
        s = self.solver
        if s.time <= self.t_start + 1e-12:
            return
        dt = min(s.dt, s.time - self.t_start)
        for k, v in self._values().items():
            self.sum[k] = self.sum.get(k, 0.0) + dt * v
            self.sumsq[k] = self.sumsq.get(k, 0.0) + dt * v * v
        self.weight += dt

    def fields(self) -> dict:
        if self.weight <= 0.0:
            return {}
        xp = self.solver.xp
        out = {}
        for k in self.sum:
            m = self.sum[k] / self.weight
            out[f"{k}_mean"] = m
            out[f"{k}_rms"] = xp.sqrt(xp.maximum(self.sumsq[k] / self.weight - m * m, 0.0))
        return out

    # sauvegarde / reprise
    def state(self) -> dict:
        h = self.solver.backend.to_host
        d = {"avg_t_start": np.array(self.t_start), "avg_weight": np.array(self.weight)}
        for k in self.sum:
            d[f"avg_sum_{k}"] = h(self.sum[k])
            d[f"avg_sumsq_{k}"] = h(self.sumsq[k])
        return d

    def load(self, d: dict):
        A = self.solver.backend.asarray
        self.weight = float(d["avg_weight"])
        for k in self._values():
            if f"avg_sum_{k}" in d:
                self.sum[k] = A(d[f"avg_sum_{k}"].copy())
                self.sumsq[k] = A(d[f"avg_sumsq_{k}"].copy())
