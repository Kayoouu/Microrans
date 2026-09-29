"""Études numériques reproductibles : précision / coût des schémas en temps.

    microrans schemes -o docs        (≈ 2 min)

1D : canal turbulent pulsé (SA, Re_τ = 395) sur une période ; erreur relative de U(T)
     par rapport à une référence SDIRK3 à 1024 pas/période.
2D : tourbillon de Taylor-Green advecté (solution exacte), maillage 32×32, t = 2 ;
     l'erreur inclut l'erreur spatiale (plancher commun à tous les schémas).
"""
from __future__ import annotations

import time
import warnings

import numpy as np

SCHEME_COLORS = {       # palette catégorielle validée (ordre fixe), voir postprocess
    "euler": "#2a78d6", "bdf2": "#eb6834", "backward": "#eb6834", "cn": "#1baf7a",
    "crankNicolson": "#1baf7a", "sdirk2": "#eda100", "sdirk3": "#e87ba4", "rk2": "#eda100",
    "rk3": "#008300", "rk4": "#e87ba4", "ab2": "#4a3aa7",
}
SCHEME_MARKERS = {"euler": "o", "bdf2": "s", "backward": "s", "cn": "^", "crankNicolson": "^",
                  "sdirk2": "D", "sdirk3": "v", "rk2": "D", "rk3": "P", "rk4": "v", "ab2": "X"}


def time_study_1d(model: str = "sa", re_tau: float = 395.0, n_cells: int = 128,
                  steps=(8, 16, 32, 64, 128),
                  schemes=("euler", "bdf2", "cn", "sdirk2", "sdirk3"), explicit="rk3"):
    from .cases import PulsatingForcing
    from .grid import channel_grid
    from .models import get_model
    from .solver import explicit_dt_limit, model_diffusivities, solve_steady, solve_unsteady
    nu = 1.0 / re_tau
    g = channel_grid(n_cells, re_tau, 0.3)
    m = get_model(model, g, nu)
    st = solve_steady(m, g, nu, forcing=1.0)
    omega = 0.01 * re_tau
    T = 2 * np.pi / omega
    f = PulsatingForcing(1.0, 10.0, omega)

    def run(scheme, n):
        t0 = time.perf_counter()
        r = solve_unsteady(m, g, nu, f, st, t_end=T, dt=T / n, scheme=scheme, inner_tol=1e-8,
                           max_inner=50)
        return r.final.U, time.perf_counter() - t0, int(r.inner_iterations.sum())

    ref = run("sdirk3", 1024)[0]
    scale = np.max(np.abs(ref))
    rows = []
    for scheme in schemes:
        for n in steps:
            U, cpu, inner = run(scheme, n)
            rows.append(dict(scheme=scheme, steps=n, dt=T / n, error=float(
                np.max(np.abs(U - ref)) / scale), cpu=cpu, inner=inner))
    if explicit:
        lim = explicit_dt_limit(g, *model_diffusivities(m, g, nu, st.U, st.state),
                                scheme=explicit)
        n = int(np.ceil(T / lim))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            U, cpu, inner = run(explicit, n)
        rows.append(dict(scheme=explicit, steps=n, dt=T / n, error=float(
            np.max(np.abs(U - ref)) / scale), cpu=cpu, inner=inner))
    return rows


def time_study_2d(n: int = 32, dts=(0.16, 0.08, 0.04, 0.02),
                  schemes=("euler", "backward", "crankNicolson", "rk2", "rk3", "rk4", "ab2")):
    from .fv2d import Settings, Solver2D
    from .mesh2d import rectangle_mesh
    L = 2 * np.pi
    U0 = np.array([1.0, 0.5])
    nu, T = 0.01, 2.0

    def exact(C, t):
        F = np.exp(-2 * nu * t)
        x, y = C[:, 0] - U0[0] * t, C[:, 1] - U0[1] * t
        return np.column_stack([U0[0] - np.cos(x) * np.sin(y) * F,
                                U0[1] + np.sin(x) * np.cos(y) * F])

    m = rectangle_mesh(0, L, 0, L, n, n, names={"left": "L", "right": "R", "bottom": "B",
                                                "top": "T"}, periodic=[("L", "R"), ("B", "T")])
    rows = []
    for scheme in schemes:
        for dt in dts:
            s = Solver2D(m, nu, {}, settings=Settings(time_scheme=scheme),
                         reference_velocity=1.0)
            s.U = exact(m.cell_centers, 0.0)
            s.F_i = np.sum(s.fvm.interp(s.U) * s.fvm.Si, axis=1)
            co = s.courant(dt)[0]
            t0 = time.perf_counter()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    s.run_transient(dt, T)
                    err = float(np.max(np.abs(s.U - exact(m.cell_centers, T))))
                except FloatingPointError:
                    err = float("inf")
            if not np.isfinite(err) or err > 1.0:
                err = float("inf")                       # instable
            rows.append(dict(scheme=scheme, dt=dt, courant=co, error=err,
                             cpu=time.perf_counter() - t0))
    return rows


def plot_time_study(rows1d, rows2d, path):
    from .postprocess import TEXT_2, _pyplot
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    panels = [(axes[0], rows1d, "1D URANS (SA, Re_τ = 395, canal pulsé)",
               "erreur relative sur U(T)"),
              (axes[1], rows2d, "2D, tourbillon de Taylor-Green advecté (32×32)",
               "erreur max sur U (espace + temps)")]
    for ax, rows, title, ylab in panels:
        names = list(dict.fromkeys(r["scheme"] for r in rows))
        for scheme in names:
            pts = [(r["cpu"], r["error"]) for r in rows if r["scheme"] == scheme
                   and np.isfinite(r["error"])]
            if not pts:
                continue
            x, y = zip(*sorted(pts))
            c = SCHEME_COLORS.get(scheme, "#4a3aa7")
            ax.loglog(x, y, color=c, lw=2, marker=SCHEME_MARKERS.get(scheme, "o"), ms=6,
                      label=scheme)
            if len(names) <= 4:                  # étiquettes directes si peu de séries
                ax.annotate(scheme, (x[-1], y[-1]), textcoords="offset points",
                            xytext=(6, 0), va="center", fontsize=9, color=TEXT_2)
        ax.set(xlabel="temps CPU (s)", ylabel=ylab, title=title)
        ax.legend(fontsize=8, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def markdown_table(rows, cols):
    head = "| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n"
    body = ""
    for r in rows:
        cells = []
        for c in cols:
            v = r[c]
            cells.append(f"{v:.2e}" if isinstance(v, float) and c in ("error", "dt")
                         else f"{v:.2f}" if isinstance(v, float) else str(v))
        body += "| " + " | ".join(cells) + " |\n"
    return head + body
