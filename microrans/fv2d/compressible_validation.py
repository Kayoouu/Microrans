"""Validation du solveur compressible : produit les tableaux et figures de
docs/compressible.md (valeurs MESURÉES, comparées aux références citées).

    python -m microrans.fv2d.compressible_validation [--out docs] [--quick] [--skip-naca]

1. Tube à choc de Sod : erreurs L1 contre la solution exacte (gasdynamics.riemann_exact),
   ordre observé, positions du choc, du contact et du milieu de la détente.
2. Rampe M = 2, θ = 10° : angle du choc, p2/p1, M2 contre la théorie du choc oblique.
3. Plaque plane laminaire M = 0.2, Re_L = 1e5 : C_f contre Blasius (2 maillages).
4. Couette compressible (dissipation visqueuse) : solution analytique.
5. NACA 0012 Euler M = 0.8, α = 1.25° : Cl, Cd sur 3 maillages en O (+ M = 0.5, α = 1.25°).
Chaque cas passe par les fichiers d'exemple (run_case) quand il en existe un.
"""
from __future__ import annotations

import argparse
import copy
import time
from pathlib import Path

import numpy as np

from ..postprocess import REF_COLOR, _pyplot
from .gasdynamics import blasius_cf, oblique_shock, riemann_exact

BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"


def _cfg(name):
    from ..cli import examples_dir
    from ..mesh2d.builder import load_config
    return load_config(examples_dir() / name), examples_dir()


def _run(cfg, base, out=None):
    from .case import run_case
    t0 = time.perf_counter()
    s, solver = run_case(cfg, base_dir=base, out_dir=out or "results/_validation",
                         verbose=False, plot=False, return_solver=True)
    return s, solver, time.perf_counter() - t0


def _cross(x, f, level, lo, hi, falling=True):
    """Abscisse (interpolée) où f croise `level` dans ]lo, hi[."""
    sel = (x > lo) & (x < hi)
    xs, fs = x[sel], f[sel]
    s = (fs[:-1] - level) * (fs[1:] - level) <= 0
    k = np.nonzero(s)[0]
    if not len(k):
        return float("nan")
    k = k[0]
    return float(xs[k] + (level - fs[k]) / (fs[k + 1] - fs[k]) * (xs[k + 1] - xs[k]))


# ------------------------------------------------------------------ 1. Sod
def sod(out: Path, ns=(100, 200, 400), fluxes=("roe", "hllc")):
    rows = []
    prof = None
    for flux in fluxes:
        for n in ns:
            cfg, base = _cfg("compressible_tube_sod.toml")
            cfg["mesh"].update(nx=n, y1=1.0 / n)
            cfg["solver"]["flux"] = flux
            cfg["output"] = {"vtk": False, "checkpoint": False}
            s, S, wt = _run(cfg, base)
            t = s["time"]
            a = np.sqrt(1e5 / 1.0)
            x = S.mesh.cell_centers[:, 0]
            r, u, p, w = riemann_exact((1.0, 0.0, 1e5), (0.125, 0.0, 1e4), (x - 0.5) / t)
            W = S.W
            xe = lambda sp: 0.5 + sp * t                                    # noqa: E731
            shock = _cross(x, W[:, 0], 0.5 * (w["rho_star_right"] + 0.125), 0.75, 0.95)
            contact = _cross(x, W[:, 0], 0.5 * (w["rho_star_left"] + w["rho_star_right"]),
                             0.6, 0.8)
            fan = _cross(x, W[:, 1], 0.5 * w["u_star"], 0.2, 0.55, falling=False)
            g = 1.4
            cl = np.sqrt(g * 1e5)
            xi_mid = (g + 1) / 2 * 0.5 * w["u_star"] - cl     # u = u*/2 dans la détente
            rows.append({"flux": flux, "n": n,
                         "L1_rho": float(np.mean(np.abs(W[:, 0] - r))),
                         "L1_u": float(np.mean(np.abs(W[:, 1] - u)) / a),
                         "L1_p": float(np.mean(np.abs(W[:, 3] - p)) / 1e5),
                         "x_shock": shock, "x_shock_exact": xe(w["right_shock"]),
                         "x_contact": contact, "x_contact_exact": xe(w["contact"]),
                         "x_fan_mid": fan, "x_fan_mid_exact": 0.5 + xi_mid * t,
                         "steps": s["steps"], "time_s": round(wt, 2),
                         "mass_drift": s["totals_final"]["mass"] / s["totals_initial"]["mass"]
                         - 1.0})
            if flux == "roe" and n == 200:
                prof = (x, W.copy(), r, u, p, a)
    for fl in fluxes:
        sub = [r for r in rows if r["flux"] == fl]
        for r0, r1 in zip(sub[:-1], sub[1:]):
            r1["order_rho"] = float(np.log(r0["L1_rho"] / r1["L1_rho"]) / np.log(
                r1["n"] / r0["n"]))
    if prof is not None:
        x, W, r, u, p, a = prof
        plt = _pyplot()
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
        xe = np.linspace(0, 1, 2001)
        re_, ue, pe, _ = riemann_exact((1.0, 0, 1e5), (0.125, 0, 1e4), (xe - 0.5) / (
            0.2 / a))
        for ax, num, ex, lab in zip(axs, (W[:, 0], W[:, 1] / a, W[:, 3] / 1e5),
                                    (re_, ue / a, pe / 1e5),
                                    ("ρ (kg/m³)", "u / √(p_L/ρ_L)", "p / p_L")):
            ax.plot(xe, ex, color=REF_COLOR, lw=1.0, label="solution exacte")
            ax.plot(x, num, "o", ms=2.5, color=BLUE, label="Roe + MUSCL, 200 mailles")
            ax.set(xlabel="x (m)", ylabel=lab)
            ax.grid(True, alpha=0.3)
        axs[0].legend(loc="upper right", fontsize=8)
        fig.suptitle("Tube à choc de Sod, t = 0.2 L/√(p_L/ρ_L)")
        fig.tight_layout()
        fig.savefig(out / "compressible_sod.png", dpi=130)
        plt.close(fig)
    return rows


# ------------------------------------------------------------------ 2. rampe
def ramp(out: Path, factors=(0.5, 1.0, 2.0)):
    ref = oblique_shock(2.0, 10.0)
    rows = []
    for fct in factors:
        cfg, base = _cfg("compressible_rampe_mach2.toml")
        for b in cfg["mesh"]["blocks"]:
            b["cells"] = [int(round(c * fct)) for c in b["cells"]]
        cfg["solver"].update(steady_scheme="implicit", cfl=5.0, max_iter=3000)
        cfg["output"] = {"vtk": False, "checkpoint": False, "forces": ["ramp"]}
        s, S, wt = _run(cfg, base)
        C = S.mesh.cell_centers
        pr = S.p / S.fs.p
        from scipy.interpolate import LinearNDInterpolator
        f = LinearNDInterpolator(C, pr)
        pm = 0.5 * (1 + ref["p_ratio"])
        ys = np.linspace(0.2, 0.7, 11)
        xs = []
        for y in ys:
            xx = np.linspace(0.55, 1.49, 3000)
            xs.append(_cross(xx, f(xx, np.full_like(xx, y)), pm, 0.55, 1.49))
        xs = np.array(xs)
        ok = np.isfinite(xs)
        beta = float(np.degrees(np.arctan(np.polyfit(xs[ok], ys[ok], 1)[0])))
        xi, yi = C[:, 0], C[:, 1]
        yr = (xi - 0.5) * np.tan(np.radians(10.0))
        yb = (xi - 0.5) * np.tan(np.radians(ref["beta_deg"]))
        zone = (xi > 0.9) & (xi < 1.45) & (yi > yr + 0.02) & (yi < yb - 0.08)
        rows.append({"cells": S.nc, "beta": beta, "beta_exact": ref["beta_deg"],
                     "p2_p1": float(pr[zone].mean()), "p2_p1_exact": ref["p_ratio"],
                     "M2": float(S.mach[zone].mean()), "M2_exact": ref["M2"],
                     "Cd_ramp": s["ramp"]["Cd"],
                     "Cd_ramp_exact": (ref["p_ratio"] - 1) * S.fs.p * np.tan(np.radians(10))
                     / s["dynamic_pressure"],
                     "iterations": s["iterations"], "converged": s["converged"],
                     "time_s": round(wt, 1)})
        if fct == 1.0:
            from ..mesh2d.plot import plot_field
            plot_field(S.mesh, S.mach, out / "compressible_rampe.png",
                       title="Rampe M = 2, θ = 10° : nombre de Mach (Euler, 5 400 cellules)",
                       cmap="viridis", label="M")
    return rows


# ------------------------------------------------------------------ 3. plaque plane
def plate(out: Path, factors=(1.0, 2.0)):
    rows, curves = [], []
    for fct in factors:
        cfg, base = _cfg("compressible_plaque_laminaire.toml")
        for b in cfg["mesh"]["blocks"]:
            b["cells"] = [int(round(c * fct)) for c in b["cells"]]
        cfg["output"] = {"vtk": False, "checkpoint": False, "forces": ["plate"]}
        cfg["solver"]["max_iter"] = 4000
        s, S, wt = _run(cfg, base)
        xf, tau, _ = S.wall_shear("plate")
        q = s["dynamic_pressure"]
        x, cf = xf[:, 0], tau / q
        o = np.argsort(x)
        x, cf = x[o], cf[o]
        Re = s["reynolds"]
        bl = blasius_cf(Re * x)
        err = cf / bl - 1
        row = {"cells": S.nc, "iterations": s["iterations"], "converged": s["converged"],
               "time_s": round(wt, 1), "Cd": s["plate"]["Cd"],
               "Cd_blasius": 1.328 / np.sqrt(Re)}
        for xs in (0.1, 0.2, 0.4, 0.6, 0.8, 0.95):
            k = int(np.argmin(np.abs(x - xs)))
            row[f"err_x{xs}"] = float(err[k])
        sel = (x > 0.1) & (x < 0.95)
        row["mean_abs_err"] = float(np.mean(np.abs(err[sel])))
        rows.append(row)
        curves.append((S.nc, x, cf, Re))
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    xr = np.linspace(0.02, 1.0, 200)
    ax.plot(xr, blasius_cf(curves[0][3] * xr), color=REF_COLOR, lw=1.0,
            label="Blasius 0.664/√Re_x")
    for (nc, x, cf, _), col in zip(curves, (BLUE, ORANGE)):
        ax.plot(x, cf, "o", ms=2.5, color=col, label=f"compressible, {nc} cellules")
    ax.set(xlabel="x / L", ylabel="C_f", ylim=(0, 0.02),
           title="Plaque plane laminaire M = 0.2, Re_L = 1e5")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "compressible_plaque.png", dpi=130)
    plt.close(fig)
    return rows


# ------------------------------------------------------------------ 4. Couette
def couette():
    from ..mesh2d.blocks import rectangle_mesh
    from .compressible import (CompressibleSettings, CompressibleSolver2D, Gas,
                               make_state)
    h, U, Tw = 1e-3, 200.0, 300.0
    mesh = rectangle_mesh(0, h / 6, 0, h, 4, 24,
                          names={"left": "l", "right": "r", "bottom": "b", "top": "t"},
                          periodic=[("l", "r")])
    gas = Gas(viscosity="constant", mu=1e-2)
    S = CompressibleSolver2D(mesh, gas, make_state(gas, 0.0, 101325.0, Tw),
                             {"b": {"type": "wall", "T": Tw},
                              "t": {"type": "wall", "T": Tw, "U": [U, 0]}},
                             CompressibleSettings(steady_scheme="implicit", cfl=5,
                                                  max_iter=600, tol=1e-9))
    t0 = time.perf_counter()
    S.run_steady()
    eta = mesh.cell_centers[:, 1] / h
    dT = gas.Pr * U ** 2 / (2 * gas.cp)
    _, tau, _ = S.wall_shear("b")
    _, qw = S.wall_heat_flux("b")
    return {"u_err_max_rel": float(np.abs(S.U[:, 0] - U * eta).max() / U),
            "T_err_max_rel_dTmax": float(np.abs(S.T - (Tw + dT * eta * (1 - eta))).max()
                                         / (dT / 4)),
            "tau_rel_err": float(tau.mean() / (gas.mu * U / h) - 1),
            "q_rel_err": float(qw.mean() / (-gas.mu * U ** 2 / (2 * h)) - 1),
            "iterations": S.iterations, "time_s": round(time.perf_counter() - t0, 2)}


# ------------------------------------------------------------------ 5. NACA 0012
def naca(out: Path, meshes=((96, 32), (192, 64), (384, 128)), mach=0.8, figure=True):
    rows = []
    for na, nr in meshes:
        cfg, base = _cfg("compressible_naca0012_transsonique.toml")
        cfg["mesh"].update(n_around=na, n_radial=nr, first_height=2e-3 * 192 / na)
        cfg["flow"]["mach"] = mach
        cfg["solver"].update(max_iter=3000, tol=1e-7)
        cfg["solver"].pop("monitor_tol", None)
        cfg["output"] = {"vtk": False, "checkpoint": False, "forces": ["airfoil"],
                         "moment_center": [0.25, 0.0]}
        s, S, wt = _run(cfg, base)
        rows.append({"mesh": f"{na}×{nr}", "cells": S.nc, "Cl": s["airfoil"]["Cl"],
                     "Cd": s["airfoil"]["Cd"], "Cm": s["airfoil"]["Cm"],
                     "iterations": s["iterations"], "converged": s["converged"],
                     "final_residual_max": max(s["final_residuals"].values()),
                     "time_s": round(wt, 1)})
        if figure and (na, nr) == (192, 64):
            from ..mesh2d.plot import plot_field
            plt = _pyplot()
            fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
            plot_field(S.mesh, S.mach, ax=axs[0], cmap="viridis", zoom=(-0.4, 1.6, -0.8, 0.8),
                       label="M", title=f"NACA 0012, M∞ = {mach}, α = 1.25° : Mach "
                       f"({S.nc} cellules)")
            sl = S.patch_slices["airfoil"]
            xf = S.mesh.face_centers[S.ni:][sl]
            cp = (S.boundary_p()[sl] - S.fs.p) / s["dynamic_pressure"]
            up = xf[:, 1] >= 0
            axs[1].plot(xf[up, 0], cp[up], "o", ms=2.5, color=BLUE, label="extrados")
            axs[1].plot(xf[~up, 0], cp[~up], "o", ms=2.5, color=ORANGE, label="intrados")
            axs[1].invert_yaxis()
            axs[1].set(xlabel="x / c", ylabel="C_p", title="Coefficient de pression")
            axs[1].grid(True, alpha=0.3)
            axs[1].legend()
            fig.tight_layout()
            fig.savefig(out / f"compressible_naca0012_M{mach:g}.png", dpi=130)
            plt.close(fig)
    return rows


# ------------------------------------------------------------------ rapport
def _table(rows, keys=None, fmt="{:.4g}"):
    if not rows:
        return ""
    keys = keys or list(rows[0])
    lines = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    for r in rows:
        lines.append("| " + " | ".join(
            fmt.format(r[k]) if isinstance(r.get(k), float) else str(r.get(k, ""))
            for k in keys) + " |")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="docs")
    ap.add_argument("--quick", action="store_true", help="maillages réduits")
    ap.add_argument("--skip-naca", action="store_true")
    ap.add_argument("--only", nargs="*", default=None,
                    help="sous-ensemble : sod ramp plate couette naca naca05")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    only = set(a.only) if a.only else {"sod", "ramp", "plate", "couette", "naca", "naca05"}
    if a.skip_naca:
        only -= {"naca", "naca05"}
    report = []
    t0 = time.perf_counter()
    if "sod" in only:
        rows = sod(out, (50, 100) if a.quick else (100, 200, 400))
        report += ["## Sod", _table(rows)]
    if "ramp" in only:
        rows = ramp(out, (0.5,) if a.quick else (0.5, 1.0, 2.0))
        report += ["## Rampe", _table(rows)]
    if "plate" in only:
        rows = plate(out, (1.0,) if a.quick else (1.0, 2.0))
        report += ["## Plaque", _table(rows)]
    if "couette" in only:
        report += ["## Couette", _table([couette()], fmt="{:.3g}")]
    if "naca" in only:
        meshes = ((96, 32),) if a.quick else ((96, 32), (192, 64), (384, 128))
        report += ["## NACA 0012 M = 0.8", _table(naca(out, meshes))]
    if "naca05" in only:
        meshes = ((96, 32),) if a.quick else ((96, 32), (192, 64))
        report += ["## NACA 0012 M = 0.5", _table(naca(out, meshes, mach=0.5, figure=False))]
    report.append(f"\nDurée totale : {time.perf_counter() - t0:.0f} s")
    text = "\n\n".join(report)
    print(text)
    (out / "compressible_validation.md").write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
