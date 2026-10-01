"""Cas de calcul compressibles (fichier TOML/JSON) : lecture, calcul, sorties.

Un cas est compressible si [physics] compressible = true ; `run_case` (case.py) lui
délègue alors tout le calcul. Grandeurs SI dimensionnelles.

    [physics]
    compressible = true
    angle_of_attack = 1.25        # ° : oriente l'écoulement amont ; Cd, Cl en axes vent
    reference_length = 1.0        # longueur de référence (corde) des coefficients

    [flow]                        # écoulement amont (= de référence), deux grandeurs parmi
    mach = 0.8                    #   pressure (Pa), temperature (K), density (kg/m³) ;
    pressure = 101325.0           #   mach ou velocity = [u, v] (m/s)
    temperature = 288.15
    gamma = 1.4                   # défauts : air sec (γ = 1.4, R = 287.058, Pr = 0.72)
    gas_constant = 287.058
    prandtl = 0.72
    viscosity = "inviscid"        # inviscid (Euler) | constant | sutherland
    mu = 1.8e-5                   # Pa·s (constant ; μ_ref de Sutherland à T_ref = 273.15 K)
    # reynolds = 1e5              # ou μ déduit de Re = ρ∞ U∞ L_ref / μ∞

    [initial]                     # facultatif (défaut : écoulement amont partout)
    rho = "where(x < 0.5, 1.0, 0.125)"   # deux grandeurs parmi rho, p, T (valeurs ou
    p = "where(x < 0.5, 1e5, 1e4)"       # formules en x, y), U = [ux, uy] ou mach
    U = [0.0, 0.0]
    # restart = "…/checkpoint.npz"

    [boundary.<patch>]            # types : voir compressible.py (farfield, inlet, outlet,
    type = "farfield"             #   supersonic_inlet, supersonic_outlet, slip_wall,
                                  #   symmetry, wall)
    [solver]
    mode = "steady"               # steady | transient
    flux = "roe"                  # roe | hllc
    order = 2
    limiter = "venkatakrishnan"   # venkatakrishnan | barth_jespersen | none
    cfl = 1.5
    steady_scheme = "rk3"         # rk3 | rk5 | implicit
    max_iter = 3000
    tol = 1e-6                    # chute des résidus (stationnaire)
    t_end = 6.3e-4                # instationnaire (s) ; dt = … pour un pas fixe

Sorties : summary.json (Cd, Cl, Cm avec q∞ = ½ρ∞U∞², état amont, Re, convergence, temps),
fields.vtk (ρ, U, p, T, Mach, Cp, entropie, vorticité), wall_<patch>.csv (x, y, p, Cp,
τ_w, Cf, y⁺, T, q), history.csv, figures (Mach, p, ρ, T, Cp pariétal, convergence) et
checkpoint.npz (reprise : même maillage exacte, sinon champs interpolés).
"""
from __future__ import annotations

import csv
import json
import os
import time
from pathlib import Path

import numpy as np

from ..mesh2d.io import write_vtk
from .compressible import (CompressibleSettings, CompressibleSolver2D, Gas, make_state)
from .report import MODES

_SETTINGS = {"flux", "order", "limiter", "venkat_k", "limiter_freeze", "entropy_fix", "cfl",
             "steady_scheme", "cfl_max", "cfl_growth", "first_order_iter",
             "viscous_factor", "max_iter", "tol", "monitor_tol", "monitor_window",
             "log_every", "linear_sweeps", "implicit_jacobian", "linear_solver",
             "linear_iter", "linear_tol", "cfl_adapt", "cfl_cuts"}


def is_compressible(cfg: dict) -> bool:
    return bool(cfg.get("physics", {}).get("compressible", False))


def gas_and_freestream(cfg: dict):
    """(Gas, State) à partir de [flow] (+ [physics] angle_of_attack, reference_length)."""
    fl = dict(cfg.get("flow", {}))
    ph = cfg.get("physics", {})
    known = {"mach", "pressure", "temperature", "density", "velocity", "gamma",
             "gas_constant", "prandtl", "viscosity", "mu", "reynolds", "T_ref",
             "sutherland_S"}
    unknown = set(fl) - known
    if unknown:
        raise ValueError(f"[flow] : clé(s) inconnue(s) {sorted(unknown)} ; attendues : "
                         f"{sorted(known)}")
    visc = str(fl.get("viscosity", "constant" if ("mu" in fl or "reynolds" in fl)
                      else "inviscid")).lower()
    gas = Gas(gamma=float(fl.get("gamma", 1.4)), R=float(fl.get("gas_constant", 287.058)),
              Pr=float(fl.get("prandtl", 0.72)), viscosity=visc,
              mu=float(fl.get("mu", 1.716e-5 if visc == "sutherland" else 1.8e-5)),
              T_ref=float(fl.get("T_ref", 273.15)), S=float(fl.get("sutherland_S", 110.4)))
    if "velocity" not in fl and "mach" not in fl:
        raise ValueError("[flow] : donner mach (ou velocity = [u, v] en m/s).")
    alpha = float(ph.get("angle_of_attack", 0.0))
    fs = make_state(gas, fl.get("mach"), fl.get("pressure"), fl.get("temperature"),
                    fl.get("density"), alpha, fl.get("velocity"))
    if "reynolds" in fl and gas.viscous:
        # μ∞ = ρ∞ U∞ L / Re ; Sutherland : μ_ref ajusté pour μ(T∞) = μ∞
        L = float(ph.get("reference_length", 1.0))
        mu_inf = fs.rho * fs.speed * L / float(fl["reynolds"])
        gas.mu = mu_inf / float(gas.mu_of(fs.T)) * gas.mu if gas.viscosity == \
            "sutherland" else mu_inf
    return gas, fs


def _initial_function(cfg: dict, gas: Gas, fs):
    """Fonction (x, y) -> (ρ, u, v, p) de [initial], ou None (écoulement amont)."""
    ini = {k: v for k, v in cfg.get("initial", {}).items()
           if k not in ("restart", "restart_mode", "restart_shift_U")}
    if not ini:
        return None
    from .compressible import _value
    unknown = set(ini) - {"rho", "density", "p", "pressure", "T", "temperature", "U", "mach"}
    if unknown:
        raise ValueError(f"[initial] : clé(s) inconnue(s) {sorted(unknown)} (compressible : "
                         "rho, p, T, U ou mach).")

    def f(x, y):
        r = ini.get("rho", ini.get("density"))
        p = ini.get("p", ini.get("pressure"))
        T = ini.get("T", ini.get("temperature"))
        r = _value(r, x, y) if r is not None else None
        p = _value(p, x, y) if p is not None else None
        T = _value(T, x, y) if T is not None else None
        if r is None and p is None and T is None:
            r, p = np.full(len(x), fs.rho), np.full(len(x), fs.p)
        elif sum(v is not None for v in (r, p, T)) < 2:
            if p is None:
                p = np.full(len(x), fs.p)
            else:
                T = np.full(len(x), fs.T) if T is None else T
        if r is None:
            r = p / (gas.R * T)
        if p is None:
            p = r * gas.R * T
        if "U" in ini:
            u, v = _value(ini["U"][0], x, y), _value(ini["U"][1], x, y)
        elif "mach" in ini:
            c = np.sqrt(gas.gamma * p / r)
            m = _value(ini["mach"], x, y)
            spd = max(fs.speed, 1e-300)
            u = m * c * (fs.u / spd if fs.speed > 0 else 1.0)
            v = m * c * (fs.v / spd if fs.speed > 0 else 0.0)
        else:
            u, v = np.full(len(x), fs.u), np.full(len(x), fs.v)
        return r, u, v, p
    return f


def build_compressible_solver(cfg: dict, base_dir=".", verbose=False, mesh=None):
    from .case import case_mesh
    if mesh is None:
        mesh = case_mesh(cfg, base_dir, verbose)
    gas, fs = gas_and_freestream(cfg)
    sc = cfg.get("solver", {})
    settings = CompressibleSettings(**{k: v for k, v in sc.items() if k in _SETTINGS})
    bcs = cfg.get("boundary", {})
    solver = CompressibleSolver2D(mesh, gas, fs, bcs, settings,
                                  _initial_function(cfg, gas, fs))
    solver.restart_info = None
    init = cfg.get("initial", {})
    if init.get("restart"):
        path = Path(init["restart"])
        if not path.is_absolute():
            path = Path(base_dir) / path
        if not path.is_file():
            raise FileNotFoundError(f"Fichier de reprise introuvable : {path}")
        solver.restart_info = load_checkpoint(
            solver, path, fields_only=init.get("restart_mode", "exact") == "fields")
        if verbose:
            ri = solver.restart_info
            print(f"Reprise ({ri['mode']}) depuis {path} : itération {ri['iteration']}, "
                  f"t = {ri['time']:.6g}")
    return solver


# ------------------------------------------------------------------ reprise
def save_checkpoint(solver, path, history=None, mode="steady") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {"format": np.array(1), "compressible": np.array(1),
            "cell_centers": np.asarray(solver.mesh.cell_centers),
            "n_faces": np.array([solver.ni, solver.nb]), "Q": solver.Q,
            "time": np.array(solver.time), "iteration": np.array(solver.iterations_total),
            "meta": np.array(json.dumps({
                "steady": mode == "steady", "history": history or [],
                "res0": None if solver._res0 is None else [float(v) for v in solver._res0]},
                default=float))}
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **data)
    os.replace(tmp, path)
    return path


def load_checkpoint(solver, path, fields_only=False) -> dict:
    """Même maillage : état conservatif exact (+ itération, temps, historique) ; sinon
    variables primitives interpolées (restart.py : Delaunay + plus proche voisin)."""
    from .restart import _interpolator, _size
    with np.load(Path(path), allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    if "compressible" not in d:
        raise ValueError(f"{path} : fichier de reprise d'un calcul incompressible.")
    meta = json.loads(str(d["meta"]))
    C = solver.mesh.cell_centers
    same = (d["cell_centers"].shape == C.shape
            and tuple(d["n_faces"]) == (solver.ni, solver.nb)
            and np.allclose(d["cell_centers"], C, rtol=0, atol=1e-9 * _size(C)))
    info = {"file": str(path), "time": float(d["time"]), "iteration": int(d["iteration"])}
    if same:
        solver.Q = d["Q"].copy()
        info["mode"] = "exact" if not fields_only else "champs"
    else:
        W = solver.primitive(d["Q"])
        f = _interpolator(d["cell_centers"], C)
        Wn = np.vstack([f(W[k]) for k in range(4)])
        Wn[0] = np.maximum(Wn[0], 1e-6 * solver.fs.rho)
        Wn[3] = np.maximum(Wn[3], 1e-6 * solver.fs.p)
        solver.Q = solver.conservative(Wn)
        info["mode"] = "interpolé"
    if info["mode"] == "exact":
        solver.time = info["time"]
        solver.iterations_total = info["iteration"]
        if meta.get("steady"):
            solver.history = list(meta.get("history") or [])
            if meta.get("res0"):
                # normalisation des résidus conservée (reprise sans saut)
                solver._res0 = np.array(meta["res0"], float)
                solver._res0_count = 10
        else:
            solver.series_restart = list(meta.get("history") or [])
    return info


# ------------------------------------------------------------------ échantillonnage
class CompressibleSampler:
    """Valeurs en des points (cellule contenante + reconstruction linéaire) : rho, Ux, Uy,
    U_mag, p, T, Mach, Cp."""

    FIELDS = ("rho", "Ux", "Uy", "U_mag", "p", "T", "Mach", "Cp")

    def __init__(self, solver, pts):
        from .sampling import locate
        self.solver = solver
        self.pts = np.asarray(pts, float).reshape(-1, 2)
        self.cell = locate(solver.mesh, self.pts)
        self.ok = self.cell >= 0
        self._c = np.where(self.ok, self.cell, 0)
        self._d = self.pts - solver.mesh.cell_centers[self._c]

    def sample(self, fields=None) -> dict:
        s = self.solver
        W = s.primitive()
        Xb = s.boundary_values(W)
        X = np.vstack([W, W[3] / (W[0] * s.gas.R)])
        gx, gy = s.gradient(X, Xb)
        # gradients limités (Barth-Jespersen, fvm.limit_grad) : pas de dépassement au
        # voisinage des chocs
        for k in range(len(X)):
            g = s.fvm.limit_grad(X[k], np.column_stack([gx[k], gy[k]]), Xb[k])
            gx[k], gy[k] = g[:, 0], g[:, 1]
        c = self._c
        v = X[:, c] + gx[:, c] * self._d[:, 0] + gy[:, c] * self._d[:, 1]
        r, u, vv, p, T = v
        fs = s.fs
        out = {"rho": r, "Ux": u, "Uy": vv, "U_mag": np.hypot(u, vv), "p": p, "T": T,
               "Mach": np.hypot(u, vv) / np.sqrt(s.gas.gamma * np.maximum(p, 1e-300)
                                                 / np.maximum(r, 1e-300))}
        if fs.speed > 0.0:
            out["Cp"] = (p - fs.p) / (0.5 * fs.rho * fs.speed ** 2)
        out = {k: np.where(self.ok, a, np.nan) for k, a in out.items()}
        if fields is not None:
            out = {k: out[k] for k in fields if k in out}
        return out


def _write_lines(solver, lines, out, plot):
    import re
    from .sampling import _plot_line, line_points
    written = []
    for i, ln in enumerate(lines or []):
        name = str(ln.get("name", f"ligne{i + 1}"))
        if not re.fullmatch(r"[\w-]{1,64}", name):
            raise ValueError(f"[[output.lines]] name = {name!r} : lettres, chiffres, _ ou - "
                             f"seulement (nom de fichier).")
        s, pts = line_points(ln["start"], ln["end"], ln.get("n", 200))
        vals = CompressibleSampler(solver, pts).sample(["rho", "Ux", "Uy", "p", "T", "Mach"])
        cols = {"s": s, "x": pts[:, 0], "y": pts[:, 1], **vals}
        path = out / f"line_{name}.csv"
        np.savetxt(path, np.column_stack(list(cols.values())), delimiter=",",
                   header=",".join(cols), comments="", encoding="utf-8")
        written.append(path)
        if plot:
            _plot_line(cols, out / f"line_{name}.png", name)
    return written


# ------------------------------------------------------------------ calcul
def run_compressible_case(cfg: dict, base_dir=".", out_dir=None, verbose=True, plot=True,
                          callback=None, mesh=None, return_solver=False):
    """Équivalent compressible de case.run_case (mêmes arguments, mêmes sorties)."""
    from .case import wind_axes
    from .sampling import parse_points
    ph, sc, oc = cfg.get("physics", {}), cfg.get("solver", {}), cfg.get("output", {})
    if ph.get("axisymmetric"):
        raise ValueError("Compressible : calcul axisymétrique non disponible (plan 2D).")
    out = Path(out_dir or oc.get("directory", "results/compressible"))
    out.mkdir(parents=True, exist_ok=True)
    solver = build_compressible_solver(cfg, base_dir, verbose, mesh=mesh)
    fs, gas = solver.fs, solver.gas
    Lref = float(ph.get("reference_length", 1.0))
    qinf = 0.5 * fs.rho * fs.speed ** 2
    if qinf <= 0.0:
        qinf = float("nan")                       # fluide au repos : pas de coefficients
    kF = 1.0 / (qinf * Lref) if qinf == qinf else float("nan")
    ed, el = wind_axes(cfg)
    center = oc.get("moment_center", (0.0, 0.0))
    walls = [p.name for p in solver.mesh.patches
             if solver.bc_types.get(p.name) in ("wall", "slip_wall")]
    force_patches = oc.get("forces", walls)
    mode = sc.get("mode", "steady")
    s = solver.settings
    Re = (fs.rho * fs.speed * Lref / float(gas.mu_of(fs.T))) if gas.viscous else None
    if verbose:
        print(f"Cas 2D compressible : {solver.nc} cellules, M∞ = {fs.mach:.4g}, "
              f"p∞ = {fs.p:.6g} Pa, T∞ = {fs.T:.5g} K, "
              + (f"Re = {Re:.4g}" if Re else "Euler (non visqueux)")
              + f", flux {s.flux}, ordre {s.order}, limiteur {s.limiter}, {MODES.get(mode, mode)}")
    t0 = time.perf_counter()
    summary = {"solver": "compressible", "mode": mode,
               "model": "euler" if not gas.viscous else "laminar",
               "n_cells": solver.nc, "flux": s.flux, "order": s.order, "limiter": s.limiter,
               "gas": {"gamma": gas.gamma, "R": gas.R, "Pr": gas.Pr,
                       "viscosity": gas.viscosity, "mu": gas.mu},
               "freestream": {k: float(v) for k, v in fs.describe().items()},
               "reynolds": Re, "reference_length": Lref,
               "dynamic_pressure": qinf if qinf == qinf else None,
               "angle_of_attack": float(ph.get("angle_of_attack", 0.0))}
    if solver.restart_info:
        summary["restart"] = solver.restart_info
    probe_pts = parse_points(oc.get("probes"))
    sampler = CompressibleSampler(solver, probe_pts) if len(probe_pts) else None
    probe_fields = ["rho", "Ux", "Uy", "p", "T", "Mach"]

    def probe_rec(sv):
        if sampler is None:
            return {}
        v = sampler.sample(probe_fields)
        return {f"probe{i + 1}_{k}": float(v[k][i]) for i in range(len(probe_pts))
                for k in probe_fields}

    def coeffs(sv):
        sv._last = {}
        rec = {}
        for name, f in sv.forces(force_patches).items():
            rec[f"Cd_{name}"] = float(f["total"] @ ed) * kF
            rec[f"Cl_{name}"] = float(f["total"] @ el) * kF
        return rec

    ckpt = out / "checkpoint.npz" if oc.get("checkpoint", True) else None
    every = 60.0 * float(oc.get("checkpoint_minutes", 5.0))
    last_save = [time.perf_counter()]

    def autosave(sv):
        if ckpt is not None and every > 0 and time.perf_counter() - last_save[0] > every:
            save_checkpoint(sv, ckpt, sv.history if mode == "steady"
                            else sv.series_restart + sv.series, mode)
            last_save[0] = time.perf_counter()

    if mode == "steady":
        def cb(sv, n):
            autosave(sv)
            return callback(sv, n) if callback else None

        def mon(sv):
            return {**coeffs(sv), **probe_rec(sv)} if (force_patches or sampler) else {}
        ok = solver.run_steady(verbose=verbose, log_every=sc.get("log_every", 100),
                               callback=cb, monitor=mon)
        last = solver.history[-1] if solver.history else {}
        summary.update(converged=bool(ok), iterations=solver.iterations_total,
                       iterations_this_run=solver.iterations,
                       final_residuals={k: v for k, v in last.items() if k != "iteration"})
        hist = _merge(solver.history, solver.monitor)
    else:
        if "t_end" not in sc:
            raise ValueError("[solver] mode = \"transient\" : donner t_end (s).")

        def probes(sv):
            return {**(coeffs(sv) if force_patches else {}), **probe_rec(sv)}

        def cb(sv, n):
            autosave(sv)
            return callback(sv, n) if callback else None
        tot0 = solver.totals()
        hist = solver.series_restart + solver.run_transient(
            float(sc["t_end"]), dt=sc.get("dt"), cfl=sc.get("cfl"), verbose=verbose,
            log_every=sc.get("log_every", 100), callback=cb,
            probes=probes if (force_patches or sampler) else None)
        tot1 = solver.totals()
        summary.update(time=solver.time, steps=solver.iterations_total,
                       steps_this_run=solver.iterations,
                       totals_initial=tot0, totals_final=tot1)
    summary["wall_time_s"] = round(time.perf_counter() - t0, 2)
    summary["backend"] = "cpu"
    if ckpt is not None:
        save_checkpoint(solver, ckpt, hist, mode)
        summary["checkpoint"] = str(ckpt)
    solver._last = {}
    for name, f in solver.forces(force_patches).items():
        summary.setdefault(name, {}).update(
            Cd=float(f["total"] @ ed * kF), Cl=float(f["total"] @ el * kF),
            Cm=solver.moment(name, center) * kF / Lref,
            Cd_pressure=float(f["pressure"] @ ed * kF),
            Cd_viscous=float(f["viscous"] @ ed * kF))
        xf, tau, yp = solver.wall_shear(name)
        pb = solver.boundary_p()[solver.patch_slices[name]]
        Tw, qw = solver.wall_heat_flux(name)
        cols = [xf, pb, (pb - fs.p) / qinf, tau, tau / qinf, yp, Tw, qw]
        if gas.viscous:
            summary[name].update(yplus_max=float(yp.max()), yplus_mean=float(yp.mean()))
        np.savetxt(out / f"wall_{name}.csv", np.column_stack(cols), delimiter=",",
                   header="x,y,p,Cp,tau_w,Cf,yplus,T_wall,q", comments="", encoding="utf-8")
    if hist:
        keys = list(dict.fromkeys(k for h in hist for k in h))
        with open(out / "history.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for h in hist:
                w.writerow(h)
    if sampler is not None:
        v = sampler.sample(probe_fields)
        summary["probes"] = [{"x": float(p[0]), "y": float(p[1]),
                              **{k: float(v[k][i]) for k in probe_fields}}
                             for i, p in enumerate(probe_pts)]
    if oc.get("lines"):
        summary["lines"] = [str(p) for p in _write_lines(solver, oc["lines"], out,
                                                         plot and oc.get("plots", True))]
    from ..postprocess import json_safe
    (out / "summary.json").write_text(json.dumps(json_safe(summary), indent=2,
                                                 ensure_ascii=False, default=float),
                                      encoding="utf-8")
    if oc.get("vtk", True):
        write_vtk(solver.mesh, out / "fields.vtk", solver.fields())
    if plot and oc.get("plots", True):
        plot_compressible(solver, hist, out, mode, force_patches, qinf)
    if verbose:                                     # avant : JSON brut de ~60 lignes
        from .report import summary_text
        print("\n" + summary_text(summary))
        print(f"Résultats dans {out.resolve()}")
    return (summary, solver) if return_solver else summary


def _merge(history, monitor):
    """Résidus (chaque itération) + efforts / sondes (toutes les 10 itérations)."""
    mon = {m["iteration"]: m for m in monitor}
    return [{**h, **{k: v for k, v in mon.get(h["iteration"], {}).items()
                     if k != "iteration"}} for h in history]


# ------------------------------------------------------------------ figures
def plot_compressible(solver, hist, out: Path, mode, force_patches, qinf):
    from ..mesh2d.plot import plot_field, plot_mesh
    from ..postprocess import _pyplot
    from .post import _zoom
    plt = _pyplot()
    zoom = _zoom(solver, 2.0)
    f = solver.fields()
    plot_mesh(solver.mesh, out / "mesh.png", zoom=zoom)
    plot_field(solver.mesh, f["Mach"], out / "Mach.png", title="Nombre de Mach", zoom=zoom,
               cmap="jet", label="M")
    plot_field(solver.mesh, f["p"], out / "p.png", title="Pression statique p (Pa)",
               zoom=zoom, cmap="RdBu_r", label="Pa")
    plot_field(solver.mesh, f["rho"], out / "rho.png", title="Masse volumique ρ (kg/m³)",
               zoom=zoom, cmap="viridis", label="kg/m³")
    plot_field(solver.mesh, f["T"], out / "T.png", title="Température T (K)", zoom=zoom,
               cmap="inferno", label="K")
    for name in force_patches:
        sl = solver.patch_slices[name]
        xf = solver.mesh.face_centers[solver.ni:][sl]
        pb = solver.boundary_p()[sl]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(xf[:, 0], (pb - solver.fs.p) / qinf, "o", ms=2.5, color="#2a78d6")
        ax.invert_yaxis()
        ax.set(xlabel="x (m)", ylabel="C_p", title=f"Coefficient de pression sur « {name} »")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(out / f"Cp_{name}.png", dpi=130)
        plt.close(fig)
    if not hist:
        return
    if mode == "steady":
        fig, ax = plt.subplots(figsize=(7, 4))
        it = [h["iteration"] for h in hist]
        for k, ls in zip(("rho", "rhoU", "rhoV", "rhoE"), ("-", "--", "-.", ":")):
            ax.semilogy(it, [max(h[k], 1e-300) for h in hist], ls=ls, label=k)
        ax.set(xlabel="itération", ylabel="résidu relatif (RMS)", title="Convergence")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / "convergence.png", dpi=130)
        plt.close(fig)
    else:
        t = [h["time"] for h in hist]
        for name in force_patches:
            if f"Cd_{name}" not in hist[-1]:
                continue
            fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
            ax[0].plot(t, [h.get(f"Cd_{name}", np.nan) for h in hist], color="#2a78d6")
            ax[0].set_ylabel("C_d")
            ax[1].plot(t, [h.get(f"Cl_{name}", np.nan) for h in hist], color="#2a78d6")
            ax[1].set(xlabel="t (s)", ylabel="C_l")
            ax[0].set_title(f"Efforts sur « {name} »")
            fig.tight_layout()
            fig.savefig(out / f"forces_{name}.png", dpi=130)
            plt.close(fig)
