"""Cas de calcul 2D décrits par un fichier TOML/JSON (esprit des dictionnaires OpenFOAM /
du fichier .cfg de SU2). Voir les exemples de `microrans/examples/`
(liste : `microrans examples`).

Sections : [mesh] (+ [domain], [[bodies]]), [physics], [initial], [turbulence],
[boundary.<patch>], [solver], [output].
"""
from __future__ import annotations

import csv
import json
import time
from dataclasses import fields as dc_fields
from pathlib import Path

import numpy as np

from ..mesh2d.builder import PRESETS, build_mesh
from ..mesh2d.io import write_vtk
from .solver import Settings, Solver2D


def _settings_from(cfg: dict) -> Settings:
    names = {f.name for f in dc_fields(Settings)}
    return Settings(**{k: v for k, v in cfg.items() if k in names})


def strouhal(times, signal, frac: float = 0.5):
    """Fréquence par passages à zéro montants de (signal − moyenne) sur la fin du signal."""
    t = np.asarray(times)
    x = np.asarray(signal)
    sel = t >= t[0] + (1 - frac) * (t[-1] - t[0])
    t, x = t[sel], x[sel] - x[sel].mean()
    z = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    if len(z) < 3:
        return float("nan"), 0
    tz = t[z] - x[z] * (t[z + 1] - t[z]) / (x[z + 1] - x[z])
    return float(1.0 / np.mean(np.diff(tz))), len(z) - 1


def build_solver(cfg: dict, base_dir=".", verbose=False, mesh=None):
    """Solveur prêt à calculer (maillage construit, ou fourni via `mesh`)."""
    mcfg = cfg.get("mesh", {})
    if mesh is not None:
        pass
    elif "preset" in mcfg:
        mesh = PRESETS[mcfg["preset"]]()
    else:
        mesh = build_mesh(cfg, base_dir=base_dir, verbose=verbose)
    ph = cfg.get("physics", {})
    if "nu" in ph:
        nu = float(ph["nu"])
    elif "reynolds" in ph:
        nu = ph.get("reference_velocity", 1.0) * ph.get("reference_length", 1.0) / ph["reynolds"]
    else:
        raise ValueError("[physics] : donner nu ou reynolds.")
    init = cfg.get("initial", {})
    solver = Solver2D(mesh, nu, cfg["boundary"], model=ph.get("model", "laminar"),
                      model_options=ph.get("model_options"),
                      body_force=ph.get("body_force", (0.0, 0.0)),
                      initial_U=init.get("U", (0.0, 0.0)),
                      turbulence_inflow=cfg.get("turbulence"),
                      settings=_settings_from(cfg.get("solver", {})),
                      reference_velocity=ph.get("reference_velocity"))
    amp = init.get("perturbation", 0.0)
    if amp:
        # tourbillon gaussien dans le sillage pour déclencher une instabilité (lâcher)
        c = np.asarray(init.get("perturbation_center", (1.5, 0.0)))
        C = mesh.cell_centers
        solver.U[:, 1] += solver.backend.asarray(
            amp * solver.U_ref * np.exp(-np.sum((C - c) ** 2, axis=1)))
    return solver


def run_case(cfg: dict, base_dir=".", out_dir=None, verbose=True, plot=True, callback=None,
             mesh=None, return_solver=False):
    """Exécute un cas complet. callback(solver, n) -> True pour arrêter (interface
    graphique) ; mesh : maillage déjà construit ; return_solver : renvoie (résumé, solveur)."""
    solver = build_solver(cfg, base_dir, verbose, mesh=mesh)
    ph, sc, oc = cfg.get("physics", {}), cfg.get("solver", {}), cfg.get("output", {})
    out = Path(out_dir or oc.get("directory", "results/case2d"))
    out.mkdir(parents=True, exist_ok=True)
    Uref = ph.get("reference_velocity", solver.U_ref)
    Lref = ph.get("reference_length", 1.0)
    qdyn = 0.5 * Uref ** 2 * Lref
    force_patches = oc.get("forces", [p.name for p in solver.mesh.patches if p.type == "wall"])
    mode = sc.get("mode", "steady")
    if verbose:
        print(f"Cas 2D : {solver.mesh.n_cells} cellules, modèle {solver.model.label}, "
              f"ν = {solver.nu:.4g}, mode {mode}")
    t0 = time.perf_counter()
    summary = {"mode": mode, "model": solver.model_name, "n_cells": solver.mesh.n_cells,
               "nu": solver.nu, "reference_velocity": Uref, "reference_length": Lref}
    if mode == "steady":
        ok = solver.run_steady(verbose=verbose, log_every=sc.get("log_every", 100),
                               callback=callback)
        summary.update(converged=bool(ok), iterations=solver.iterations)
        hist = solver.history
    else:
        def probes(s):
            rec = {}
            for name, f in s.forces(force_patches).items():
                rec[f"Cd_{name}"] = f["total"][0] / qdyn
                rec[f"Cl_{name}"] = f["total"][1] / qdyn
            return rec
        vtk_every = oc.get("vtk_every", 0)

        def cb(s, n):
            if vtk_every and n % vtk_every == 0:
                write_vtk(s.mesh, out / f"fields_{n:06d}.vtk", s.fields())
            return callback(s, n) if callback else None
        hist = solver.run_transient(sc["dt"], sc["t_end"], verbose=verbose,
                                    log_every=sc.get("log_every", 100), probes=probes,
                                    callback=cb)
        t = np.array([h["time"] for h in hist])
        for name in (force_patches if len(hist) > 2 else []):
            cd = np.array([h[f"Cd_{name}"] for h in hist])
            cl = np.array([h[f"Cl_{name}"] for h in hist])
            f, n = strouhal(t, cl)
            late = t >= t[0] + 0.5 * (t[-1] - t[0])
            summary[name] = {"Cd_mean": float(cd[late].mean()),
                             "Cl_rms": float(np.std(cl[late])),
                             "Cl_amplitude": float(0.5 * (cl[late].max() - cl[late].min())),
                             "strouhal": f * Lref / Uref, "periods_used": n}
    summary["wall_time_s"] = round(time.perf_counter() - t0, 2)
    summary["backend"] = solver.backend.name
    solver.to_cpu()                     # post-traitement sur CPU (no-op si déjà CPU)
    for name, f in solver.forces(force_patches).items():
        summary.setdefault(name, {}).update(
            Cd=float(f["total"][0] / qdyn), Cl=float(f["total"][1] / qdyn),
            Cd_pressure=float(f["pressure"][0] / qdyn), Cd_viscous=float(f["viscous"][0] / qdyn))
        xf, tau, yp = solver.wall_shear(name)
        summary[name].update(yplus_max=float(yp.max()), yplus_mean=float(yp.mean()))
        # distribution pariétale (comme les « XY plots » de Fluent) : Cf, Cp, y+
        pb = solver.boundary_p(solver.p)[solver.patch_slices[name]]
        np.savetxt(out / f"wall_{name}.csv",
                   np.column_stack([xf, tau, tau / (0.5 * Uref ** 2), pb / (0.5 * Uref ** 2), yp]),
                   delimiter=",", header="x,y,tau_w,Cf,Cp,yplus", comments="")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                      encoding="utf-8")
    if hist:
        keys = list(dict.fromkeys(k for h in hist for k in h))
        with open(out / "history.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for h in hist:
                w.writerow(h)
    if oc.get("vtk", True):
        write_vtk(solver.mesh, out / "fields.vtk", solver.fields())
    if plot and oc.get("plots", True):
        from .post import plot_case
        plot_case(solver, hist, out, mode, force_patches, qdyn)
    if verbose:
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        print(f"Résultats dans {out.resolve()}")
    return (summary, solver) if return_solver else summary
