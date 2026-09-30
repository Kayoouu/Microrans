"""Cas de calcul 2D décrits par un fichier TOML/JSON (esprit des dictionnaires OpenFOAM /
du fichier .cfg de SU2). Voir les exemples de `microrans/examples/`
(liste : `microrans examples`).

Sections : [mesh] (+ [domain], [[bodies]]), [physics], [initial], [turbulence],
[boundary.<patch>], [solver], [output].

Incidence : [physics] angle_of_attack = α (degrés) fait tourner la vitesse des entrées /
champs lointains et la vitesse initiale ; Cd et Cl sont alors donnés dans les axes de
l'écoulement (traînée selon U∞), Cm autour de [output] moment_center (défaut (0, 0)).
Balayage (polaire) : voir sweep.py.

Axisymétrique : [physics] axisymmetric = true (x = axe, y = rayon ; patch de type « axis »
sur l'axe ; [mesh] cut_axis = true pour ne garder que la moitié y > 0 d'un maillage autour
d'un corps). Efforts et flux de chaleur totaux (sur 360°) ; Cd rapporté à
[physics] reference_area (défaut : maître-couple π L²/4, L = reference_length = diamètre).

Sorties supplémentaires ([output]) : probes = [[x, y], …] (sondes : Ux, Uy, p (, T) à
chaque itération / pas de temps dans history.csv) ; [[output.lines]] name, start, end, n
(profils line_<name>.csv/.png) ; average_from = t (instationnaire : moyennes et écarts-types
Ux_mean, Ux_rms, p_mean…) ; animate = "vorticity" (| U_mag | p | T…), animate_every = N
(instationnaire : animation_<grandeur>.gif, voir animation.py).

Scalaires passifs : [scalars.<nom>] diffusivity (ou schmidt), Sc_t, source, initial, scheme ;
valeurs aux frontières [boundary.<patch>.scalars] <nom> = valeur (défaut : 0 en entrée,
gradient nul ailleurs), flux entrant [boundary.<patch>.scalar_flux] <nom> = q. Bilan par
frontière dans summary.json (« scalars »). Fluide non newtonien : [physics.viscosity]
model = power_law | carreau | cross | herschel_bulkley | bingham | casson (voir
rheology.py) ; ν de [physics] facultatif (défaut : ν(γ̇_ref), γ̇_ref = U_ref / L_ref).

Zones poreuses : [[porous]] name, region = "rectangle" (x0, x1, y0, y1) | "circle" (center,
radius) | "expression" (condition en x, y), darcy = d ou [d1, d2] (1/m², ou permeability
K = 1/d en m²), forchheimer = f ou [f1, f2] (1/m), angle (° : axes principaux de la zone).
Perte de charge par unité de longueur : ν d U + ½ f U² (vitesse superficielle).

Reprise : [initial] restart = "…/checkpoint.npz" (reprise exacte sur le même maillage,
interpolation sinon) ; [output] checkpoint = true (défaut) écrit checkpoint.npz à la fin
du calcul, à l'arrêt demandé et toutes les `checkpoint_minutes` (défaut 5) minutes.
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
from .restart import load_checkpoint, save_checkpoint
from .sampling import Sampler, TimeAverage, parse_points, write_lines
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


def _energy_from(cfg):
    """Section [energy] : Pr, Pr_t, beta, T_ref, gravity, T0, delta_T, source ;
    ou [physics] rayleigh = Ra (convection naturelle adimensionnée, voir exemples)."""
    e = cfg.get("energy")
    if e is None or e.get("enabled", True) is False:
        return None
    e = dict(e)
    e.pop("enabled", None)
    return e


def wind_axes(cfg: dict):
    """Directions de traînée et de portance (vecteurs unitaires) pour l'incidence du cas."""
    a = np.radians(float(cfg.get("physics", {}).get("angle_of_attack", 0.0)))
    return np.array([np.cos(a), np.sin(a)]), np.array([-np.sin(a), np.cos(a)])


def _apply_incidence(cfg: dict) -> tuple[dict, list]:
    """(conditions aux limites, vitesse initiale) tournées de l'incidence du cas."""
    bcs = cfg["boundary"]
    U0 = cfg.get("initial", {}).get("U", (0.0, 0.0))
    alpha = float(cfg.get("physics", {}).get("angle_of_attack", 0.0))
    if alpha == 0.0:
        return bcs, U0
    ed, el = wind_axes(cfg)
    R = np.column_stack([ed, el])                   # repère écoulement -> repère maillage

    def rot(u, where):
        if not all(isinstance(c, (int, float)) for c in u):
            raise ValueError(f"{where} : angle_of_attack exige une vitesse constante "
                             f"(pas d'expression en x, y) ; reçu {u}.")
        return [float(c) for c in R @ np.asarray(u, float)]
    out = {}
    for name, spec in bcs.items():
        spec = dict(spec)
        if spec["type"] in ("inlet", "farfield"):
            spec["U"] = rot(spec["U"], f"[boundary.{name}]")
        out[name] = spec
    return out, rot(U0, "[initial] U")


def case_mesh(cfg: dict, base_dir=".", verbose=False):
    mcfg = cfg.get("mesh", {})
    if "preset" in mcfg:
        return PRESETS[mcfg["preset"]]()
    return build_mesh(cfg, base_dir=base_dir, verbose=verbose)


def build_solver(cfg: dict, base_dir=".", verbose=False, mesh=None):
    """Solveur prêt à calculer (maillage construit, ou fourni via `mesh`)."""
    if mesh is None:
        mesh = case_mesh(cfg, base_dir, verbose)
    ph = cfg.get("physics", {})
    visc = ph.get("viscosity")
    if visc and str(visc.get("model", "newtonian")).lower() != "newtonian":
        visc = {"gamma_ref": ph.get("reference_velocity", 1.0) / ph.get("reference_length", 1.0),
                **visc}
    else:
        visc = None
    if "nu" in ph:
        nu = float(ph["nu"])
    elif "reynolds" in ph:
        nu = ph.get("reference_velocity", 1.0) * ph.get("reference_length", 1.0) / ph["reynolds"]
    elif visc:
        from .rheology import reference_viscosity
        nu = reference_viscosity({k: v for k, v in visc.items() if k != "gamma_ref"},
                                 visc["gamma_ref"])
    else:
        raise ValueError("[physics] : donner nu ou reynolds.")
    init = cfg.get("initial", {})
    axi = bool(ph.get("axisymmetric", False))
    if axi and float(ph.get("angle_of_attack", 0.0)) != 0.0:
        raise ValueError("Axisymétrique : pas d'incidence possible (l'écoulement amont doit "
                         "être parallèle à l'axe).")
    bcs, U0 = _apply_incidence(cfg)
    solver = Solver2D(mesh, nu, bcs, model=ph.get("model", "laminar"),
                      model_options=ph.get("model_options"),
                      body_force=ph.get("body_force", (0.0, 0.0)),
                      initial_U=U0,
                      turbulence_inflow=cfg.get("turbulence"),
                      settings=_settings_from(cfg.get("solver", {})),
                      reference_velocity=ph.get("reference_velocity"),
                      energy=_energy_from(cfg), axisymmetric=axi, viscosity=visc,
                      scalars=cfg.get("scalars"), porous=cfg.get("porous"))
    solver.restart_info = None
    if init.get("restart"):
        path = Path(init["restart"])
        if not path.is_absolute():
            path = Path(base_dir) / path
        if not path.is_file():
            raise FileNotFoundError(f"Fichier de reprise introuvable : {path}")
        solver.restart_info = load_checkpoint(
            solver, path, fields_only=init.get("restart_mode", "exact") == "fields",
            shift_U=init.get("restart_shift_U"))
        if verbose:
            ri = solver.restart_info
            print(f"Reprise ({ri['mode']}) depuis {path} : t = {ri['time']:.6g}, "
                  f"itération {ri['iteration']}"
                  + (f" ; variables absentes (valeurs amont) : {ri['ignored']}"
                     if ri["ignored"] else ""))
        return solver
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
    ph, sc, oc = cfg.get("physics", {}), cfg.get("solver", {}), cfg.get("output", {})
    out = Path(out_dir or oc.get("directory", "results/case2d"))
    out.mkdir(parents=True, exist_ok=True)
    fmg = None
    levels = int(sc.get("fmg_levels", 0))
    if levels > 0 and sc.get("mode", "steady") == "steady" and not cfg.get("initial", {}).get(
            "restart"):
        from .fmg import fmg_initialize
        if verbose:
            print(f"Démarrage multigrille : {levels} niveau(x) grossier(s)")
        t_f = time.perf_counter()
        path, fmg = fmg_initialize(cfg, levels, base_dir, out, verbose, callback)
        if path is not None:
            cfg = {**cfg, "initial": {**cfg.get("initial", {}), "restart": str(path)}}
            fmg = {"levels": fmg, "wall_time_s": round(time.perf_counter() - t_f, 2)}
    solver = build_solver(cfg, base_dir, verbose, mesh=mesh)
    Uref = ph.get("reference_velocity", solver.U_ref)
    Lref = ph.get("reference_length", 1.0)
    qdyn = 0.5 * Uref ** 2 * Lref
    axi = solver.axisymmetric
    # coefficient = effort · kF ; axisymétrique : effort sur 360° (2π × par radian) rapporté
    # à ½ρU²·A_ref ; la résultante radiale est nulle par symétrie (Cl = Cm = 0)
    if axi:
        a_ref = float(ph.get("reference_area", np.pi * Lref ** 2 / 4.0))
        kF = 2.0 * np.pi / (0.5 * Uref ** 2 * a_ref)
    else:
        kF = 1.0 / qdyn
    ed, el = wind_axes(cfg)
    center = oc.get("moment_center", (0.0, 0.0))
    force_patches = oc.get("forces", [p.name for p in solver.mesh.patches if p.type == "wall"])
    mode = sc.get("mode", "steady")
    if verbose:
        print(f"Cas 2D : {solver.mesh.n_cells} cellules, modèle {solver.model.label}, "
              f"ν = {solver.nu:.4g}, mode {mode}")
    t0 = time.perf_counter()
    summary = {"mode": mode, "model": solver.model_name, "n_cells": solver.mesh.n_cells,
               "axisymmetric": axi,
               "nu": solver.nu, "reference_velocity": Uref, "reference_length": Lref,
               "angle_of_attack": float(ph.get("angle_of_attack", 0.0))}
    if solver.restart_info:
        summary["restart"] = {k: solver.restart_info[k] for k in ("file", "mode", "time",
                                                                  "iteration")}
    if fmg:
        summary["fmg"] = fmg
    # sondes : Ux, Uy, p (, T) en des points fixes, à chaque itération / pas de temps
    probe_pts = parse_points(oc.get("probes"))
    sampler = Sampler(solver, probe_pts) if len(probe_pts) else None
    probe_fields = (["Ux", "Uy", "p"] + (["T"] if solver.energy is not None else [])
                    + list(solver.scalars))
    if sampler is not None and not sampler.ok.all() and verbose:
        print(f"  ATTENTION : sondes hors du domaine : {probe_pts[~sampler.ok].tolist()}")

    def probe_rec(s):
        if sampler is None:
            return {}
        v = sampler.sample(probe_fields)
        return {f"probe{i + 1}_{k}": float(v[k][i]) for i in range(len(probe_pts))
                for k in probe_fields}
    if mode != "steady" and oc.get("average_from") is not None:
        solver.averager = TimeAverage(solver, float(oc["average_from"]))
        ra = getattr(solver, "_restart_avg", None)
        if ra is not None and float(ra["avg_t_start"]) == solver.averager.t_start:
            solver.averager.load(ra)
    ckpt = out / "checkpoint.npz" if oc.get("checkpoint", True) else None
    every = 60.0 * float(oc.get("checkpoint_minutes", 5.0))
    last_save = [time.perf_counter()]

    def autosave(s):
        # sauvegarde périodique (temps réel écoulé) : un plantage ou une coupure ne fait
        # perdre que les dernières minutes de calcul
        if ckpt is not None and every > 0 and time.perf_counter() - last_save[0] > every:
            save_checkpoint(s, ckpt, s.history if mode == "steady"
                            else s.series_restart + s.series)
            last_save[0] = time.perf_counter()

    if mode == "steady":
        def cb_steady(s, n):
            autosave(s)
            return callback(s, n) if callback else None
        ok = solver.run_steady(verbose=verbose, log_every=sc.get("log_every", 100),
                               callback=cb_steady, probes=probe_rec if sampler else None)
        summary.update(converged=bool(ok), iterations=solver.iterations_total,
                       iterations_this_run=solver.iterations)
        hist = solver.history
    else:
        def probes(s):
            rec = {}
            for name, f in s.forces(force_patches).items():
                rec[f"Cd_{name}"] = float(f["total"] @ s.backend.asarray(ed)) * kF
                rec[f"Cl_{name}"] = 0.0 if axi else float(f["total"] @ s.backend.asarray(el)) * kF
            rec.update(probe_rec(s))
            return rec
        vtk_every = oc.get("vtk_every", 0)

        n0 = len(solver.series_restart)
        recorder = None
        if oc.get("animate"):
            from .animation import Recorder
            est = max((float(sc["t_end"]) - solver.time) / float(sc["dt"]), 1.0)
            recorder = Recorder(solver, oc["animate"],
                                oc.get("animate_every", max(1, int(round(est / 100)))))

        def cb(s, n):
            if s.averager is not None:
                s.averager.update()
            if recorder is not None:
                recorder.record(s, n0 + n)
            if vtk_every and (n0 + n) % vtk_every == 0:
                write_vtk(s.mesh, out / f"fields_{n0 + n:06d}.vtk", s.fields())
            autosave(s)
            return callback(s, n) if callback else None
        hist = solver.series_restart + solver.run_transient(
            sc["dt"], sc["t_end"], verbose=verbose, log_every=sc.get("log_every", 100),
            probes=probes, callback=cb)
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
    if fmg:
        summary["wall_time_total_s"] = round(summary["wall_time_s"] + fmg["wall_time_s"], 2)
    # vitesse moyenne (pondérée par le volume ; conduite périodique : vitesse débitante)
    summary["U_mean"] = [float(v) for v in solver.backend.to_host(
        solver.xp.sum(solver.U * solver.fvm.V[:, None], axis=0) / solver.xp.sum(solver.fvm.V))]
    summary["backend"] = solver.backend.name
    if solver.rheology is not None:
        nl = solver.backend.to_host(solver.nu_lam)
        summary["viscosity"] = {"model": solver.rheology.describe(),
                                "nu_min": float(nl.min()), "nu_max": float(nl.max()),
                                "cells_at_nu_max": int(np.sum(nl >= 0.999 * solver.rheology.nu_max))}
    if solver.porous is not None:
        summary["porous"] = solver.porous_report()
    if solver.scalars:
        # bilan : flux sortants par frontière (convection + diffusion) et source totale
        summary["scalars"] = {k: solver.scalar_fluxes(k) for k in solver.scalars}
    if ckpt is not None:
        save_checkpoint(solver, ckpt, hist)
        summary["checkpoint"] = str(ckpt)
    solver.to_cpu()                     # post-traitement sur CPU (no-op si déjà CPU)
    for name, f in solver.forces(force_patches).items():
        summary.setdefault(name, {}).update(
            Cd=float(f["total"] @ ed * kF), Cl=0.0 if axi else float(f["total"] @ el * kF),
            Cm=0.0 if axi else solver.moment(name, center) * kF / Lref,
            Cd_pressure=float(f["pressure"] @ ed * kF),
            Cd_viscous=float(f["viscous"] @ ed * kF))
        xf, tau, yp = solver.wall_shear(name)
        summary[name].update(yplus_max=float(yp.max()), yplus_mean=float(yp.mean()))
        # distribution pariétale (comme les « XY plots » de Fluent) : Cf, Cp, y+ (, T, q)
        pb = solver.boundary_p(solver.p)[solver.patch_slices[name]]
        cols = [xf, tau, tau / (0.5 * Uref ** 2), pb / (0.5 * Uref ** 2), yp]
        header = "x,y,tau_w,Cf,Cp,yplus"
        if solver.energy is not None:
            Tw, q = solver.wall_heat_flux(name)
            e = solver.energy
            alpha = solver.nu / float(e["Pr"])
            nu_ref = float(e.get("delta_T", 1.0)) * alpha / Lref
            mag = solver.fvm.magSb[solver.patch_slices[name]]
            summary[name].update(heat_flux=float(np.sum(q * mag)) * (2 * np.pi if axi else 1),
                                 Nu_mean=float(np.sum(q * mag) / np.sum(mag) / nu_ref))
            cols += [Tw, q, q / nu_ref]
            header += ",T_wall,q,Nu"
            if oc.get("nusselt") == "bulk":
                # conduite : Nu local = q L / (α (T_paroi − T_mélange(x))), L = diamètre
                # hydraulique (reference_length)
                Tb = solver.bulk_temperature(xf[:, 0])
                cols += [Tb, q * Lref / (alpha * (Tw - Tb))]
                header += ",T_bulk,Nu_bulk"
        np.savetxt(out / f"wall_{name}.csv", np.column_stack(cols), delimiter=",",
                   header=header, comments="", encoding="utf-8")
    if hist:
        keys = list(dict.fromkeys(k for h in hist for k in h))
        with open(out / "history.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for h in hist:
                w.writerow(h)
    if mode != "steady" and recorder is not None:
        from .post import _body_size, _zoom
        _, Lb = _body_size(solver)
        far = solver.mesh.wall_distance > 0.2 * Lb if Lb else None
        gif = recorder.write(solver.mesh, out / f"animation_{recorder.key}.gif",
                             zoom=_zoom(solver, 14.0),
                             mirror=(-1 if recorder.key in ("vorticity", "Uy") else 1)
                             if axi else None,
                             fps=int(oc.get("animate_fps", 15)), far_mask=far)
        if gif is not None:
            summary["animation"] = str(gif)
    if sampler is not None:
        v = sampler.sample(probe_fields)
        summary["probes"] = [{"x": float(p[0]), "y": float(p[1]),
                              **{k: float(v[k][i]) for k in probe_fields}}
                             for i, p in enumerate(probe_pts)]
    if oc.get("lines"):
        summary["lines"] = [str(p) for p in write_lines(solver, oc["lines"], out,
                                                        plot and oc.get("plots", True))]
    # écrit après les sondes, profils et animation (sinon absents du fichier)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                      encoding="utf-8")
    if oc.get("vtk", True):
        write_vtk(solver.mesh, out / "fields.vtk", solver.fields())
    if plot and oc.get("plots", True):
        from .post import plot_case
        plot_case(solver, hist, out, mode, force_patches, qdyn)
    if verbose:
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        print(f"Résultats dans {out.resolve()}")
    return (summary, solver) if return_solver else summary
