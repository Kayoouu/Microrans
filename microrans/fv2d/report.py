"""Résumé lisible d'un calcul (interface, fin de calcul) à partir de summary.json : phrases
en français, 4 chiffres significatifs (avant : JSON brut à 16 chiffres)."""
from __future__ import annotations

from pathlib import Path

MODES = {"steady": "stationnaire", "transient": "instationnaire"}
MODELS = {"laminar": "laminaire", "sa": "Spalart-Allmaras", "ke": "k-ε", "kw": "k-ω",
          "sst": "k-ω SST", "sst_gamma": "k-ω SST + transition γ", "euler": "Euler (non "
          "visqueux)"}
_DONE = {"mode", "model", "n_cells", "axisymmetric", "nu", "reference_velocity",
         "reference_length", "angle_of_attack", "converged", "iterations",
         "iterations_this_run", "wall_time_s", "wall_time_total_s", "U_mean", "backend",
         "checkpoint", "solver", "flux", "order", "limiter", "gas", "freestream", "reynolds",
         "dynamic_pressure", "final_residuals", "time", "steps", "steps_this_run",
         "totals_initial", "totals_final", "viscosity", "porous", "actuator_disks", "scalars",
         "probes", "lines", "animation", "restart", "fmg"}


def _g(v) -> str:
    if isinstance(v, bool) or v is None:
        return {True: "oui", False: "non", None: "—"}[v]
    if isinstance(v, int):
        return f"{v:,}".replace(",", " ")
    if isinstance(v, float):
        return f"{v:.4g}" if v == v else "non défini"
    if isinstance(v, (list, tuple)):
        return "(" + ", ".join(_g(x) for x in v) + ")"
    return str(v)


def _walls(s: dict) -> dict:
    return {k: v for k, v in s.items() if isinstance(v, dict) and ("Cd" in v or "Cd_mean" in v)}


def summary_text(s: dict) -> str:
    """Texte du résumé ; les clés inconnues sont listées à la fin (rien n'est perdu)."""
    comp = s.get("solver") == "compressible"
    out = []
    mode = MODES.get(s.get("mode"), s.get("mode", "?"))
    model = MODELS.get(str(s.get("model")).lower(), s.get("model", "?"))
    head = f"Calcul {'compressible ' if comp else ''}{mode}, {model}"
    if s.get("axisymmetric"):
        head += ", axisymétrique"
    out.append(f"{head} — {_g(s.get('n_cells'))} cellules.")
    if comp:
        out.append(f"Schéma : flux {s.get('flux')}, ordre {s.get('order')}, limiteur "
                   f"{s.get('limiter')}.")
    t = s.get("wall_time_total_s", s.get("wall_time_s"))
    dur = f" ({t:.1f} s)" if isinstance(t, (int, float)) else ""
    if s.get("converged") is True:
        out.append(f"Convergé en {_g(s.get('iterations'))} itérations{dur}.")
    elif s.get("converged") is False:
        out.append(f"NON CONVERGÉ après {_g(s.get('iterations'))} itérations{dur} : "
                   "résultats à vérifier (augmenter max_iter, voir les résidus).")
    elif "time" in s:
        steps = f" en {_g(s['steps'])} pas" if "steps" in s else ""
        out.append(f"Temps simulé t = {_g(s['time'])}{steps}{dur}.")
    elif dur:
        out.append(f"Terminé{dur}.")
    if comp:
        fs = s.get("freestream") or {}
        if fs:
            out.append(f"Amont : Mach {_g(fs.get('mach'))}, vitesse {_g(fs.get('speed'))} "
                       f"m/s, p = {_g(fs.get('p'))} Pa, T = {_g(fs.get('T'))} K, "
                       f"ρ = {_g(fs.get('rho'))} kg/m³"
                       + (f", Re = {_g(s['reynolds'])}" if s.get("reynolds") else "") + ".")
        res = s.get("final_residuals")
        if res:
            out.append("Résidus finaux : " + ", ".join(f"{k} {_g(v)}" for k, v in res.items())
                       + ".")
        ti, tf = s.get("totals_initial"), s.get("totals_final")
        if ti and tf and ti.get("mass"):
            out.append(f"Conservation : masse {_g(tf['mass'] / ti['mass'] - 1.0)} (écart "
                       "relatif final / initial).")
    elif s.get("nu"):
        U, L = s.get("reference_velocity", 1.0), s.get("reference_length", 1.0)
        out.append(f"ν = {_g(s['nu'])} m²/s, U = {_g(U)}, L = {_g(L)} → Re = U L / ν = "
                   f"{_g(U * L / s['nu'])}.")
    if s.get("angle_of_attack"):
        out.append(f"Incidence α = {_g(s['angle_of_attack'])}°.")
    walls = _walls(s)
    if walls:
        out.append("")
        if comp:
            out.append("Efforts (coefficients : force / (½ ρ U² L), par unité de profondeur) :")
        elif s.get("axisymmetric"):
            out.append("Efforts (coefficients : force sur 360° / (½ U² A_ref), A_ref = π L²/4 "
                       "par défaut) :")
        else:
            out.append("Efforts (coefficients : force / (½ U² L), par unité de profondeur) :")
        def c(v, w):                                # bruit d'arrondi (1e-17…) : 0
            ref = max(abs(w.get(k, 0.0) or 0.0) for k in ("Cd", "Cl", "Cm", "Cd_mean"))
            return "0" if isinstance(v, float) and abs(v) < 1e-9 * ref else _g(v)

        yp = any("yplus_max" in w for w in walls.values())
        w0 = max(9, max(len(n) for n in walls) + 1)
        out.append(f"  {'frontière':<{w0}}{'Cd':>10}{'Cl':>10}{'Cm':>10}"
                   + (f"{'y⁺ max':>9}" if yp else ""))
        for name, w in walls.items():
            out.append(f"  {name:<{w0}}{c(w.get('Cd'), w):>10}{c(w.get('Cl'), w):>10}"
                       f"{c(w.get('Cm'), w):>10}"
                       + (f"{_g(w.get('yplus_max')):>9}" if yp else ""))
        for name, w in walls.items():
            if "Cd_pressure" in w:
                visc = w.get("Cd_viscous") or 0.0
                sign = "−" if visc < 0 else "+"
                out.append(f"  {name} : Cd = {c(w.get('Cd_pressure'), w)} (pression) {sign} "
                           f"{c(abs(visc), w)} (frottement).")
            if "strouhal" in w:
                out.append(f"  {name}, sur la 2e moitié du calcul : Cd moyen "
                           f"{_g(w.get('Cd_mean'))}, Cl efficace {_g(w.get('Cl_rms'))}, "
                           f"amplitude de Cl {_g(w.get('Cl_amplitude'))}, Strouhal "
                           f"{_g(w.get('strouhal'))} ({_g(w.get('periods_used'))} périodes).")
            if "Nu_mean" in w:
                out.append(f"  {name} : Nusselt moyen {_g(w['Nu_mean'])}, flux de chaleur "
                           f"{_g(w.get('heat_flux'))}.")
            if "torque" in w:
                out.append(f"  {name} : couple {_g(w['torque'])}.")
            extra = {k: v for k, v in w.items() if k not in (
                "Cd", "Cl", "Cm", "Cd_pressure", "Cd_viscous", "yplus_max", "yplus_mean",
                "strouhal", "Cd_mean", "Cl_rms", "Cl_amplitude", "periods_used", "Nu_mean",
                "heat_flux", "torque")}
            if extra:
                out.append(f"  {name} : " + ", ".join(f"{k} = {_g(v)}" for k, v in extra.items()))
    other = []
    v = s.get("viscosity")
    if isinstance(v, dict):
        other.append(f"Viscosité : {v.get('model')} ; ν de {_g(v.get('nu_min'))} à "
                     f"{_g(v.get('nu_max'))} m²/s ; {_g(v.get('cells_at_nu_max'))} cellule(s) "
                     "à la borne ν max.")
    for z in s.get("porous") or []:
        other.append(f"Zone poreuse « {z.get('name')} » : {_g(z.get('cells'))} cellules, "
                     f"vitesse moyenne {_g(z.get('U_mean'))}, dissipation "
                     f"{_g(z.get('dissipation'))}.")
    for d in s.get("actuator_disks") or []:
        other.append(f"Disque « {d.get('name')} » : poussée {_g(d.get('thrust'))}, couple "
                     f"{_g(d.get('torque'))}, puissance {_g(d.get('power'))}, vitesse au "
                     f"disque {_g(d.get('U_disk'))}.")
    for name, b in (s.get("scalars") or {}).items():
        flux = ", ".join(f"{p} {_g(x['flux'])}" for p, x in b.items()
                         if isinstance(x, dict) and "flux" in x)
        other.append(f"Scalaire « {name} » : flux sortants {flux} ; source totale "
                     f"{_g(b.get('source_total'))}.")
    for p in s.get("probes") or []:
        vals = ", ".join(f"{k} = {_g(x)}" for k, x in p.items() if k not in ("x", "y"))
        other.append(f"Sonde ({_g(p.get('x'))}, {_g(p.get('y'))}) : {vals}.")
    if s.get("U_mean") and not comp:
        other.append(f"Vitesse moyenne dans le domaine : {_g(s['U_mean'])}.")
    if other:
        out.append("")
        out.extend(other)
    files = [Path(p).name for p in (s.get("lines") or [])]
    if s.get("animation"):
        files.append(Path(s["animation"]).name)
    if files:
        out.append("")
        out.append("Fichiers : " + ", ".join(files) + ".")
    rest = {k: v for k, v in s.items() if k not in _DONE and k not in walls}
    if rest:
        out.append("")
        out.append("Autres valeurs : " + "; ".join(f"{k} = {_g(v)}" for k, v in rest.items()))
    out.append("")
    out.append("Toutes les valeurs : summary.json dans le dossier de résultats.")
    return "\n".join(out)
