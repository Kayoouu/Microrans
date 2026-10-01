"""Sorties : fichiers CSV/JSON et figures matplotlib (PNG)."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .flow import FlowField
from .reference import log_law, reichardt

# Couleur fixe par modèle (suit l'entité, jamais le rang). Palette validée CVD.
MODEL_COLORS = {
    "sa": "#2a78d6",
    "ke": "#eb6834",
    "kw": "#1baf7a",
    "sst": "#eda100",
    "laminar": "#4a3aa7",
}
REF_COLOR = "#52514e"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
PHASE_RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
              "#2a78d6", "#256abf", "#1c5cab", "#104281"]


# ---------------------------------------------------------------------- utilitaires
def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def json_safe(o):
    """Types NumPy -> Python ; NaN et ±inf -> None (null) : JSON valide pour les autres
    outils (« NaN » n'est pas du JSON standard)."""
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return json_safe(o.tolist())
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        o = o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def write_json(path: Path, data: dict):
    path.write_text(json.dumps(json_safe(data), indent=2, ensure_ascii=False,
                               default=_json_default), encoding="utf-8")


def write_columns(path: Path, columns: dict[str, np.ndarray]):
    names = list(columns)
    data = np.column_stack([np.asarray(columns[k], dtype=float) for k in names])
    np.savetxt(path, data, delimiter=",", header=",".join(names), comments="", fmt="%.10e",
               encoding="utf-8")


def write_table(path: Path, rows: list[dict]):
    if not rows:
        return
    keys = list(rows[0])
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": TEXT, "text.color": TEXT,
        "xtick.color": TEXT_2, "ytick.color": TEXT_2, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "lines.linewidth": 1.6, "font.size": 9,
        "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
        "legend.fontsize": 8,
    })
    return plt


def _color(name: str) -> str:
    return MODEL_COLORS.get(name, "#4a3aa7")


def profile_columns(sol) -> dict[str, np.ndarray]:
    """Profils d'une solution en unités « externes » et « de paroi »."""
    g, nu = sol.grid, sol.nu
    u_tau = max(sol.u_tau, 1e-300)
    flow = FlowField.from_velocity(g, sol.U)
    cols = {
        "y": g.y, "d": g.wall_distance, "y_plus": g.wall_distance * u_tau / nu,
        "U": sol.U, "U_plus": sol.U / u_tau, "nut_over_nu": sol.nut / nu,
        "dUdy": flow.dudy,
    }
    for k, v in sol.state.items():
        cols[k] = v
    for k, v in sol.model.extra_fields(sol.state, flow).items():
        cols[k] = v
    if "k" in sol.state:
        cols["k_plus"] = sol.state["k"] / u_tau ** 2
    return cols


def _lower_half(sol):
    n = sol.grid.n
    return slice(1, n // 2 + 1)


# ---------------------------------------------------------------------------- RANS
def save_rans(result, outdir: Path, plot: bool = True) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    sol = result.solution
    summary = result.summary
    write_json(outdir / "summary.json", summary)
    write_columns(outdir / "profiles.csv", profile_columns(sol))
    if sol.residuals:
        names = list(sol.residuals[0])
        write_columns(outdir / "residuals.csv",
                      {"iteration": np.arange(1, len(sol.residuals) + 1),
                       **{k: [r[k] for r in sol.residuals] for k in names}})
    if plot:
        plot_rans(result, outdir / f"rans_{result.model_name}.png")
    return summary


def plot_rans(result, path: Path):
    plt = _pyplot()
    sol = result.solution
    cols = profile_columns(sol)
    lh = _lower_half(sol)
    yp = cols["y_plus"][lh]
    c = _color(result.model_name)
    re_tau = sol.u_tau * sol.grid.h / sol.nu

    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    fig.suptitle(f"RANS canal plan — {sol.model.label} — Re_τ = {re_tau:.0f}", fontsize=11)

    a = ax[0, 0]
    ref_y = np.logspace(-1, np.log10(yp.max()), 200)
    a.semilogx(ref_y, reichardt(ref_y), color=REF_COLOR, ls="--", lw=1.0, label="Reichardt")
    a.semilogx(ref_y[ref_y > 20], log_law(ref_y[ref_y > 20]), color=REF_COLOR, ls=":",
               lw=1.0, label="log : ln(y⁺)/0.41 + 5.2")
    a.semilogx(yp, cols["U_plus"][lh], color=c, label=result.model_name.upper())
    a.set(xlabel="y⁺", ylabel="U⁺", title="Vitesse moyenne")
    a.legend(loc="upper left")

    a = ax[0, 1]
    a.plot(sol.grid.y / sol.grid.h, cols["nut_over_nu"], color=c)
    a.set(xlabel="y/h", ylabel="ν_t/ν", title="Viscosité turbulente")

    a = ax[1, 0]
    if "k_plus" in cols:
        a.semilogx(yp, cols["k_plus"][lh], color=c, label="k⁺")
        a.set(ylabel="k⁺", title="Énergie cinétique turbulente")
    elif "nu_tilde" in cols:
        a.semilogx(yp, cols["nu_tilde"][lh] / sol.nu, color=c, label="ν̃/ν")
        a.semilogx(yp, 0.41 * yp, color=REF_COLOR, ls="--", lw=1.0, label="κ y⁺")
        a.set(ylabel="ν̃/ν", title="Variable de Spalart-Allmaras")
        a.legend(loc="upper left")
    a.set(xlabel="y⁺")

    a = ax[1, 1]
    if sol.residuals:
        names = list(sol.residuals[0])
        it = np.arange(1, len(sol.residuals) + 1)
        styles = ["-", "--", "-."]
        for i, k in enumerate(names):
            a.semilogy(it, [max(r[k], 1e-300) for r in sol.residuals], color=c,
                       ls=styles[i % 3], label=k)
        a.legend()
    a.set(xlabel="itération", ylabel="variation relative max", title="Convergence")

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def save_rans_comparison(results: list, outdir: Path, plot: bool = True,
                         reference: tuple | None = None) -> list[dict]:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = [r.summary for r in results]
    write_table(outdir / "summary.csv", rows)
    if plot:
        plt = _pyplot()
        fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
        re_tau = results[0].re_tau
        fig.suptitle(f"RANS canal plan — comparaison des modèles — Re_τ = {re_tau:g}",
                     fontsize=11)
        yp_max = max(profile_columns(r.solution)["y_plus"].max() for r in results)
        ref_y = np.logspace(-1, np.log10(yp_max), 200)
        ax[0].semilogx(ref_y, reichardt(ref_y), color=REF_COLOR, ls="--", lw=1.0,
                       label="Reichardt")
        if reference is not None:
            ax[0].semilogx(reference[0], reference[1], color=TEXT, ls="none", marker="o",
                           ms=3, mfc="none", label="référence (fichier)")
        for r in results:
            cols = profile_columns(r.solution)
            lh = _lower_half(r.solution)
            c = _color(r.model_name)
            lbl = r.solution.model.label
            ax[0].semilogx(cols["y_plus"][lh], cols["U_plus"][lh], color=c, label=lbl)
            ax[1].plot(r.solution.grid.y[lh] / r.solution.grid.h,
                       cols["nut_over_nu"][lh], color=c, label=lbl)
            if "k_plus" in cols:
                ax[2].semilogx(cols["y_plus"][lh], cols["k_plus"][lh], color=c, label=lbl)
        ax[0].set(xlabel="y⁺", ylabel="U⁺", title="Vitesse moyenne")
        ax[1].set(xlabel="y/h", ylabel="ν_t/ν", title="Viscosité turbulente (demi-canal)")
        ax[2].set(xlabel="y⁺", ylabel="k⁺", title="k⁺ (modèles à 2 équations)")
        ax[0].legend(loc="upper left")
        ax[1].legend(loc="upper left")
        if ax[2].lines:
            ax[2].legend(loc="upper left")
        fig.tight_layout()
        fig.savefig(outdir / "comparison.png", dpi=130)
        plt.close(fig)
    return rows


# --------------------------------------------------------------------------- URANS
def save_urans(result, outdir: Path, plot: bool = True) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    h = result.history
    summary = result.summary
    write_json(outdir / "summary.json", summary)
    write_columns(outdir / "history.csv", {
        "t": h.times, "forcing": h.forcing, "tau_wall_bottom": h.tau_bottom,
        "tau_wall_top": h.tau_top, "U_bulk": h.bulk_velocity,
        "U_centerline": h.centerline_velocity, "inner_iterations": h.inner_iterations})
    g = result.grid
    phase_cols = {"y": g.y}
    for i, ph in enumerate(result.phase):
        deg = int(round(np.degrees(ph))) % 360
        phase_cols[f"U_phase{deg:03d}"] = result.phase_U[i]
    for i, ph in enumerate(result.phase):
        deg = int(round(np.degrees(ph))) % 360
        phase_cols[f"nut_phase{deg:03d}"] = result.phase_nut[i]
    write_columns(outdir / "phase_profiles.csv", phase_cols)
    lam = result.laminar_stokes_harmonic()
    write_columns(outdir / "harmonic.csv", {
        "y": g.y, "y_plus": g.wall_distance / result.nu,
        "amplitude": np.abs(result.harmonic_U),
        "phase_deg": np.degrees(np.angle(result.harmonic_U / -1.0j)),
        "amplitude_laminar_stokes": np.abs(lam),
        "phase_deg_laminar_stokes": np.degrees(np.angle(lam / -1.0j))})
    write_columns(outdir / "final_profiles.csv", profile_columns(h.final))
    if plot:
        plot_urans(result, outdir / f"urans_{result.model_name}.png")
    return summary


def plot_urans(result, path: Path):
    plt = _pyplot()
    h, g = result.history, result.grid
    c = _color(result.model_name)
    spp = result.steps_per_period
    last = slice(-3 * spp - 1, None)
    t = (h.times[last] - h.times[last][0]) / result.period
    tau = 0.5 * (h.tau_bottom + h.tau_top)
    lh = slice(1, g.n // 2 + 1)
    yp = g.wall_distance[lh] / result.nu
    s = result.summary

    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(f"URANS canal pulsé — {result.steady.model.label} — Re_τ = {result.re_tau:g}, "
                 f"ω⁺ = {s['omega_plus']:.4g}, A = {result.forcing.amplitude:g}", fontsize=11)
    gs = fig.add_gridspec(3, 3)
    # Séries temporelles : un axe par grandeur (jamais de double axe)
    series = [(h.forcing[last], "f(t)", REF_COLOR), (tau[last], "τ_w", c),
              (h.bulk_velocity[last], "U_b", c)]
    for i, (y, lbl, col) in enumerate(series):
        a = fig.add_subplot(gs[i, 0])
        a.plot(t, y, color=col)
        a.set_ylabel(lbl)
        if i == 0:
            a.set_title("3 dernières périodes")
        if i < 2:
            a.tick_params(labelbottom=False)
        else:
            a.set_xlabel("t / T")

    a = fig.add_subplot(gs[:, 1])
    for i, ph in enumerate(result.phase):
        a.plot(result.phase_U[i][lh], g.y[lh] / g.h, color=PHASE_RAMP[i % 8],
               label=f"ωt = {np.degrees(ph):.0f}°")
    a.plot(result.steady.U[lh], g.y[lh] / g.h, color=TEXT, ls="--", lw=1.0,
           label="RANS stationnaire")
    a.set(xlabel="⟨U⟩ (phase)", ylabel="y/h", title="Profils moyennés en phase")
    a.legend(loc="upper left")

    lam = result.laminar_stokes_harmonic()
    a_ref = result.forcing.amplitude / result.forcing.omega
    a = fig.add_subplot(gs[0:2, 2])
    a.semilogx(yp, np.abs(result.harmonic_U[lh]) / a_ref, color=c, label="URANS")
    a.semilogx(yp, np.abs(lam[lh]) / a_ref, color=REF_COLOR, ls="--", lw=1.0,
               label="Stokes laminaire")
    a.set(ylabel="|Û| / (A/ω)", title="1er harmonique de U")
    a.tick_params(labelbottom=False)
    a.legend(loc="upper left")
    a = fig.add_subplot(gs[2, 2])
    a.semilogx(yp, np.degrees(np.angle(result.harmonic_U[lh] / -1.0j)), color=c)
    a.semilogx(yp, np.degrees(np.angle(lam[lh] / -1.0j)), color=REF_COLOR, ls="--", lw=1.0)
    a.set(xlabel="y⁺", ylabel="phase / f (°)")

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def save_urans_comparison(results: list, outdir: Path, plot: bool = True) -> list[dict]:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = [r.summary for r in results]
    write_table(outdir / "summary.csv", rows)
    if plot:
        plt = _pyplot()
        fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
        r0 = results[0]
        s0 = r0.summary
        fig.suptitle(f"URANS canal pulsé — comparaison — Re_τ = {r0.re_tau:g}, "
                     f"ω⁺ = {s0['omega_plus']:.4g}, A = {r0.forcing.amplitude:g}", fontsize=11)
        g = r0.grid
        lh = slice(1, g.n // 2 + 1)
        a_ref = r0.forcing.amplitude / r0.forcing.omega
        lam = r0.laminar_stokes_harmonic()
        yp = g.wall_distance[lh] / r0.nu
        ax[1].semilogx(yp, np.abs(lam[lh]) / a_ref, color=REF_COLOR, ls="--", lw=1.0,
                       label="Stokes laminaire")
        ax[2].semilogx(yp, np.degrees(np.angle(lam[lh] / -1.0j)), color=REF_COLOR,
                       ls="--", lw=1.0, label="Stokes laminaire")
        for r in results:
            h = r.history
            c = _color(r.model_name)
            lbl = r.steady.model.label
            spp = r.steps_per_period
            last = slice(-spp - 1, None)
            ph = (h.times[last] - h.times[last][0]) / r.period * 360.0
            ax[0].plot(ph, 0.5 * (h.tau_bottom[last] + h.tau_top[last]), color=c, label=lbl)
            gg = r.grid
            lhr = slice(1, gg.n // 2 + 1)
            ypr = gg.wall_distance[lhr] / r.nu
            ax[1].semilogx(ypr, np.abs(r.harmonic_U[lhr]) / a_ref, color=c, label=lbl)
            ax[2].semilogx(ypr, np.degrees(np.angle(r.harmonic_U[lhr] / -1.0j)), color=c,
                           label=lbl)
        ax[0].set(xlabel="phase dans la dernière période (°)", ylabel="τ_w",
                  title="Frottement pariétal")
        ax[1].set(xlabel="y⁺", ylabel="|Û| / (A/ω)", title="Amplitude du 1er harmonique")
        ax[2].set(xlabel="y⁺", ylabel="phase / f (°)", title="Phase du 1er harmonique")
        for a in ax:
            a.legend(loc="best")
        fig.tight_layout()
        fig.savefig(outdir / "comparison.png", dpi=130)
        plt.close(fig)
    return rows
