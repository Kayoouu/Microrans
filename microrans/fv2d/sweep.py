"""Balayage d'un paramètre du cas et polaire Cl(α), Cd(α), Cm(α).

    microrans polar naca0012_polaire --alpha -4 12 2
    microrans sweep cylindre_re20 --param physics.reynolds --values 10 20 40

Chaque point est un calcul complet (sous-dossier) ; par défaut il démarre de la solution du
point précédent (continuation, comme les polaires de Fluent / SU2). En incidence, la
variation de U∞ est ajoutée au champ précédent (le champ lointain est ainsi déjà orienté).
Mesures (nombre d'itérations SIMPLE, continuation / départ du repos) : NACA 0012 SA
Re = 1e6, 2° → 4° : 864 / 902 ; 10° → 12° : 561 / 950 ; laminaire Re = 500, 0° → 3° :
486 / 317. Le gain n'est donc pas systématique (SIMPLE converge en un nombre d'itérations
fixé surtout par le maillage). Au-delà du décrochage, la solution peut dépendre du sens du
balayage (hystérésis, comme en soufflerie) et un calcul stationnaire peut ne pas converger
(décollement instationnaire) : ces points sont marqués « non convergé ».
"""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

import numpy as np

from .case import case_mesh, run_case

ALPHA = "physics.angle_of_attack"
_MESH_KEYS = ("mesh.", "domain.", "bodies")


def parse_values(spec) -> list[float]:
    """Liste, « début:fin:pas » (fin incluse) ou « a, b, c »."""
    if isinstance(spec, (list, tuple)):
        return [float(v) for v in spec]
    spec = str(spec).strip()
    if ":" in spec:
        a, b, h = (float(x) for x in spec.split(":"))
        if h == 0 or (b - a) * h < 0:
            raise ValueError(f"Plage invalide : {spec}")
        n = int(np.floor((b - a) / h + 1e-9)) + 1
        if n > 10000:
            raise ValueError(f"Plage {spec} : {n} points (maximum 10 000).")
        return [round(a + i * h, 12) for i in range(n)]
    return [float(x) for x in spec.replace(";", ",").split(",") if x.strip()]


def set_key(cfg: dict, key: str, value):
    """cfg["physics"]["reynolds"] = value pour key = "physics.reynolds"."""
    parts = key.split(".")
    d = cfg
    for p in parts[:-1]:
        d = d.setdefault(p, {})
    d[parts[-1]] = value


def _label(key: str, v: float) -> str:
    name = "alpha" if key == ALPHA else key.split(".")[-1]
    return f"{name}_{v:+g}".replace("+", "p").replace("-", "m").replace(".", "_")


def run_sweep(cfg: dict, key: str, values, base_dir=".", out_dir=None, continuation=True,
              verbose=True, plot=True, callback=None, mesh=None, on_point=None):
    """Exécute le balayage ; renvoie la liste des lignes du tableau (une par valeur).

    callback(solver, n) -> True : arrêt (arrête aussi le balayage) ; on_point(ligne) :
    appelé après chaque point (interface graphique)."""
    values = parse_values(values)
    out = Path(out_dir or cfg.get("output", {}).get("directory", "results/balayage"))
    out.mkdir(parents=True, exist_ok=True)
    remesh = key.startswith(_MESH_KEYS)
    if mesh is None and not remesh:
        mesh = case_mesh(cfg, base_dir, verbose)          # un seul maillage pour tout
    stop = [False]

    def cb(s, n):
        if callback and callback(s, n):
            stop[0] = True
            return True
        return None

    rows, prev, prev_v = [], None, None
    for v in values:
        c = copy.deepcopy(cfg)
        set_key(c, key, v)
        c.setdefault("output", {})
        c["output"]["plots"] = False
        if continuation and prev is not None:
            c.setdefault("initial", {})["restart"] = str(prev)
            c["initial"]["restart_mode"] = "fields"
            if key == ALPHA:
                c["initial"]["restart_shift_U"] = _freestream_change(cfg, prev_v, v)
        sub = out / _label(key, v)
        if verbose:
            print(f"=== {key} = {v:g} ===")
        summary = run_case(c, base_dir=base_dir, out_dir=sub, verbose=verbose, plot=False,
                           callback=cb, mesh=None if remesh else mesh)
        row = {key: v, "converged": summary.get("converged", True),
               "iterations": summary.get("iterations_this_run", summary.get("iterations"))}
        for name, d in summary.items():
            if isinstance(d, dict) and "Cd" in d:
                for k in ("Cl", "Cd", "Cm", "Cd_pressure", "Cd_viscous", "Cd_mean",
                          "Cl_rms", "strouhal", "Nu_mean"):
                    if k in d:
                        row[f"{k}_{name}"] = d[k]
        rows.append(row)
        prev, prev_v = sub / "checkpoint.npz", v
        if on_point:
            on_point(row)
        if stop[0]:
            break
    write_table(rows, out)
    if plot and rows:
        plot_sweep(rows, key, out / ("polaire.png" if key == ALPHA else "balayage.png"))
    if verbose:
        print(f"Tableau : {out / 'balayage.csv'}")
    return rows


def _freestream_change(cfg, a_old, a_new):
    """U∞(α_new) − U∞(α_old), U∞ = vitesse (non tournée) du 1er champ lointain / entrée."""
    u0 = cfg.get("initial", {}).get("U", (0.0, 0.0))
    for spec in cfg["boundary"].values():
        if spec["type"] in ("farfield", "inlet") and all(
                isinstance(c, (int, float)) for c in spec.get("U", ())):
            u0 = spec["U"]
            break
    u0 = np.asarray(u0, float)

    def rot(a):
        a = np.radians(a)
        return np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]) @ u0
    return [float(x) for x in rot(a_new) - rot(a_old)]


def write_table(rows, out: Path):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "balayage.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    (out / "balayage.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False),
                                       encoding="utf-8")


COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]


def plot_sweep(rows, key, path, fig=None):
    """Polaire (α) : Cl(α), Cd(α), Cl(Cd), Cm(α). Autre paramètre : coefficients en
    fonction du paramètre. Points non convergés : marqueurs creux."""
    if fig is None:
        from matplotlib.figure import Figure          # sans pyplot : sûr dans l'interface
        fig = Figure(figsize=(10, 7.5), layout="constrained")
    x = np.array([r[key] for r in rows])
    conv = np.array([bool(r["converged"]) for r in rows])
    patches = sorted({k.split("_", 1)[1] for r in rows for k in r
                      if k.startswith("Cl_") and not k.startswith("Cl_rms")})
    xl = "incidence α (°)" if key == ALPHA else key
    panels = ([("Cl", xl), ("Cd", xl), ("polar", "Cd"), ("Cm", xl)] if key == ALPHA
              else [("Cl", xl), ("Cd", xl), ("Cm", xl)])
    axes = fig.subplots(2, 2).ravel()
    for ax, (what, xlabel) in zip(axes, panels):
        for i, p in enumerate(patches):
            col = COLORS[i % len(COLORS)]
            if what == "polar":
                xs = np.array([r.get(f"Cd_{p}", np.nan) for r in rows])
                ys = np.array([r.get(f"Cl_{p}", np.nan) for r in rows])
                ylab = "Cl"
            else:
                xs = x
                ys = np.array([r.get(f"{what}_{p}", np.nan) for r in rows])
                ylab = what
            ax.plot(xs, ys, "-", color=col, lw=1.5, label=p)
            ax.plot(xs[conv], ys[conv], "o", color=col, ms=5)
            if np.any(~conv):
                ax.plot(xs[~conv], ys[~conv], "o", mfc="none", color=col, ms=7,
                        label="non convergé" if i == 0 else None)
        ax.set(xlabel=xlabel, ylabel=ylab)
        ax.grid(True, alpha=0.3)
    for ax in axes[len(panels):]:
        ax.set_visible(False)
    if len(patches) > 1 or np.any(~conv):
        axes[0].legend(loc="best")
    fig.suptitle("Polaire" if key == ALPHA else f"Balayage de {key}")
    if path is not None:
        fig.savefig(path, dpi=130)
    return fig
