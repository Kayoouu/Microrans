"""Sauvegarde et reprise d'un calcul 2D (fichier `checkpoint.npz`).

- Même maillage : reprise EXACTE (champs, flux aux faces, temps, niveaux de temps
  précédents pour BDF2 / AB2, historique des résidus ou des efforts) — comme les
  répertoires de temps d'OpenFOAM ou les fichiers .dat de Fluent.
- Maillage différent : les champs aux cellules sont interpolés (linéaire par triangulation
  des centres, plus proche voisin hors de l'enveloppe), les flux sont recalculés — comme
  `mapFields` (OpenFOAM) ou « Interpolate » (Fluent). Sert à démarrer un maillage fin depuis
  un calcul grossier, ou un autre modèle depuis une solution laminaire.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

FORMAT = 1


def _host(solver, a):
    return solver.backend.to_host(a)


def save_checkpoint(solver, path, history=None) -> Path:
    """Écrit l'état du solveur. Écriture atomique : un fichier plus ancien n'est jamais
    remplacé par un fichier incomplet (arrêt brutal pendant l'écriture)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    h = lambda a: _host(solver, a)                                        # noqa: E731
    data = {
        "format": np.array(FORMAT),
        "cell_centers": np.asarray(solver.mesh.cell_centers),
        "n_faces": np.array([solver.fvm.ni, solver.fvm.nb]),
        "U": h(solver.U), "p": h(solver.p), "F_i": h(solver.F_i), "F_b": h(solver.F_b),
        "time": np.array(solver.time), "dt": np.array(solver.dt),
        "iteration": np.array(getattr(solver, "iterations_total", 0)),
        "meta": np.array(json.dumps({
            "model": solver.model_name, "steady": bool(solver.steady),
            "time_scheme": solver.settings.time_scheme,
            "state": list(solver.state), "energy": solver.energy is not None,
            "history": history if history is not None else solver.history,
        }, default=float)),
    }
    for k, v in solver.state.items():
        data["state_" + k] = h(v)
    if solver.energy is not None:
        data["T"] = h(solver.T)
    hist = getattr(solver, "_hist", None) or {}
    for lvl, U in enumerate(hist.get("U", [])[:2]):           # niveaux n et n-1 (BDF2)
        data[f"hist{lvl}_U"] = h(U)
        data[f"hist{lvl}_F"] = h(hist["F"][lvl])
        data[f"hist{lvl}_dt"] = np.array(hist["dt"][lvl])
        for k, v in hist["state"][lvl].items():
            data[f"hist{lvl}_state_{k}"] = h(v)
        if hist["T"][lvl] is not None:
            data[f"hist{lvl}_T"] = h(hist["T"][lvl])
    if getattr(solver, "averager", None) is not None:        # moyennes temporelles
        data.update(solver.averager.state())
    if hist.get("R"):                                         # AB2
        data["histR"] = h(hist["R"][0])
        data["histR_dt"] = np.array(hist["dt"][0])
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **data)
    os.replace(tmp, path)
    return path


def read_checkpoint(path) -> dict:
    with np.load(Path(path), allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    if int(d.get("format", 0)) > FORMAT:
        raise ValueError(f"{path} : fichier de reprise d'une version plus récente de microrans.")
    d["meta"] = json.loads(str(d["meta"]))
    return d


def same_mesh(solver, d) -> bool:
    C = np.asarray(solver.mesh.cell_centers)
    return (d["cell_centers"].shape == C.shape
            and tuple(d["n_faces"]) == (solver.fvm.ni, solver.fvm.nb)
            and np.allclose(d["cell_centers"], C, rtol=0, atol=1e-9 * _size(C)))


def _size(C):
    return float(np.max(np.ptp(C, axis=0))) if len(C) else 1.0


def _interpolator(src, dst):
    """Interpolation linéaire (triangulation de Delaunay des centres source), plus proche
    voisin en dehors de l'enveloppe convexe ; renvoie f(valeurs) -> valeurs sur dst."""
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import Delaunay, cKDTree

    tri = Delaunay(src)
    near = cKDTree(src).query(dst)[1]

    def f(v):
        out = LinearNDInterpolator(tri, v)(dst)
        bad = ~np.all(np.isfinite(out.reshape(len(dst), -1)), axis=1)
        out[bad] = v[near[bad]]
        return out
    return f


def load_checkpoint(solver, path, fields_only=False, shift_U=None) -> dict:
    """Initialise `solver` depuis un fichier de reprise ; renvoie un résumé de l'opération
    ({"mode": "exact" | "champs" | "interpolé", ...}).

    fields_only : seuls les champs aux cellules servent d'état initial (nouveau calcul :
    temps, historique et flux repartent de zéro) — continuation d'un balayage de paramètres,
    ou conditions aux limites modifiées. shift_U : vitesse uniforme ajoutée au champ lu
    (continuation en incidence : variation de U∞ entre deux points)."""
    d = read_checkpoint(path)
    meta = d["meta"]
    xp = solver.xp
    A = solver.backend.asarray
    same = same_mesh(solver, d)
    exact = same and not fields_only
    info = {"file": str(path),
            "mode": "exact" if exact else ("champs" if same else "interpolé"),
            "time": float(d["time"]), "iteration": int(d["iteration"]), "ignored": []}
    if same:
        conv = lambda v: v                                               # noqa: E731
    else:
        conv = _interpolator(d["cell_centers"], np.asarray(solver.mesh.cell_centers))
    U = conv(d["U"]).reshape(-1, 2).copy()
    if shift_U is not None:
        U += np.asarray(shift_U, float)
    solver.U = A(U)
    solver.p = A(conv(d["p"]).copy())
    for k in solver.state:
        if "state_" + k in d:
            v = conv(d["state_" + k])
            if k in ("k", "epsilon", "omega"):                # grandeurs positives
                v = np.maximum(v, 1e-12 * np.max(np.abs(v)) + 1e-300)
            solver.state[k] = A(v.copy())
        else:
            info["ignored"].append(k)           # autre modèle : valeur amont conservée
    if solver.energy is not None and "T" in d:
        solver.T = A(conv(d["T"]).copy())
    if exact:
        solver.F_i, solver.F_b = A(d["F_i"].copy()), A(d["F_b"].copy())
        solver.time, solver.dt = float(d["time"]), float(d["dt"])
        solver.iterations_total = int(d["iteration"])
        solver.history = list(meta.get("history") or []) if meta.get("steady") else []
        solver.series_restart = [] if meta.get("steady") else list(meta.get("history") or [])
        hist = {"U": [], "F": [], "state": [], "dt": [], "R": [], "T": []}
        lvl = 0
        while f"hist{lvl}_U" in d:
            hist["U"].append(A(d[f"hist{lvl}_U"].copy()))
            hist["F"].append(A(d[f"hist{lvl}_F"].copy()))
            hist["dt"].append(float(d[f"hist{lvl}_dt"]))
            hist["state"].append({k: A(d[f"hist{lvl}_state_{k}"].copy())
                                  for k in solver.state if f"hist{lvl}_state_{k}" in d})
            hist["T"].append(A(d[f"hist{lvl}_T"].copy()) if f"hist{lvl}_T" in d else None)
            lvl += 1
        if any(len(s) != len(solver.state) for s in hist["state"]):
            hist = {"U": [], "F": [], "state": [], "dt": [], "R": [], "T": []}
        if "histR" in d:
            hist["R"] = [A(d["histR"].copy())]
            if not hist["dt"]:
                hist["dt"] = [float(d["histR_dt"])]
        # toujours transmis (même vide, schémas RK sans historique) : le pas adaptatif
        # reprend alors sa croissance lissée depuis le dernier Δt, comme sans arrêt
        solver._restart_hist = hist
        if "avg_weight" in d:
            solver._restart_avg = {k: v for k, v in d.items() if k.startswith("avg_")}
    else:
        # flux recalculés depuis la vitesse interpolée ; le 1er pas de pression les rend
        # conservatifs
        fvm = solver.fvm
        solver.F_i = xp.sum(fvm.interp(solver.U) * fvm.Si, axis=1)
        solver.F_b = xp.sum(solver.boundary_U(solver.U) * fvm.Sb, axis=1)
        if not same:
            solver.time = float(d["time"])
    solver.update_nut()
    return info
