"""Figures d'un calcul 2D : champs, convergence, efforts."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..mesh2d.plot import plot_field, plot_mesh
from ..postprocess import MODEL_COLORS, REF_COLOR, _pyplot


def _zoom(solver):
    walls = [p.name for p in solver.mesh.patches if p.type == "wall"]
    if not walls:
        return None
    pts = np.vstack([solver.mesh.points[solver.mesh.patch_face_nodes(w)].reshape(-1, 2)
                     for w in walls])
    (x0, y0), (x1, y1) = pts.min(axis=0), pts.max(axis=0)
    L = max(x1 - x0, y1 - y0)
    bb = solver.mesh.bbox()
    if L > 0.5 * max(bb[2] - bb[0], bb[3] - bb[1]):
        return None                       # parois = frontières du domaine : vue globale
    return (x0 - 1.0 * L, x1 + 4.0 * L, y0 - 1.5 * L, y1 + 1.5 * L)


def plot_case(solver, hist, out: Path, mode, force_patches, qdyn):
    plt = _pyplot()
    zoom = _zoom(solver)
    f = solver.fields()
    plot_mesh(solver.mesh, out / "mesh.png", zoom=zoom)
    plot_field(solver.mesh, f["U_mag"], out / "U.png", title="|U|", zoom=zoom, cmap="viridis")
    plot_field(solver.mesh, f["p"], out / "p.png", title="p (cinématique)", zoom=zoom,
               cmap="RdBu_r")
    w = f["vorticity"]
    lim = np.percentile(np.abs(w), 98)
    omz = solver.grad_U(solver.U)
    plot_field(solver.mesh, omz[:, 1, 0] - omz[:, 0, 1], out / "vorticity.png",
               title="vorticité ω_z", zoom=zoom, cmap="RdBu_r", vmin=-lim, vmax=lim)
    if "nut_over_nu" in f:
        plot_field(solver.mesh, f["nut_over_nu"], out / "nut.png", title="ν_t/ν", zoom=zoom,
                   cmap="magma")
    if not hist:
        return
    color = MODEL_COLORS.get(solver.model_name, "#4a3aa7")
    if mode == "steady":
        fig, ax = plt.subplots(figsize=(7, 4))
        keys = [k for k in hist[0] if k not in ("iteration",)]
        styles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]
        it = [h["iteration"] for h in hist]
        for i, k in enumerate(keys):
            ax.semilogy(it, [max(h.get(k, np.nan), 1e-300) for h in hist], ls=styles[i % 6],
                        color=color if k.startswith("U") else REF_COLOR, label=k)
        ax.set(xlabel="itération", ylabel="résidu normalisé", title="Convergence")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / "convergence.png", dpi=130)
        plt.close(fig)
    else:
        t = np.array([h["time"] for h in hist])
        for name in force_patches:
            fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
            ax[0].plot(t, [h[f"Cd_{name}"] for h in hist], color=color)
            ax[0].set_ylabel("C_d")
            ax[1].plot(t, [h[f"Cl_{name}"] for h in hist], color=color)
            ax[1].set(xlabel="t", ylabel="C_l")
            ax[0].set_title(f"Efforts sur '{name}'")
            fig.tight_layout()
            fig.savefig(out / f"forces_{name}.png", dpi=130)
            plt.close(fig)
