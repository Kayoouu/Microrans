"""Figures d'un calcul 2D (3D : plan x, y ou z = cte, défaut z médian) : champs,
convergence, efforts."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..mesh2d.plot import plot_field, plot_mesh
from ..postprocess import MODEL_COLORS, REF_COLOR, _pyplot, log_values


def _body_size(solver):
    """Sommets des parois (projetés sur x, y en 3D) et taille du corps."""
    walls = [p.name for p in solver.mesh.patches if p.type == "wall"]
    if not walls:
        return None, None
    m = solver.mesh
    nodes = np.concatenate([m.patch_face_nodes(w).ravel() for w in walls])
    pts = m.points[nodes[nodes >= 0], :2]
    return pts, float(max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1])))


def _zoom(solver, downstream: float = 4.0):
    pts, L = _body_size(solver)
    if pts is None:
        return None
    (x0, y0), (x1, y1) = pts.min(axis=0), pts.max(axis=0)
    P = solver.mesh.points
    bb = (*P[:, :2].min(axis=0), *P[:, :2].max(axis=0))
    if L > 0.5 * max(bb[2] - bb[0], bb[3] - bb[1]):
        return None                       # parois = frontières du domaine : vue globale
    return (x0 - 1.0 * L, x1 + downstream * L, y0 - 1.5 * L, y1 + 1.5 * L)


def plot_case(solver, hist, out: Path, mode, force_patches, qdyn, slice_axis="z",
              slice_value=None):
    """Figures du calcul ; 3D : champs dans le plan slice_axis = slice_value ([output],
    défaut : plan z médian), champs complets dans fields.vtk."""
    plt = _pyplot()
    zoom = _zoom(solver, 4.0 if mode == "steady" else 14.0)
    f = solver.fields()
    mir = 1 if solver.axisymmetric else None          # image miroir par rapport à l'axe
    mesh, sel, plane = solver.mesh, slice(None), ""
    axis = str(slice_axis or "z").lower() if solver.dim == 3 else "z"
    if solver.dim == 3:
        # L5 : avant, toujours le plan z médian (conduite carrée : bande de 2 mailles le
        # long de l'axe, la section x = cte est la figure utile)
        from ..mesh3d.slice import slice_mesh
        try:
            mesh = slice_mesh(solver.mesh, axis, slice_value)
        except ValueError as exc:
            mesh = None                               # pas de coupe : figures de champs omises
            import warnings
            warnings.warn(f"Figures de champs omises : {exc} Champs complets : fields.vtk "
                          "(ParaView).", stacklevel=2)
        else:
            sel, plane = mesh.cells, f" (plan {mesh.axis} = {mesh.value:.4g})"
            if axis != "z":
                zoom = None                           # zoom calculé dans le plan x, y
    if mesh is not None:
        if solver.dim == 2:
            plot_mesh(mesh, out / "mesh.png", zoom=zoom)
        plot_field(mesh, f["U_mag"][sel], out / "U.png", title="|U|" + plane, zoom=zoom,
                   cmap="viridis", mirror=mir)
        plot_field(mesh, f["p"][sel], out / "p.png", title="p (cinématique)" + plane,
                   zoom=zoom, cmap="RdBu_r", mirror=mir)
        omz = solver.grad_U(solver.U)                 # composante normale au plan
        i, j = {"x": (2, 1), "y": (0, 2), "z": (1, 0)}[axis]
        omz = (omz[:, i, j] - omz[:, j, i])[sel]
        # échelle de couleur calée hors couche limite (sinon le sillage paraît délavé)
        _, L = _body_size(solver)
        far = (solver.mesh.wall_distance[sel] > 0.2 * L if L
               else np.ones(len(omz), dtype=bool))
        lim = np.percentile(np.abs(omz[far]) if np.any(far) else np.abs(omz), 99)
        plot_field(mesh, omz, out / "vorticity.png", title=f"vorticité ω_{axis}" + plane,
                   zoom=zoom, cmap="RdBu_r", vmin=-lim, vmax=lim, mirror=mir and -1)
        if "nut_over_nu" in f:
            plot_field(mesh, f["nut_over_nu"][sel], out / "nut.png", title="ν_t/ν" + plane,
                       zoom=zoom, cmap="magma", mirror=mir)
        for name, v in solver.scalars.items():
            plot_field(mesh, v[sel], out / f"scalar_{name}.png",
                       title=f"scalaire {name}" + plane, zoom=zoom, cmap="viridis",
                       mirror=mir)
        if solver.rheology is not None:
            plot_field(mesh, np.log10(solver.nu_lam)[sel], out / "viscosity.png",
                       title="log₁₀ ν (non newtonien)" + plane, zoom=zoom, cmap="magma",
                       mirror=mir)
    if not hist:
        return
    color = MODEL_COLORS.get(solver.model_name, "#4a3aa7")
    if mode == "steady":
        fig, ax = plt.subplots(figsize=(7, 4))
        keys = [k for k in hist[0] if k not in ("iteration",)]
        styles = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]
        it = [h["iteration"] for h in hist]
        for i, k in enumerate(keys):
            ax.semilogy(it, log_values([h.get(k, np.nan) for h in hist]), ls=styles[i % 6],
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
