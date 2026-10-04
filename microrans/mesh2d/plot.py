"""Tracé de maillages et de champs 2D (matplotlib)."""
from __future__ import annotations

import numpy as np

PATCH_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7",
                "#e34948"]


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _figsize(mesh, zoom, width=10.0, extra=0.0):
    if zoom:
        dx, dy = zoom[1] - zoom[0], zoom[3] - zoom[2]
    else:
        b = mesh.bbox()
        dx, dy = b[2] - b[0], b[3] - b[1]
    h = float(np.clip((width - extra) * dy / max(dx, 1e-300) + 0.8, 3.0, 9.0))
    return (width, h)


def plot_mesh(mesh, path=None, ax=None, title=None, zoom=None, linewidth=0.3, show_patches=True):
    """Trace les cellules et colore les patches frontières. zoom = (x0, x1, y0, y1)."""
    from matplotlib.collections import LineCollection, PolyCollection
    plt = _plt()
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=_figsize(mesh, zoom, 9.0))
    ax.grid(False)
    polys = [mesh.points[row[:nv]] for row, nv in zip(mesh.cell_nodes, mesh.cell_nv)]
    ax.add_collection(PolyCollection(polys, facecolors="#f4f3ef", edgecolors="#52514e",
                                     linewidths=linewidth))
    if show_patches:
        for k, (name, ptype, edges) in enumerate(mesh.all_boundary_patches()):
            segs = mesh.points[edges]
            ax.add_collection(LineCollection(segs, colors=PATCH_COLORS[k % len(PATCH_COLORS)],
                                             linewidths=1.8, label=f"{name} ({ptype})"))
        ax.legend(loc="upper right", fontsize=8, frameon=True)
    ax.set_aspect("equal")
    if zoom:
        ax.set_xlim(zoom[0], zoom[1])
        ax.set_ylim(zoom[2], zoom[3])
    else:
        ax.autoscale_view()
    if title is None:                       # coupe d'un maillage 3D : titre fourni
        q = mesh.quality()
        types = ", ".join(f"{v} {k}" for k, v in q["cell_types"].items())
        title = f"{q['n_cells']} cellules ({types})"
    ax.set_title(title, fontsize=10)
    if own and path:
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
    return ax


def plot_field(mesh, values, path=None, ax=None, title=None, cmap="viridis", zoom=None,
               vmin=None, vmax=None, label=None, vectors=None, vector_stride=None,
               mirror=None, vector_ref=None):
    """Champ aux cellules (couleur par cellule), optionnellement des vecteurs.

    mirror : ±1 (calcul axisymétrique) — trace aussi l'image miroir par rapport à l'axe
    y = 0, avec les valeurs multipliées par ce signe (−1 pour u_r, la vorticité).
    vector_ref : vitesse de référence (coupe 3D : |U| complet maximal). Composantes dans
    le plan < 0.1 % de vector_ref : pas de flèches (avant, dans une coupe 3D, du bruit à
    1e-16 était dessiné en flèches pleine longueur, U17) ; sinon échelle automatique et,
    si elles sont petites devant |U|, le titre donne leur maximum en % de |U|."""
    from matplotlib.collections import PolyCollection
    plt = _plt()
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=_figsize(mesh, zoom, 10.0, extra=1.5))
    ax.grid(False)
    polys = [mesh.points[row[:nv]] for row, nv in zip(mesh.cell_nodes, mesh.cell_nv)]
    values = np.asarray(values)
    if mirror:
        polys = polys + [p * np.array([1.0, -1.0]) for p in polys]
        values = np.concatenate([values, mirror * values])
    pc = PolyCollection(polys, array=values, cmap=cmap, edgecolors="face",
                        linewidths=0.05)
    if vmin is not None or vmax is not None:
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    cb = ax.figure.colorbar(pc, ax=ax, shrink=0.85)
    if label:
        cb.set_label(label)
    if vectors is not None:
        C = mesh.cell_centers
        sel = np.arange(0, mesh.n_cells, vector_stride or max(mesh.n_cells // 1500, 1))
        if zoom:
            m = ((C[sel, 0] > zoom[0]) & (C[sel, 0] < zoom[1]) & (C[sel, 1] > zoom[2])
                 & (C[sel, 1] < zoom[3]))
            sel = sel[m]
        if vector_ref:
            mag = np.hypot(vectors[sel, 0], vectors[sel, 1])
            ratio = float(mag.max()) / vector_ref if mag.size else 0.0
            if ratio < 1e-3:
                sel = sel[:0]
                title = (title or "") + " — vitesse dans le plan < 0.1 % de |U| : pas de flèches"
            elif ratio < 0.5:
                title = (title or "") + f" — flèches agrandies (max {100 * ratio:.0f} % de |U|)"
        if len(sel):
            ax.quiver(C[sel, 0], C[sel, 1], vectors[sel, 0], vectors[sel, 1],
                      color="#0b0b0b", scale_units="xy", angles="xy", width=0.0015)
        if mirror and len(sel):
            ax.quiver(C[sel, 0], -C[sel, 1], vectors[sel, 0], -vectors[sel, 1],
                      color="#0b0b0b", scale_units="xy", angles="xy", width=0.0015)
    ax.set_aspect("equal")
    if zoom:
        ax.set_xlim(zoom[0], zoom[1])
        ax.set_ylim(zoom[2], zoom[3])
    else:
        ax.autoscale_view()
    if title:
        ax.set_title(title, fontsize=10)
    if own and path:
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
    return ax
