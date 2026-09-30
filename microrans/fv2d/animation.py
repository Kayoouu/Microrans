"""Animation d'un calcul instationnaire (GIF).

[output] animate = "vorticity" (ou U_mag, Ux, Uy, p, T, nut_over_nu, k, omega…),
animate_every = N pas de temps (défaut : ~100 images sur le calcul), animate_fps = 15.
Les champs sont mémorisés pendant le calcul (float32, quelques Mo), puis les images sont
dessinées à la fin avec une échelle de couleurs FIXE (sinon l'animation « clignote ») et
assemblées par Pillow (fourni avec Matplotlib).
"""
from __future__ import annotations

import numpy as np

LABELS = {"vorticity": "vorticité ω_z", "U_mag": "|U|", "Ux": "U_x", "Uy": "U_y",
          "p": "pression p", "T": "température T", "nut_over_nu": "ν_t / ν"}
SIGNED = ("vorticity", "Uy", "p")            # échelle symétrique, palette divergente


def field_values(solver, key):
    s = solver
    if key == "vorticity":
        g = s.grad_U(s.U)
        v = g[:, 1, 0] - g[:, 0, 1]
    elif key == "Ux":
        v = s.U[:, 0]
    elif key == "Uy":
        v = s.U[:, 1]
    elif key == "U_mag":
        v = s.xp.linalg.norm(s.U, axis=1)
    elif key == "p":
        v = s.p
    elif key == "T" and s.energy is not None:
        v = s.T
    elif key in s.scalars:
        v = s.scalars[key]
    elif key == "viscosity" and s.rheology is not None:
        v = s.nu_lam
    elif key in s.state:
        v = s.state[key]
    elif key == "nut_over_nu":
        v = s.nut / s.nu
    else:
        raise ValueError(f"animate = '{key}' : grandeur inconnue (choix : vorticity, U_mag, "
                         f"Ux, Uy, p, T, nut_over_nu, viscosity, "
                         f"{', '.join(list(s.state) + list(s.scalars))}).")
    return np.asarray(s.backend.to_host(v), dtype=np.float32)


class Recorder:
    def __init__(self, solver, key, every: int = 1, max_frames: int = 400):
        field_values(solver, key)                  # vérifie la grandeur dès le départ
        self.key, self.every, self.max_frames = key, max(int(every), 1), max_frames
        self.frames, self.times = [], []

    def record(self, solver, n):
        if n % self.every == 0 and len(self.frames) < self.max_frames:
            self.frames.append(field_values(solver, self.key))
            self.times.append(float(solver.time))

    def write(self, mesh, path, zoom=None, mirror=None, fps: int = 15, far_mask=None):
        """Écrit le GIF ; renvoie son chemin (None si moins de 2 images)."""
        import matplotlib

        from ..postprocess import SURFACE, TEXT, TEXT_2
        if len(self.frames) < 2:
            return None
        style = {"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "text.color": TEXT,
                 "axes.labelcolor": TEXT, "xtick.color": TEXT_2, "ytick.color": TEXT_2,
                 "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
                 "axes.grid": False}
        with matplotlib.rc_context(style):
            return self._write(mesh, path, zoom, mirror, fps, far_mask)

    def _write(self, mesh, path, zoom, mirror, fps, far_mask):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.collections import PolyCollection
        from matplotlib.figure import Figure
        from PIL import Image

        F = np.stack(self.frames)
        sel = F[:, far_mask] if far_mask is not None and np.any(far_mask) else F
        if self.key in SIGNED:
            lim = float(np.percentile(np.abs(sel), 99)) or 1.0
            clim, cmap = (-lim, lim), "RdBu_r"
        else:
            clim, cmap = (float(np.percentile(sel, 0.5)), float(np.percentile(sel, 99.5))), \
                "viridis"
        polys = [mesh.points[row[:nv]] for row, nv in zip(mesh.cell_nodes, mesh.cell_nv)]
        if mirror:
            polys = polys + [q * np.array([1.0, -1.0]) for q in polys]
        pts = np.concatenate(polys)
        x0, x1, y0, y1 = zoom or (pts[:, 0].min(), pts[:, 0].max(),
                                  pts[:, 1].min(), pts[:, 1].max())
        w = 8.0
        h = float(np.clip(w * (y1 - y0) / max(x1 - x0, 1e-300), 2.5, 7.0)) + 0.6
        fig = Figure(figsize=(w + 1.2, h), dpi=90, layout="constrained")
        FigureCanvasAgg(fig)
        ax = fig.add_subplot()
        pc = PolyCollection(polys, cmap=cmap, edgecolors="face", linewidths=0.05)
        pc.set_clim(*clim)
        ax.add_collection(pc)
        fig.colorbar(pc, ax=ax, shrink=0.85)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        title = ax.set_title("")
        images = []
        for v, t in zip(self.frames, self.times):
            pc.set_array(np.concatenate([v, (mirror or 1) * v]) if mirror else v)
            title.set_text(f"{LABELS.get(self.key, self.key)}    t = {t:.4g}")
            fig.canvas.draw()
            rgba = np.asarray(fig.canvas.buffer_rgba())
            images.append(Image.fromarray(rgba[..., :3].copy()).convert(
                "P", palette=Image.ADAPTIVE, colors=128))
        images[0].save(path, save_all=True, append_images=images[1:],
                       duration=int(1000 / max(fps, 1)), loop=0, optimize=True)
        return path
