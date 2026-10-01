"""Briques de l'interface : saisie numérique scientifique, liaison formulaire ↔ cas,
zone de tracé matplotlib, exécution en tâche de fond."""
from __future__ import annotations

import contextlib
import io
import time
import traceback

from PySide6.QtCore import QObject, QRegularExpression, Qt, QTimer, Signal, Slot
from PySide6.QtGui import QRegularExpressionValidator
from PySide6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLineEdit, QSpinBox,
                               QVBoxLayout, QWidget)

# ----------------------------------------------------------------------------- saisie


# nombre réel, point OU virgule décimale, exposant facultatif ; états partiels (« 1e »,
# « -0, ») acceptés pendant la frappe. QDoubleValidator suit la langue du système : en
# français il supprimait le point en silence (« 0.5 » tapé → 5).
_NUMBER = QRegularExpression(r"^[+-]?(\d*([.,]\d*)?)([eE][+-]?\d*)?$")


class SciEdit(QLineEdit):
    """Nombre réel en notation libre (1e-6, 0.25 ou 0,25...) ; valeur vide = None."""
    valueChanged = Signal()

    def __init__(self, value=None, allow_empty=False, placeholder="", parent=None):
        super().__init__(parent)
        self.setValidator(QRegularExpressionValidator(_NUMBER, self))
        self.allow_empty = allow_empty
        self.setPlaceholderText(placeholder)
        self.set_value(value)
        self.textChanged.connect(lambda _: self.valueChanged.emit())

    def value(self):
        t = self.text().strip().replace(",", ".")
        if not t:
            return None
        try:
            return float(t)
        except ValueError:
            return None

    def set_value(self, v):
        self.blockSignals(True)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            t = f"{v:g}"
            if float(t) != v:
                t = repr(float(v))          # pas d'arrondi (avant : 6 chiffres, réécrits)
        else:
            t = "" if v is None else str(v)
        self.setText(t)
        self.blockSignals(False)


class Vec2(QWidget):
    """Vecteur à 2 composantes (x, y), ou 3 (x, y, z) après set_dim(3) : cas 3D."""
    valueChanged = Signal()

    def __init__(self, value=(0.0, 0.0), parent=None, z_default=0.0):
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.z_default = z_default          # 3e composante absente (cas 2D rouvert en 3D)
        self.x, self.y = SciEdit(value[0]), SciEdit(value[1])
        self.z = SciEdit(value[2] if len(value) > 2 else z_default)
        for w in (self.x, self.y, self.z):
            lay.addWidget(w)
            w.valueChanged.connect(self.valueChanged.emit)
        self.dim = 2
        self.set_dim(len(value))

    def set_dim(self, n: int):
        self.dim = 3 if n == 3 else 2
        self.z.setVisible(self.dim == 3)

    def value(self):
        v = [self.x.value() or 0.0, self.y.value() or 0.0]
        z = self.z.value()
        return v + [self.z_default if z is None else z] if self.dim == 3 else v

    def set_value(self, v):
        v = v if v is not None else (0.0, 0.0)

        def f(c):
            return float(c) if isinstance(c, (int, float)) else c
        self.x.set_value(f(v[0]))
        self.y.set_value(f(v[1]))
        self.z.set_value(f(v[2]) if len(v) > 2 else self.z_default)


def combo(options, parent=None) -> QComboBox:
    """options : [(valeur, libellé), ...]. La liste fermée ne prend pas la largeur du plus
    long libellé (sinon la colonne des réglages déborde) ; la liste ouverte les montre en
    entier."""
    c = QComboBox(parent)
    for val, label in options:
        c.addItem(label, val)
    c.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
    c.setMinimumContentsLength(14)
    fm = c.fontMetrics()
    c.view().setMinimumWidth(max((fm.horizontalAdvance(lab) for _, lab in options),
                                 default=0) + 40)
    return c


def set_combo(c: QComboBox, value) -> bool:
    """Sélectionne l'option de donnée `value` ; False si elle n'existe pas (sélection
    inchangée)."""
    c.blockSignals(True)
    if value is None:                        # option « aucune » (donnée None)
        i = next((k for k in range(c.count()) if c.itemData(k) is None), -1)
    else:
        i = c.findData(value)
    if i < 0 and isinstance(value, str):
        low = [str(c.itemData(k)).lower() for k in range(c.count())]
        i = low.index(value.lower()) if value.lower() in low else -1
    if i >= 0:
        c.setCurrentIndex(i)
    c.blockSignals(False)
    return i >= 0


_EXTRA = Qt.UserRole + 1                     # option ajoutée pour une valeur hors liste


# ----------------------------------------------------------------------------- liaison
class Binder(QObject):
    """Associe des widgets à des chemins dans le dictionnaire du cas (ex. ('solver', 'dt'))."""
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.items = []

    def _add(self, path, widget, default, kind):
        self.items.append((tuple(path), widget, default, kind))
        sig = {"sci": "valueChanged", "vec": "valueChanged", "int": "valueChanged",
               "combo": "currentIndexChanged", "check": "toggled", "text": "textChanged",
               "points": "textChanged"}[kind]
        getattr(widget, sig).connect(lambda *_: self.changed.emit())
        return widget

    def sci(self, path, default=None, allow_empty=False, placeholder=""):
        return self._add(path, SciEdit(default, allow_empty, placeholder), default, "sci")

    def vec(self, path, default=(0.0, 0.0), follow_dim=False, z_default=0.0):
        """follow_dim : 2 ou 3 composantes selon la dimension du cas (vitesse, forces) ;
        z_default : 3e composante quand le cas n'en donne pas (raffinement : 1)."""
        w = Vec2(default, z_default=z_default)
        w.follow_dim = follow_dim
        return self._add(path, w, default, "vec")

    def int(self, path, default=0, lo=0, hi=10 ** 7):
        w = QSpinBox()
        w.setRange(lo, hi)
        w.setValue(default)
        return self._add(path, w, default, "int")

    def combo(self, path, options, default=None):
        w = combo(options)
        set_combo(w, default if default is not None else options[0][0])
        return self._add(path, w, default, "combo")

    def check(self, path, default=False, text=""):
        w = QCheckBox(text)
        w.setChecked(bool(default))
        return self._add(path, w, default, "check")

    def text(self, path, default=""):
        w = QLineEdit(str(default))
        return self._add(path, w, default, "text")

    def points(self, path):
        """Liste de points [[x, y], ...] saisie « x y ; x y » (sondes)."""
        return self._add(path, QLineEdit(), None, "points")

    # -------------------------------------------------------------------------
    @staticmethod
    def _get(cfg, path):
        d = cfg
        for p in path[:-1]:
            d = d.get(p, {}) if isinstance(d, dict) else {}
        return d.get(path[-1]) if isinstance(d, dict) else None

    @staticmethod
    def _set(cfg, path, value):
        d = cfg
        for p in path[:-1]:
            d = d.setdefault(p, {})
        if value is None:
            d.pop(path[-1], None)
        else:
            d[path[-1]] = value

    def load(self, cfg):
        for path, w, default, kind in self.items:
            v = self._get(cfg, path)
            v = default if v is None else v
            w.blockSignals(True)
            if kind == "vec":
                ok_len = (2, 3) if getattr(w, "follow_dim", False) else (w.dim,)
                simple = v is None or (len(v) in ok_len and all(isinstance(c, (int, float))
                                                                for c in v))
                # valeur avancée (ex. multi-grading) : non éditable ici, conservée telle quelle
                w.setEnabled(simple)
                w.setToolTip("" if simple else "Valeur avancée : modifiez-la dans l'onglet TOML")
                if simple:
                    w.set_value(v)
            elif kind == "sci":
                w.set_value(v)
            elif kind == "int":
                w.setValue(int(v or 0))
            elif kind == "combo":
                for k in reversed(range(w.count())):
                    if w.itemData(k, _EXTRA):
                        w.removeItem(k)
                if not set_combo(w, v) and v is not None:
                    # valeur du fichier absente de la liste : gardée telle quelle (avant :
                    # remplacée en silence, ex. linearUpwindLimited -> linearUpwind)
                    w.addItem(f"{v} (valeur du fichier)", v)
                    w.setItemData(w.count() - 1, True, _EXTRA)
                    w.setCurrentIndex(w.count() - 1)
            elif kind == "check":
                w.setChecked(bool(v))
            elif kind == "points":
                w.setText(_points_text(v))
            else:
                w.setText("" if v is None else str(v))
            w.blockSignals(False)

    def store(self, cfg, only=None):
        for path, w, default, kind in self.items:
            if only is not None and path not in only:
                continue
            if not w.isEnabled():
                continue                    # groupe masqué/désactivé : clé inchangée
            if kind in ("sci", "vec"):
                v = w.value()
            elif kind == "int":
                v = int(w.value())
            elif kind == "combo":
                v = w.currentData()
            elif kind == "check":
                v = bool(w.isChecked())
            elif kind == "points":
                v = _points_value(w.text())
            else:
                v = w.text().strip() or None
            self._set(cfg, path, v)


def _points_text(v) -> str:
    """[[x, y], ...] -> « x y ; x y » (3D : « x y z ») ; texte laissé tel quel."""
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    try:
        return " ; ".join(" ".join(f"{float(c):g}" for c in p) for p in v)
    except (TypeError, ValueError):
        return str(v)


def _points_value(text):
    """« x y ; x y » (3D : « x y z ; … ») -> [[x, y], ...] ; texte illisible conservé
    (erreur claire au calcul)."""
    from ..fv2d.sampling import parse_points
    t = text.strip()
    if not t:
        return None
    for dim in (2, 3):
        try:
            return [[float(c) for c in p] for p in parse_points(t, dim)]
        except ValueError:
            continue
    return t


# ----------------------------------------------------------------------------- tracés
class PlotCanvas(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        from matplotlib.figure import Figure
        self.fig = Figure(figsize=(7, 5), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas, 1)

    def axes(self, n=1, **kw):
        self.fig.clear()
        if n == 1:
            return self.fig.add_subplot(111, **kw)
        return self.fig.subplots(1, n, **kw)

    def draw(self):
        self.canvas.draw_idle()

    def message(self, text):
        import textwrap
        ax = self.axes()
        ax.axis("off")
        text = "\n".join(textwrap.fill(par, 60) for par in text.split("\n"))
        ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=11, color="#52514e")
        self.draw()


def apply_plot_style():
    import matplotlib as mpl

    from ..postprocess import GRID, SURFACE, TEXT, TEXT_2
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": TEXT, "text.color": TEXT,
        "xtick.color": TEXT_2, "ytick.color": TEXT_2, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.6, "font.size": 9, "axes.titlesize": 10,
        "legend.frameon": False, "legend.fontsize": 8})


# ----------------------------------------------------------------------------- tâches
class _SignalStream(io.TextIOBase):
    def __init__(self, signal):
        self.signal = signal
        self.buf = ""

    def write(self, s):
        self.buf += s
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            self.signal.emit(line)
        return len(s)

    def flush(self):
        if self.buf:
            self.signal.emit(self.buf)
            self.buf = ""


class Worker(QObject):
    """Exécute fn(worker) hors du fil graphique. fn peut appeler worker.report(dict) et
    consulter worker.stop_requested ; stdout/stderr vont dans le journal."""
    progress = Signal(object)
    log = Signal(str)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn):
        super().__init__()
        self.fn = fn
        self.stop_requested = False
        self._last = 0.0

    def report(self, data, force=False):
        now = time.perf_counter()
        if force or now - self._last > 0.25:
            self._last = now
            self.progress.emit(data)

    @Slot()
    def run(self):
        stream = _SignalStream(self.log)
        try:
            with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("always")
                    warnings.showwarning = (lambda msg, cat, *a, **k:
                                            print(f"Attention : {msg}"))
                    res = self.fn(self)
            stream.flush()
            self.done.emit(res)
        except Exception as exc:                    # noqa: BLE001 — remonté à l'interface
            stream.flush()
            self.failed.emit(f"{type(exc).__name__} : {exc}\n\n{traceback.format_exc()}")


def debounce(parent, ms, fn) -> QTimer:
    t = QTimer(parent)
    t.setSingleShot(True)
    t.setInterval(ms)
    t.timeout.connect(fn)
    return t
