"""Briques de l'interface : saisie numérique scientifique, liaison formulaire ↔ cas,
zone de tracé matplotlib, exécution en tâche de fond."""
from __future__ import annotations

import contextlib
import io
import time
import traceback

from PySide6.QtCore import QObject, QTimer, Signal, Slot
from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLineEdit, QSpinBox,
                               QVBoxLayout, QWidget)

# ----------------------------------------------------------------------------- saisie


class SciEdit(QLineEdit):
    """Nombre réel en notation libre (1e-6, 0.25...) ; valeur vide = None."""
    valueChanged = Signal()

    def __init__(self, value=None, allow_empty=False, placeholder="", parent=None):
        super().__init__(parent)
        v = QDoubleValidator(self)
        v.setNotation(QDoubleValidator.ScientificNotation)
        self.setValidator(v)
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
        self.setText("" if v is None else (f"{v:g}" if isinstance(v, (int, float)) else str(v)))
        self.blockSignals(False)


class Vec2(QWidget):
    valueChanged = Signal()

    def __init__(self, value=(0.0, 0.0), parent=None):
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.x, self.y = SciEdit(value[0]), SciEdit(value[1])
        for w in (self.x, self.y):
            lay.addWidget(w)
            w.valueChanged.connect(self.valueChanged.emit)

    def value(self):
        return [self.x.value() or 0.0, self.y.value() or 0.0]

    def set_value(self, v):
        v = v if v is not None else (0.0, 0.0)
        self.x.set_value(float(v[0]) if isinstance(v[0], (int, float)) else v[0])
        self.y.set_value(float(v[1]) if isinstance(v[1], (int, float)) else v[1])


def combo(options, parent=None) -> QComboBox:
    """options : [(valeur, libellé), ...]"""
    c = QComboBox(parent)
    for val, label in options:
        c.addItem(label, val)
    return c


def set_combo(c: QComboBox, value):
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
               "combo": "currentIndexChanged", "check": "toggled", "text": "textChanged"}[kind]
        getattr(widget, sig).connect(lambda *_: self.changed.emit())
        return widget

    def sci(self, path, default=None, allow_empty=False, placeholder=""):
        return self._add(path, SciEdit(default, allow_empty, placeholder), default, "sci")

    def vec(self, path, default=(0.0, 0.0)):
        return self._add(path, Vec2(default), default, "vec")

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
                simple = v is None or (len(v) == 2 and all(isinstance(c, (int, float))
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
                set_combo(w, v)
            elif kind == "check":
                w.setChecked(bool(v))
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
            else:
                v = w.text().strip() or None
            self._set(cfg, path, v)


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
        ax = self.axes()
        ax.axis("off")
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
