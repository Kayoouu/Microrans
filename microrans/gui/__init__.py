"""Interface graphique (PySide6) : python -m microrans.gui, microrans gui ou microrans-gui."""
import os
import sys
from pathlib import Path


def log_path() -> Path:
    base = Path(os.environ.get("MICRORANS_RESULTS", Path.home() / "microrans_resultats"))
    return base / "microrans-gui.log"


def _install_crash_log():
    """Journal de diagnostic : erreurs Python non rattrapées et plantages bas niveau
    (faulthandler) écrits dans ~/microrans_resultats/microrans-gui.log. Indispensable pour
    l'exécutable Windows, qui n'a pas de console où afficher les erreurs."""
    import datetime
    import faulthandler
    import platform
    import threading
    import traceback

    path = log_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fh = open(path, "a", encoding="utf-8", buffering=1)
    except OSError:
        return None
    from .. import __version__
    fh.write(f"\n=== {datetime.datetime.now():%Y-%m-%d %H:%M:%S} microrans {__version__} "
             f"{platform.platform()} Python {platform.python_version()} "
             f"{'exécutable' if getattr(sys, 'frozen', False) else 'Python'} ===\n")
    faulthandler.enable(file=fh, all_threads=True)
    # application fenêtrée sans console : sys.stdout / sys.stderr valent None
    if sys.stdout is None:
        sys.stdout = fh
    if sys.stderr is None:
        sys.stderr = fh

    def hook(exc_type, exc, tb):
        text = "".join(traceback.format_exception(exc_type, exc, tb))
        fh.write("ERREUR NON RATTRAPÉE\n" + text)
        try:
            from PySide6.QtWidgets import QApplication, QMessageBox
            if QApplication.instance() is not None:
                QMessageBox.critical(None, "microrans : erreur",
                                     f"{exc_type.__name__} : {exc}\n\nDétails enregistrés dans :"
                                     f"\n{path}")
        except Exception:                           # noqa: BLE001
            pass

    sys.excepthook = hook
    threading.excepthook = lambda args: hook(args.exc_type, args.exc_value,
                                             args.exc_traceback)
    return fh


def main(argv=None) -> int:
    fh = _install_crash_log()
    try:
        from .app import run
    except ImportError as exc:                      # PySide6 absent
        print("L'interface graphique nécessite PySide6 : pip install \"microrans[gui]\"\n"
              f"({exc})")
        return 1
    code = run(argv)
    if fh is not None:
        fh.write(f"fermeture normale (code {code})\n")
    return code
