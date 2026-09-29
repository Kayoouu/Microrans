"""Interface graphique (PySide6) : python -m microrans.gui, microrans gui ou microrans-gui."""


def main(argv=None) -> int:
    try:
        from .app import run
    except ImportError as exc:                      # PySide6 absent
        print("L'interface graphique nécessite PySide6 : pip install \"microrans[gui]\"\n"
              f"({exc})")
        return 1
    return run(argv)
