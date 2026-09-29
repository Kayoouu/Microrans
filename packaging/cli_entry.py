"""Point d'entrée de l'exécutable en ligne de commande (PyInstaller).

Double-clic sur `microrans(.exe)` sans argument : on ouvre l'interface graphique au lieu
d'afficher l'aide dans une console qui se refermerait aussitôt (piège classique pour les
débutants). Avec des arguments (`microrans run2d ...`), comportement ligne de commande.
"""
import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()
    if len(sys.argv) == 1 and getattr(sys, "frozen", False):
        if sys.platform == "win32":
            import ctypes
            ctypes.windll.kernel32.FreeConsole()        # ferme la fenêtre noire
            sys.stdout = sys.stderr = None              # → journal microrans-gui.log
        from microrans.gui import main as gui_main
        sys.exit(gui_main([]))
    from microrans.cli import main
    sys.exit(main())
