"""Point d'entrée de l'exécutable graphique (PyInstaller)."""
import multiprocessing
import sys

from microrans.gui import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main())
