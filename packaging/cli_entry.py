"""Point d'entrée de l'exécutable en ligne de commande (PyInstaller)."""
import multiprocessing
import sys

from microrans.cli import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main())
