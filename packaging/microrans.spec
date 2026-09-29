# -*- mode: python -*-
# Construction : pyinstaller packaging/microrans.spec  (depuis la racine du dépôt)
# Produit dist/microrans/ avec deux programmes partageant les mêmes bibliothèques :
#   microrans-gui(.exe)  interface graphique
#   microrans(.exe)      ligne de commande (microrans run2d ..., microrans rans ...)
import os
from PyInstaller.utils.hooks import collect_submodules

ROOT = os.path.abspath(os.path.join(SPECPATH, ".."))
hidden = collect_submodules("microrans") + ["matplotlib.backends.backend_qtagg"]
datas = [(os.path.join(ROOT, "microrans", "examples"), "microrans/examples"),
         (os.path.join(ROOT, "microrans", "gui", "icon.png"), "microrans/gui")]
ICON = os.path.join(SPECPATH, "microrans.ico")
VERSION = os.path.join(SPECPATH, "version_info.txt")   # ressource de version (Windows)
excludes = ["tkinter", "numba", "llvmlite", "IPython", "pytest", "PySide6.QtWebEngineCore",
            "PySide6.QtWebEngineWidgets", "PySide6.Qt3DCore", "PySide6.QtQuick",
            "PySide6.QtQml", "PySide6.QtMultimedia", "PySide6.QtCharts",
            "PySide6.QtDataVisualization", "PySide6.QtPdf", "PySide6.QtBluetooth",
            "PySide6.QtPositioning", "PySide6.QtSql", "PySide6.QtTest"]


def analysis(script):
    return Analysis([os.path.join(SPECPATH, script)], pathex=[ROOT], datas=datas,
                    hiddenimports=hidden, excludes=excludes, noarchive=False)


a_gui = analysis("gui_entry.py")
a_cli = analysis("cli_entry.py")
# pas de compression UPX (source fréquente de fausses alertes antivirus), icône et
# informations de version renseignées
exe_gui = EXE(PYZ(a_gui.pure), a_gui.scripts, [], exclude_binaries=True,
              name="microrans-gui", console=False, upx=False, icon=ICON, version=VERSION)
exe_cli = EXE(PYZ(a_cli.pure), a_cli.scripts, [], exclude_binaries=True,
              name="microrans", console=True, upx=False, icon=ICON, version=VERSION)
coll = COLLECT(exe_gui, a_gui.binaries, a_gui.datas, exe_cli, a_cli.binaries, a_cli.datas,
               name="microrans", upx=False)
