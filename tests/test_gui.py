"""Interface graphique : auto-test complet hors écran (maillage, calcul, tracés)."""
import os
import subprocess
import sys

import pytest

pytest.importorskip("PySide6.QtWidgets", exc_type=ImportError)


def test_gui_selftest(tmp_path):
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", MICRORANS_RESULTS=str(tmp_path))
    shot = tmp_path / "gui.png"
    r = subprocess.run([sys.executable, "-m", "microrans.gui", "--selftest", str(shot)],
                       env=env, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "SELFTEST OK" in r.stdout
    assert shot.exists()
