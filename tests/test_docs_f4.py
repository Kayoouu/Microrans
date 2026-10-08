"""Lot F4 (audit 2) : textes et documentation. Les extraits des docs sont exécutés, et les
messages cités par le dépannage confrontés aux vrais messages."""
import os
import re
import shlex
from pathlib import Path

import pytest

from microrans.cli import main
from microrans.tomlio import loads

ROOT = Path(__file__).resolve().parents[1]


def _section(path, title):
    text = (ROOT / path).read_text(encoding="utf-8")
    start = text.index(title)
    end = text.find("\n#", start + len(title))
    return text[start:end if end > 0 else None]


def test_3d_in_general_descriptions():
    """T1 : `microrans --help` (« Écoulements 2D en volumes finis »), description du paquet
    (« RANS/URANS 1D et 2D ») et « À propos » sans la 3D."""
    from microrans.cli import build_parser
    assert "2D et 3D" in build_parser().description
    meta = loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert "3D" in meta["project"]["description"]
    pytest.importorskip("PySide6.QtWidgets", exc_type=ImportError)
    from microrans.gui.app import about_text
    assert "2D et 3D" in about_text() and "<b>3D</b>" in about_text()


def test_readme_3d_case_runs(tmp_path):
    """D5 : l'extrait du README (§ 3, « Cas 3D ») laissait croire que `patch_types` suffit
    pour des faces z en symétrie ; il faut [boundary.back] et [boundary.front]. Les deux
    variantes du README sont calculées (maillage réduit, 2 itérations)."""
    from microrans.fv2d.case import run_case
    blocks = re.findall(r"```toml\n(.*?)```", _section("README.md", "### Cas 3D"), re.S)
    assert len(blocks) == 2
    for variant in ("périodique", "symétrie"):
        cfg = loads(blocks[0])
        cfg["mesh"].update(n_around=16, n_radial=8)
        cfg["mesh"]["extrude"]["nz"] = 2
        if variant == "symétrie":
            cfg["mesh"]["extrude"].pop("periodic")
            cfg["boundary"].update(loads(blocks[1])["boundary"])
        cfg["solver"] = {"max_iter": 2}
        s = run_case(cfg, out_dir=tmp_path / variant, verbose=False, plot=False)
        assert s["dimension"] == 3 and s["iterations"] == 2, variant


def test_tutorial_3d_commands_run(tmp_path, monkeypatch):
    """D6 : la 3D tenait en un paragraphe du tutoriel ; § 7 « Un cas 3D » : les commandes
    données s'exécutent (maillages réduits)."""
    sec = _section("docs/tutoriel.md", "## 7. Un cas 3D")
    cmds = [line for block in re.findall(r"```\n(.*?)```", sec, re.S)
            for line in block.splitlines() if line.startswith("microrans ")]
    assert "microrans run2d conduite_carree_3d.toml" in cmds
    ext = next(c for c in cmds if "mesh.extrude.nz=1" in c)
    monkeypatch.chdir(tmp_path)
    assert main(["examples", "conduite_carree_3d"]) == 0
    assert main(["examples", "cylindre_re20"]) == 0
    fast = {"conduite_carree_3d.toml": ["mesh.ny=8", "mesh.nz=8", "solver.max_iter=3"],
            "cylindre_re20.toml": ["mesh.n_around=16", "mesh.n_radial=8",
                                   "solver.max_iter=3"]}
    for cmd in ("microrans run2d conduite_carree_3d.toml", ext):
        argv = shlex.split(cmd)[1:]
        if "--set" not in argv:
            argv.append("--set")
        argv[argv.index("--set") + 1:argv.index("--set") + 1] = fast[argv[1]]
        assert main(argv + ["--no-plot", "-q"]) in (0, 1), cmd     # 1 : non convergé


def test_troubleshooting_quotes_real_3d_messages():
    """D7 : pas d'entrée de dépannage pour les messages 3D (M13 à M17) ; chaque extrait cité
    doit être le vrai message."""
    from microrans.fv2d.validate import check_case
    from microrans.mesh2d.mesh import unknown_patch_type
    table = _section("docs/depannage.md", "## 1. Le cas est refusé")
    real = [
        unknown_patch_type("symetrie", "walls"),
        check_case({"mesh": {"type": "file", "path": "x.msh", "periodic": [["a", "b"]]},
                    "physics": {"nu": 0.01}, "boundary": {}}, mesh_only=True)[0],
        next(w for w in check_case({"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0,
                                             "y1": 1, "nx": 2, "ny": 2},
                                    "physics": {"nu": 0.01},
                                    "boundary": {"left": {"type": "wall"}},
                                    "output": {"slice_axis": "x"}}) if "slice" in w),
    ]
    quoted = re.findall(r"^\| `([^`]*)`", table, re.M)
    for msg in real:
        hits = [q for q in quoted if all(part in msg for part in q.split("…") if part.strip())]
        assert hits, msg


def test_gui_3d_case_message():
    """T2 : à l'ouverture d'un cas 3D, « figures dans un plan z = constante » (périmé depuis
    les coupes x et y)."""
    pytest.importorskip("PySide6.QtWidgets", exc_type=ImportError)
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication

    from microrans.cli import examples_dir
    from microrans.gui.app import MainWindow
    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    try:
        w.open_case(examples_dir() / "conduite_carree_3d.toml")
        log = w.log_view.toPlainText()
        assert "Cas 3D : figures dans un plan x, y ou z = constante" in log
    finally:
        w.close()
        app.processEvents()
