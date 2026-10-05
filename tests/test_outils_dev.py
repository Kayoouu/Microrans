"""Outils de développement (tools/dev) : contrôles avant poussée et comparaison de sorties."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools" / "dev"))

import ab  # noqa: E402
import avant_push  # noqa: E402
from _commun import BRANCHE, arbres_code  # noqa: E402


def test_push_command_rules():
    ok = f"git push -u origin {BRANCHE} 2>&1 | tail -1"
    assert avant_push.problemes_commande(ok) == []
    assert any("forcée" in p for p in avant_push.problemes_commande(
        f"git push --force origin {BRANCHE}"))
    assert any("forcée" in p for p in avant_push.problemes_commande(
        f"git push origin +{BRANCHE}"))
    assert any("tags" in p for p in avant_push.problemes_commande(
        f"git push origin {BRANCHE} --tags"))
    assert any("branche" in p for p in avant_push.problemes_commande("git push"))
    # poussée reconnue en position de commande seulement : pas le mot « push » d'un
    # message de commit, ni une poussée citée dans un texte ou un heredoc
    g = "git" + " push"
    for cmd, pousse in ((f"cd x && git -C /d push -u origin {BRANCHE}", True),
                        (f"for i in 1 2; do {g} -u origin {BRANCHE} && break; done", True),
                        ('git commit -m "push des résultats" && ls', False),
                        (f"echo 'avant `{g}` : branche'", False),
                        (f"python3 - <<'X'\nprint('{g}')\nX", False)):
        assert bool(avant_push._PUSH.search(cmd)) == pousse, cmd


def test_forbidden_added_lines():
    nom = "claude-" + "op" + "us-5"          # écrit en morceaux : pas de nom dans le dépôt
    diff = "\n".join([
        "+++ b/docs/a.md", "+modèle : " + nom, "-ligne retirée " + nom,
        "+++ b/tests/test_a.py", "+import tomllib", "+x = 1",
        "+++ b/microrans/b.py", "+import tomllib  # permis hors des tests",
        "+++ b/docs/c.md", "+Une fable de La Fontaine, " + "Op" + "us 5 cité",
    ])
    out = avant_push.lignes_interdites(diff)
    assert len(out) == 3
    assert out[0].startswith("docs/a.md") and "modèle" in out[0]
    assert out[1].startswith("tests/test_a.py") and "tomllib" in out[1]
    assert out[2].startswith("docs/c.md")


def test_compare_output_folders(tmp_path):
    a, b = tmp_path / "A", tmp_path / "B"
    for d, cd, t in ((a, 0.1, 1.0), (b, 0.1, 2.0)):
        d.mkdir()
        (d / "summary.json").write_text(json.dumps(
            {"wall_time_s": t, "airfoil": {"Cd": cd}, "out": str(d)}), encoding="utf-8")
        (d / "fig.png").write_bytes(b"x")
        np.savez(d / "checkpoint.npz", U=np.arange(3.0))
    (a / "history.csv").write_text("iteration,Ux\n1,0.5\n", encoding="utf-8")
    (b / "history.csv").write_text("iteration,Ux\n1,0.25\n", encoding="utf-8")
    res = ab.comparer_dossiers(a, b)
    assert res["summary.json"] == []          # temps et chemin du dossier ignorés
    assert res["checkpoint.npz"] == []
    assert res["fig.png"] is None
    assert res["history.csv"] == ["colonne Ux : écart max 0.25"]


def test_code_tree_ids():
    if not (Path(__file__).resolve().parents[1] / ".git").exists():
        pytest.skip("copie sans dépôt git (tools/dev/suite.py)")
    ids = arbres_code("HEAD")
    assert set(ids) == {"microrans", "tests"}
    assert all(len(v) == 40 for v in ids.values())
