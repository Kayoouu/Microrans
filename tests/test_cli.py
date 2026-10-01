"""Ligne de commande : --set à plusieurs niveaux, erreurs sans trace Python brute
(audit utilisateur M1, L3)."""
import pytest

from microrans.cli import _apply_set, main


def test_set_nested_lists_and_values():
    cfg = {"boundary": {"lid": {"type": "wall", "U": [1, 0]}},
           "bodies": [{"type": "circle", "radius": 0.5}]}
    _apply_set(cfg, ["boundary.lid.U=[2,0]", "bodies.0.radius=0.3", "physics.nu=0,01",
                     "physics.model=sst", "output.vtk=false"])
    assert cfg["boundary"]["lid"]["U"] == [2, 0]
    assert cfg["bodies"][0]["radius"] == 0.3
    assert cfg["physics"] == {"nu": 0.01, "model": "sst"}
    assert cfg["output"]["vtk"] is False


@pytest.mark.parametrize("item,expected", [
    ("solver.max_iter", "écrire SECTION.CLÉ=VALEUR"),
    ("max_iter=5", "clé « max_iter » incomplète"),
    ("bodies.3.radius=1", "indice 0 à 0 attendu"),
    ("boundary.lid.type.x=1", "« boundary.lid.type » est une valeur, pas une section"),
])
def test_set_errors(item, expected):
    cfg = {"boundary": {"lid": {"type": "wall"}}, "bodies": [{"type": "circle"}]}
    with pytest.raises(ValueError, match=expected):
        _apply_set(cfg, [item])


def _case(tmp_path, text):
    f = tmp_path / "cas.toml"
    f.write_text(text, encoding="utf-8")
    return str(f)


def test_errors_without_traceback(tmp_path, capsys):
    case = _case(tmp_path, '[mesh]\ntype = "file"\npath = "absent.msh"\n[physics]\nnu = 0.01\n'
                           '[boundary.wall]\ntype = "wall"\n')
    assert main(["run2d", case, "-o", str(tmp_path / "o"), "--no-plot", "-q"]) == 2
    err = capsys.readouterr().err
    assert "Erreur : [mesh] path : fichier de maillage introuvable" in err
    assert "Traceback" not in err
    assert main(["run2d", case, "--set", "solver.max_iter", "-q"]) == 2
    assert "écrire SECTION.CLÉ=VALEUR" in capsys.readouterr().err


def test_internal_error_message(monkeypatch, capsys):
    import microrans.cli as cli

    def boom(args):
        raise TypeError("défaut simulé")
    monkeypatch.setattr(cli, "cmd_examples", boom)
    assert main(["examples"]) == 3
    err = capsys.readouterr().err
    assert "Erreur interne inattendue (TypeError : défaut simulé)" in err
    assert "microrans --debug" in err and "Traceback" not in err
    assert main(["--debug", "examples"]) == 3
    assert "Traceback" in capsys.readouterr().err


def test_version(capsys):
    with pytest.raises(SystemExit):
        main(["--version"])
    assert capsys.readouterr().out.startswith("microrans ")


def test_examples_catalog_is_complete(capsys):
    """Audit U6 : chaque exemple fourni est classé (titre, durée) ; un nouvel exemple non
    catalogué apparaîtrait dans « Autres », ce test le signale."""
    from microrans.catalog import GROUPS, OTHERS, catalog
    from microrans.cli import examples_dir
    groups = catalog(examples_dir())
    assert OTHERS not in [g for g, _ in groups]
    listed = [stem for _, items in GROUPS for stem, *_ in items]
    assert len(listed) == len(set(listed)) == len(list(examples_dir().glob("*.toml")))
    assert all(stem.startswith("mesh_") for stem, *_ in GROUPS[-1][1])
    assert main(["examples"]) == 0
    out = capsys.readouterr().out
    assert "Commencer ici" in out and "cavite_re100" in out
