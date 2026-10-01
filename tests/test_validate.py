"""Vérification du fichier de cas (fv2d/validate.py) : fautes de frappe signalées avec une
suggestion, aucun faux avertissement sur les exemples fournis (audit utilisateur C4)."""
import copy
import tomllib
from dataclasses import fields
from pathlib import Path

import pytest

from microrans.fv2d.validate import Key, Table, _solver_table, check_case

EX = Path(__file__).resolve().parent.parent / "microrans" / "examples"


def _ex(name):
    return tomllib.loads((EX / f"{name}.toml").read_text(encoding="utf-8"))


@pytest.mark.parametrize("path", sorted(EX.glob("*.toml")), ids=lambda p: p.stem)
def test_examples_have_no_warning(path):
    cfg = tomllib.loads(path.read_text(encoding="utf-8"))
    assert check_case(cfg) == []
    assert check_case(cfg, mesh_only=True) == []


def test_schema_covers_solver_settings_and_is_documented():
    from microrans.fv2d.compressible import CompressibleSettings
    from microrans.fv2d.solver import Settings
    solver = _solver_table().keys
    for cls in (Settings, CompressibleSettings):
        missing = {f.name for f in fields(cls)} - set(solver)
        assert not missing, f"{cls.__name__} : clés absentes de validate.py : {missing}"

    def undocumented(t: Table, path):
        for k, v in t.keys.items():
            if isinstance(v, Table):
                yield from undocumented(v, path + [k])
            elif not v.doc:
                yield ".".join(path + [k])
    from microrans.fv2d.validate import SCHEMA
    assert list(undocumented(Table({**SCHEMA.keys, "solver": Table(solver)}), [])) == []
    assert all(isinstance(v, (Key, Table)) for v in solver.values())


def _warn(cfg):
    return "\n".join(check_case(cfg))


def test_typos_are_reported_with_suggestion():
    cav = _ex("cavite_re100")
    c = copy.deepcopy(cav)
    c["solver"]["max_iters"] = c["solver"].pop("max_iter")
    assert "[solver] max_iters : clé inconnue, ignorée — vouliez-vous dire « max_iter » ?" \
        in _warn(c)
    c = copy.deepcopy(cav)
    c["solveur"] = c.pop("solver")
    assert "Section [solveur] inconnue, ignorée — vouliez-vous dire [solver] ?" in _warn(c)
    c = copy.deepcopy(cav)
    c["boundary"]["lid"]["vitesse"] = [1.0, 0.0]
    c["physics"]["Nu"] = c["physics"].pop("nu")
    w = _warn(c)
    assert "vitesse : clé inconnue, ignorée — vouliez-vous dire « U » ?" in w
    assert "Nu : clé inconnue, ignorée — vouliez-vous dire « nu » ?" in w
    c = copy.deepcopy(cav)
    c["nu"] = 0.01                                  # écrit avant toute section
    assert "à placer sous [physics]" in _warn(c)


def test_keys_without_effect_are_reported():
    c = _ex("compressible_rampe_mach2")
    c["physics"].update(nu=0.01, model="sst")
    w = _warn(c)
    assert "[physics] nu : sans effet avec le solveur compressible" in w
    assert "[flow] mu ou reynolds" in w and "pas de modèle de turbulence" in w
    c = _ex("cavite_re100")
    c["mesh"]["h_max"] = 0.1
    assert "[mesh] h_max : sans effet pour un maillage « rectangle »" in _warn(c)


def test_impossible_structure_is_an_error():
    c = _ex("cylindre_re20")
    c["bodies"] = c["bodies"][0]                    # [bodies] au lieu de [[bodies]]
    with pytest.raises(ValueError, match=r"écrire \[\[bodies\]\]"):
        check_case(c)
    c = _ex("cavite_re100")
    c["solver"] = 3
    with pytest.raises(ValueError, match=r"\[solver\] : table attendue"):
        check_case(c)


@pytest.mark.filterwarnings("default::microrans.fv2d.validate.CaseWarning")
def test_command_line_prints_warning(tmp_path, capsys):
    from microrans.cli import main
    from microrans.tomlio import dumps
    c = _ex("cavite_re100")
    c["mesh"].update(nx=6, ny=6)
    c["solver"].update(max_iter=2, relax_u=0.5)
    f = tmp_path / "cas.toml"
    f.write_text(dumps(c), encoding="utf-8")
    main(["run2d", str(f), "-o", str(tmp_path / "out"), "--no-plot", "-q"])
    err = capsys.readouterr().err
    assert "ATTENTION : [solver] relax_u : clé inconnue, ignorée — vouliez-vous dire " \
           "« relax_U » ?" in err
