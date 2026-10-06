"""Vérification du fichier de cas (fv2d/validate.py) : fautes de frappe signalées avec une
suggestion, aucun faux avertissement sur les exemples fournis (audit utilisateur C4)."""
import copy
from dataclasses import fields
from pathlib import Path

import pytest

from microrans.fv2d.validate import Key, Table, _solver_table, check_case
from microrans.tomlio import loads as toml_loads

EX = Path(__file__).resolve().parent.parent / "microrans" / "examples"


def _ex(name):
    return toml_loads((EX / f"{name}.toml").read_text(encoding="utf-8"))


@pytest.mark.parametrize("path", sorted(EX.glob("*.toml")), ids=lambda p: p.stem)
def test_examples_have_no_warning(path):
    cfg = toml_loads(path.read_text(encoding="utf-8"))
    assert check_case(cfg, mesh_only=True) == []
    if path.stem.startswith("mesh_"):               # audit M4 : message immédiat
        with pytest.raises(ValueError, match="seulement un maillage"):
            check_case(cfg)
    else:
        assert check_case(cfg) == []


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
    assert "vitesse : clé inconnue, ignorée — vouliez-vous dire « U » ?" in _warn(c)
    c["physics"]["Nu"] = c["physics"].pop("nu")     # plus de viscosité : erreur + la cause
    with pytest.raises(ValueError) as exc:
        check_case(c)
    assert "donner la viscosité nu" in str(exc.value)
    assert "Remarque : [physics] Nu : clé inconnue, ignorée — vouliez-vous dire « nu » ?" \
        in str(exc.value)
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


def _set(path, value):
    def f(c):
        d = c
        for p in path[:-1]:
            d = d.setdefault(p, {})
        d[path[-1]] = value
    return f


def _pop(path):
    def f(c):
        d = c
        for p in path[:-1]:
            d = d[p]
        d.pop(path[-1])
    return f


# (modification du cas de la cavité, extrait attendu du message) : campagne « entrées
# invalides » de l'audit (C5, M1, M2)
_BAD = [
    (_set(("physics", "nu"), 0.0), "[physics] nu = 0.0 : doit être > 0"),
    (_set(("physics", "nu"), -0.01), "[physics] nu = -0.01 : doit être > 0"),
    (_set(("mesh", "x1"), -1.0), "[mesh] x1 = -1.0 doit être > x0 = 0.0"),
    (_set(("mesh", "nx"), 0), "[mesh] nx = 0 : doit être ≥ 1"),
    (_set(("mesh", "nx"), 10.5), "[mesh] nx = 10.5 : nombre entier attendu."),
    (_set(("mesh", "nx"), 10.0), "[mesh] nx = 10.0 : nombre entier attendu (écrire 10, sans "
                                 "« .0 »)."),
    (_set(("mesh", "nx"), "dix"), "[mesh] nx = 'dix' : nombre entier attendu"),
    (_set(("solver", "max_iter"), "beaucoup"), "[solver] max_iter = 'beaucoup' : nombre entier"),
    (_set(("solver", "tol"), -1.0), "[solver] tol = -1.0 : doit être ≥ 0"),
    (_set(("solver", "relax_U"), 0.0), "[solver] relax_U = 0.0 : doit être > 0"),
    (_set(("solver", "relax_U"), 1.5), "[solver] relax_U = 1.5 : doit être ≤ 1"),
    (_set(("solver", "mode"), "unsteady"), "vouliez-vous dire « transient » ?"),
    (_set(("solver", "mode"), "transient"), "[solver] dt manquant (mode transient"),
    (lambda c: c["solver"].update(mode="transient", dt=0.0, t_end=1.0),
     "[solver] dt = 0.0 : doit être > 0"),
    (lambda c: c["solver"].update(mode="transient", dt=-0.01, t_end=1.0),
     "[solver] dt = -0.01 : doit être > 0"),
    (_set(("boundary", "lid", "U"), [1.0, 0.0, 0.0]), "[boundary.lid] U = [1.0, 0.0, 0.0] : "
                                                      "2 composantes attendues"),
    (_set(("boundary", "lid", "U"), "vite"), "[boundary.lid] U = 'vite' : 2 composantes"),
    (_set(("physics", "reynolds"), 5000.0), "nu = 0.01 et reynolds = 5000.0 donnés ensemble"),
    (_pop(("physics", "nu")), "[physics] : donner la viscosité nu (m²/s) ou le nombre de "
                              "Reynolds reynolds"),
    (_pop(("mesh",)), "[mesh] manquante"),
    (_set(("physics", "angle_of_attack"), "dix"), "[physics] angle_of_attack = 'dix' : nombre "
                                                  "attendu"),
    (_set(("solver", "algorithm"), "PISO"), "[solver] algorithm = 'PISO' inconnu. Choix : "
                                            "SIMPLE, SIMPLEC, coupled"),
    (_set(("solver", "convection_U"), "linearupwind"), "vouliez-vous dire « linearUpwind » ?"),
    (_set(("physics", "model_options"), {"ft3": True}), "option inconnue pour le modèle"),
    (_pop(("boundary", "lid", "type")), "[boundary.lid] type manquant"),
    (_set(("physics", "viscosity"), {"model": "carreau", "nu0": 1.0, "n": 0.5}),
     "[physics.viscosity] loi carreau : paramètre(s) nu_inf, lambda manquant(s), sans valeur "
     "par défaut (attendus : nu0, nu_inf, lambda, n)."),
    (_set(("physics", "viscosity"), {"model": "power_law", "K": 0.1, "n": -1.0}),
     "[physics.viscosity] n = -1.0 : doit être > 0"),
    (_set(("physics", "viscosity"), {"model": "carreaux"}), "vouliez-vous dire « carreau » ?"),
    (lambda c: (c["physics"].update(model="ke"), c["solver"].update(
        wall_treatment="wall_function")), "lois de paroi incompatibles avec le modèle ke"),
    (lambda c: (c["physics"].update(model="sst_gamma"), c["solver"].update(
        wall_treatment="wall_function")), "incompatibles avec le modèle sst_gamma"),
    (lambda c: c["physics"].update(model="sst", viscosity={
        "model": "power_law", "K": 0.1, "n": 0.5}), "non newtonien en laminaire uniquement"),
    (_set(("boundary", "lid", "omega"), "abc"), "[boundary.lid] omega = 'abc' : nombre attendu"),
    (_set(("boundary", "lid", "q"), "1,5"), "[boundary.lid] q = '1,5' : nombre attendu"),
]


@pytest.mark.parametrize("mod,expected", _BAD, ids=[e[:40] for _, e in _BAD])
def test_impossible_values_refused(mod, expected):
    c = _ex("cavite_re100")
    if "model_options" in expected or "option inconnue" in expected:
        c["physics"]["model"] = "sa"
    mod(c)
    with pytest.raises(ValueError) as exc:
        check_case(c)
    assert expected in str(exc.value)


def test_all_errors_reported_at_once():
    c = _ex("cavite_re100")
    c["physics"]["nu"] = -1.0
    c["solver"]["tol"] = -1.0
    with pytest.raises(ValueError) as exc:
        check_case(c)
    assert "nu = -1.0" in str(exc.value) and "tol = -1.0" in str(exc.value)


def test_compressible_values():
    c = _ex("compressible_rampe_mach2")
    c["flow"]["mach"] = -2.0
    with pytest.raises(ValueError, match=r"\[flow\] mach = -2.0 : doit être ≥ 0"):
        check_case(c)


def test_empty_file():
    with pytest.raises(ValueError, match="vide"):
        check_case({})


def test_condition_for_missing_patch():
    """Audit C5 : une condition pour une frontière absente du maillage était ignorée en
    silence ; une faute de frappe sur le nom donnait seulement « conditions manquantes »."""
    from microrans.fv2d.case import build_solver
    from microrans.fv2d.validate import CaseWarning
    c = _ex("cavite_re100")
    c["mesh"].update(nx=6, ny=6)
    c["boundary"]["inlet"] = {"type": "inlet", "U": [1.0, 0.0]}
    with pytest.warns(CaseWarning, match=r"\[boundary.inlet\] : aucune frontière « inlet » "
                                         r"dans le maillage \(frontières : "):
        build_solver(c)
    c = _ex("cavite_re100")
    c["mesh"].update(nx=6, ny=6)
    c["boundary"]["lidd"] = c["boundary"].pop("lid")
    with pytest.raises(ValueError, match=r"manquantes pour : lid — ajouter \[boundary.lid\] "
                                         r"avec type = un de : wall, inlet.* ; "
                                         r"\[boundary.lidd\] ne correspond à aucune frontière"
                                         r" — vouliez-vous dire « lid » \?"):
        build_solver(c)


def test_huge_structured_mesh_announced():
    """Audit M5 : 4 millions de cellules, aucun retour pendant des minutes."""
    c = _ex("cavite_re100")
    c["mesh"].update(nx=2000, ny=2000)
    w = check_case(c, mesh_only=True)
    assert len(w) == 1 and w[0].startswith("[mesh] 4 000 000 cellules : prévoir ~4.2 Go")
    c["mesh"].update(nx=200, ny=200)
    assert check_case(c) == []


def _two_bodies(kind, x2, r2=0.5, domain_x1=10.0):
    m = {"type": kind, "h_max": 2.0, "h_surface": 0.1}
    if kind == "hybrid":
        m["layers"] = {"n": 4, "first_height": 0.01, "ratio": 1.2}
    return {"mesh": m,
            "domain": {"x0": -5.0, "x1": domain_x1, "y0": -5.0, "y1": 5.0},
            "bodies": [{"type": "circle", "center": [0.0, 0.0], "radius": 0.5, "name": "b1"},
                       {"type": "circle", "center": [x2, 0.0], "radius": r2, "name": "b2"}]}


def test_overlapping_bodies_detected_before_meshing():
    """Audit U5 (mesuré) : corps identiques en hybride → « non manifold » après 6 à 20 s ;
    corps imbriqués ou trop proches (écart < 2 × épaisseur des couches, ici 0.054) →
    maillage « réussi » mais faux (couches dans un solide, frontières _b1_top) ; non
    structuré : corps recouvert disparu en silence."""
    with pytest.raises(ValueError, match="« b1 » et « b2 » se recouvrent ou se touchent : "
                                         "impossible en maillage hybride"):
        check_case(_two_bodies("hybrid", 0.0), mesh_only=True)
    with pytest.raises(ValueError, match="se recouvrent"):
        check_case(_two_bodies("hybrid", 0.1, 0.2), mesh_only=True)     # b2 dans b1
    with pytest.raises(ValueError, match=r"écart 0.08 < 2 × épaisseur des couches de paroi "
                                         r"\(0.0537\) : couches en collision"):
        check_case(_two_bodies("hybrid", 1.08), mesh_only=True)
    assert check_case(_two_bodies("hybrid", 1.12), mesh_only=True) == []
    w = check_case(_two_bodies("unstructured", 0.0), mesh_only=True)
    assert len(w) == 1 and "ils formeront un seul obstacle" in w[0]
    assert check_case(_two_bodies("unstructured", 1.08), mesh_only=True) == []


def test_body_at_domain_edge():
    with pytest.raises(ValueError, match="« b2 » : à 0 du bord du domaine"):
        check_case(_two_bodies("hybrid", 9.8), mesh_only=True)
    assert check_case(_two_bodies("unstructured", 9.8), mesh_only=True) == []  # coupé : permis
    w = check_case(_two_bodies("unstructured", 20.0), mesh_only=True)
    assert w == ["[[bodies]] « b2 » est entièrement hors du domaine de calcul [domain] : il "
                 "sera ignoré."]


def test_wall_function_in_laminar_warned():
    c = _ex("cavite_re100")
    c["solver"]["wall_treatment"] = "wall_function"
    assert check_case(c) == ["[solver] wall_treatment = \"wall_function\" : loi de paroi "
                             "turbulente, sans objet en laminaire (garder \"resolved\")."]


def test_reference_document_up_to_date():
    """Audit D4 : docs/reference_cas.md est généré à partir de SCHEMA ; il doit être
    régénéré (python -m microrans.fv2d.validate > docs/reference_cas.md) quand une clé
    change, sinon la documentation ment."""
    from microrans.fv2d.validate import reference_markdown
    doc = Path(__file__).resolve().parent.parent / "docs" / "reference_cas.md"
    assert doc.read_text(encoding="utf-8") == reference_markdown()
    text = reference_markdown()
    for key in ("sst_gamma", "linearUpwindLimited", "`max_iter` | itérations maximales "
                "(stationnaire) | 3000 | 5000"):
        assert key in text
