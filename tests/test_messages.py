"""Lot F2 (audit 2) : messages en français, avec les choix possibles et une suggestion, à la
place de messages Python bruts ou d'erreurs internes. Un test par point d'audit."""
import copy
import os
import warnings

import pytest

from microrans.cli import examples_dir
from microrans.mesh2d.builder import load_config


def _ex(name):
    return copy.deepcopy(load_config(examples_dir() / f"{name}.toml"))


def _cavity(n=4):
    c = _ex("cavite_re100")
    c["mesh"].update(nx=n, ny=n)
    c["solver"]["max_iter"] = 2
    return c


def test_periodic_with_unknown_patch_names_it():
    """M13 : avant « 'lidd' is not in list »."""
    from microrans.mesh2d.builder import build_mesh
    c = _cavity()
    c["mesh"]["periodic"] = [["walls", "lidd"]]
    with pytest.raises(ValueError, match=r"frontière « lidd » inexistante\. Frontières du "
                                         r"maillage : walls, lid — vouliez-vous dire « lid »"):
        build_mesh(c)
    c = _ex("conduite_carree_3d")                   # même message en 3D
    c["mesh"]["periodic"] = [["inlet", "outlett"]]
    with pytest.raises(ValueError, match=r"« outlett » inexistante.*« outlet » \?"):
        build_mesh(c)


def test_unknown_patch_type_lists_choices():
    """M15 : avant « Type de patch inconnu 'symetrie' (z0) », sans choix ni suggestion."""
    from microrans.mesh2d.builder import build_mesh
    c = _cavity()
    c["mesh"]["patch_types"] = {"walls": "symetrie"}
    with pytest.raises(ValueError, match=r"« symetrie » \(frontière walls\)\. Choix : wall, "
                                         r"patch, symmetry, empty — vouliez-vous dire "
                                         r"« symmetry » \?"):
        build_mesh(c)


def test_key_in_wrong_section_names_the_right_one():
    """M16 : avant « clé inconnue, ignorée » sans dire où la mettre."""
    from microrans.fv2d.validate import check_case
    c = _cavity()
    c["physics"]["moment_center"] = [0.25, 0.0]
    c["physics"]["max_iter"] = 10
    w = check_case(c)
    assert "[physics] moment_center : clé inconnue, ignorée — clé de [output] : la " \
           "déplacer" in w
    assert any("max_iter" in x and "[solver]" in x for x in w)


def test_missing_boundary_condition_message():
    """M17 : avant « … pour les patches ['front'] » (liste Python) ; faces d'extrusion :
    la condition d'un écoulement 2D est proposée."""
    from microrans.fv2d.case import build_solver
    c = _cavity()
    c["mesh"]["extrude"] = {"z0": 0.0, "z1": 0.5, "nz": 1}
    c["boundary"]["lid"]["U"] = [1.0, 0.0, 0.0]
    with pytest.raises(ValueError) as e:
        build_solver(c)
    msg = str(e.value)
    assert "manquantes pour : back, front — ajouter [boundary.back] et [boundary.front]" in msg
    assert "type = \"symmetry\"" in msg and "['" not in msg


def test_expert_solver_keys_checked_before_the_run():
    """M14 / M18 : texte → « Erreur interne inattendue » ; −1 accepté en silence ;
    cn_theta = −1 → divergence au lieu d'un refus. Tout est signalé ensemble."""
    from microrans.fv2d.validate import check_case
    c = _ex("compressible_rampe_mach2")
    c["solver"].update(venkat_k="0.3", linear_iter="x", limiter_freeze=-1, cfl_growth=0.5,
                       entropy_fix=-1, linear_tol=2.0, implicit_jacobian="roee", order=3)
    with pytest.raises(ValueError) as e:
        check_case(c)
    msg = str(e.value)
    for part in ("venkat_k = '0.3' : nombre attendu", "linear_iter = 'x' : nombre entier",
                 "limiter_freeze = -1 : doit être ≥ 0", "cfl_growth = 0.5 : doit être ≥ 1",
                 "entropy_fix = -1", "linear_tol = 2.0 : doit être ≤ 1",
                 "implicit_jacobian = 'roee' inconnu", "order = 3 : 1 ou 2"):
        assert part in msg, part
    c = _cavity()
    c["solver"].update(cn_theta=-1, nonorth_limit=1.5, threads="4", numba="oui", dt="abc")
    with pytest.raises(ValueError) as e:
        check_case(c)
    msg = str(e.value)
    for part in ("cn_theta = -1 : entre 0.5", "nonorth_limit = 1.5 : doit être ≤ 1",
                 "threads = '4' : nombre entier", "numba = 'oui' : true ou false",
                 "dt = 'abc' : nombre attendu"):          # stationnaire : dt non lu
        assert part in msg, part
    for name in sorted(os.listdir(examples_dir())):   # aucun exemple refusé à tort
        if name.endswith(".toml") and not name.startswith("mesh_"):
            check_case(load_config(examples_dir() / name))


def test_minor_messages():
    """M19 : trailing_edge ignoré en silence sur un cercle ; noms internes dans le message
    de --steps-per-period."""
    from microrans.cases import run_pulsating_channel
    from microrans.fv2d.validate import check_case
    c = _ex("cylindre_re20")
    c["bodies"][0]["trailing_edge"] = "open"
    assert any("trailing_edge sans effet (type = \"circle\"" in w for w in check_case(c))
    with pytest.raises(ValueError, match=r"--steps-per-period = 50 : doit être un multiple "
                                         r"de 8 .* 48 ou 56"):
        run_pulsating_channel(steps_per_period=50)


def test_case_file_encodings_and_syntax(tmp_path):
    """M20 : BOM (Bloc-notes, PowerShell) → « Invalid statement (at line 1, column 1) » ;
    Latin-1 → message du codec ; syntaxe → message anglais seul."""
    text = "[physics]\n# viscosité\nnu = 0.01\n"
    f = tmp_path / "bom.toml"
    f.write_bytes(b"\xef\xbb\xbf" + text.encode("utf-8"))
    with warnings.catch_warnings():
        warnings.simplefilter("error")             # pas de remarque pour un BOM
        assert load_config(f) == {"physics": {"nu": 0.01}}
    f = tmp_path / "latin1.toml"
    f.write_bytes(text.encode("latin-1"))
    with pytest.warns(UserWarning, match="n'est pas en UTF-8"):
        assert load_config(f) == {"physics": {"nu": 0.01}}
    f = tmp_path / "faux.toml"
    f.write_text("[physics]\nnu = 0,01\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"faux\.toml : syntaxe TOML incorrecte \(ligne 2, "
                                         r"colonne \d+\).*point décimal"):
        load_config(f)


def test_output_folder_is_a_file_and_binary_gmsh(tmp_path, capsys):
    """M20 : -o vers un fichier → « [Errno 17] File exists » ; Gmsh binaire → erreur
    interne (KeyError 'Nodes')."""
    from microrans.cli import main
    from microrans.mesh2d.io import read_gmsh
    (tmp_path / "f").write_text("x", encoding="utf-8")
    rc = main(["run2d", "cavite_re100", "-o", str(tmp_path / "f"), "--no-plot", "-q",
               "--set", "mesh.nx=4", "mesh.ny=4", "solver.max_iter=1"])
    assert rc == 2
    assert "existe déjà et n'est pas un dossier : choisir un autre dossier de sortie" \
        in capsys.readouterr().err
    msh = tmp_path / "b.msh"
    msh.write_bytes(b"$MeshFormat\n4.1 1 8\n\x01\x00\x00\x00\n$EndMeshFormat\n")
    with pytest.raises(ValueError, match="fichier Gmsh binaire, non lu. Le réexporter en "
                                         "texte"):
        read_gmsh(msh)


def test_single_cell_mesh_refused_clearly():
    """M21 : avant « Factor is exactly singular »."""
    from microrans.fv2d.case import build_solver
    with pytest.raises(ValueError, match="1 cellule\\(s\\) sans face intérieure"):
        build_solver(_cavity(1))


def test_folder_with_an_example_name(tmp_path, monkeypatch):
    """M22 : un dossier de résultats nommé comme un exemple le masquait (« un fichier est
    attendu, pas un dossier »)."""
    from microrans.cli import resolve_example
    monkeypatch.chdir(tmp_path)
    (tmp_path / "cavite_re100").mkdir()
    with pytest.warns(UserWarning, match="est aussi un dossier ici : exemple "
                                         "cavite_re100.toml utilisé"):
        assert resolve_example("cavite_re100") == examples_dir() / "cavite_re100.toml"
    (tmp_path / "pas_un_exemple").mkdir()
    with pytest.raises(ValueError, match="est un dossier, pas un fichier de cas"):
        resolve_example("pas_un_exemple")


def test_probes_with_decimal_comma():
    """U16 : « 0,25 0,75 » lu comme 4 nombres (refus au lancement seulement)."""
    from microrans.fv2d.sampling import parse_points
    from microrans.fv2d.validate import check_case
    assert parse_points("0,25 0,75 ; 1 0,5").tolist() == [[0.25, 0.75], [1.0, 0.5]]
    assert parse_points("0.25, 0.75").tolist() == [[0.25, 0.75]]
    c = _cavity()
    c["output"]["probes"] = "0,25,0,75"
    with pytest.raises(ValueError, match=r"\[output\] probes = '0,25,0,75' : Points "
                                         r"invalides"):
        check_case(c)


def test_mean_velocity_without_rounding_noise():
    """L6 : « (17.64, -4.447e-16, -2.495e-39) » affiché dans le résumé 3D."""
    from microrans.fv2d.report import summary_text
    txt = summary_text({"converged": True, "iterations": 3, "dimension": 3,
                        "U_mean": [15.35, 4.759e-17, -6.668e-47]})
    assert "Vitesse moyenne dans le domaine : (15.35, 0, 0)." in txt
