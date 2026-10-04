"""Incidence de l'écoulement, polaire et balayage de paramètres."""
import copy
import csv

import numpy as np
import pytest

from microrans.cli import main
from microrans.fv2d.case import run_case
from microrans.fv2d.sweep import parse_values, run_sweep


def test_parse_values():
    assert parse_values("-4:12:4") == [-4.0, 0.0, 4.0, 8.0, 12.0]
    assert parse_values("0:1:0.25") == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert parse_values("10, 20; 40") == [10.0, 20.0, 40.0]
    assert parse_values([1, 2]) == [1.0, 2.0]
    with pytest.raises(ValueError):
        parse_values("0:10:-1")


def _ogrid_case(body, re, n_around, n_radial, first, radius=20.0):
    return {"mesh": {"type": "ogrid", "n_around": n_around, "n_radial": n_radial,
                     "farfield_radius": radius, "first_height": first},
            "bodies": [body],
            "physics": {"reynolds": re, "reference_velocity": 1.0, "reference_length": 1.0},
            "initial": {"U": [1.0, 0.0]},
            "boundary": {body["name"]: {"type": "wall"},
                         "farfield": {"type": "farfield", "U": [1.0, 0.0]}},
            "solver": {"max_iter": 3000, "tol": 1e-6},
            "output": {"plots": False, "vtk": False}}


def test_incidence_rotates_flow_and_force_axes(tmp_path):
    """Cylindre (symétrie de révolution) : Cd dans les axes de l'écoulement indépendant de
    l'incidence, portance nulle."""
    cfg = _ogrid_case({"type": "circle", "center": [0, 0], "radius": 0.5, "name": "cyl"},
                      20, 48, 32, 0.02)
    res = {}
    for a in (0.0, 30.0):
        c = copy.deepcopy(cfg)
        c["physics"]["angle_of_attack"] = a
        res[a] = run_case(c, out_dir=tmp_path / str(a), verbose=False, plot=False)["cyl"]
    assert res[30.0]["Cd"] == pytest.approx(res[0.0]["Cd"], rel=1e-3)
    assert abs(res[30.0]["Cl"]) < 2e-3


def test_symmetric_airfoil_polar(tmp_path):
    """Profil symétrique, maillage symétrique : Cl(−α) = −Cl(α), Cd(−α) = Cd(α),
    Cm(−α) = −Cm(α) ; la continuation ne change pas la solution (hors décrochage)."""
    cfg = _ogrid_case({"type": "naca", "code": "0012", "chord": 1.0, "name": "airfoil"},
                      500, 48, 24, 5e-3, radius=10.0)
    cfg["output"]["moment_center"] = [0.25, 0.0]
    cfg["solver"]["tol"] = 1e-5
    rows = run_sweep(cfg, "physics.angle_of_attack", "-3:3:3", out_dir=tmp_path / "c",
                     verbose=False)
    cl, cd, cm = ([r[f"{k}_airfoil"] for r in rows] for k in ("Cl", "Cd", "Cm"))
    assert all(r["converged"] for r in rows)
    assert cl[0] == pytest.approx(-cl[2], rel=1e-3) and cl[2] > 0.05
    assert abs(cl[1]) < 1e-5 * cl[2] + 1e-8
    assert cd[0] == pytest.approx(cd[2], rel=1e-3) and cd[1] < cd[2]
    assert cm[0] == pytest.approx(-cm[2], rel=1e-2, abs=1e-6)
    assert (tmp_path / "c" / "polaire.png").exists()
    cold = run_sweep(cfg, "physics.angle_of_attack", "0:3:3", out_dir=tmp_path / "f",
                     continuation=False, verbose=False, plot=False)
    assert rows[2]["Cl_airfoil"] == pytest.approx(cold[1]["Cl_airfoil"], rel=1e-3)


def test_cli_sweep_reynolds(tmp_path):
    case = tmp_path / "cav.toml"
    case.write_text("""
[mesh]
type = "rectangle"
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
nx = 12
ny = 12
names = { left = "walls", right = "walls", bottom = "walls", top = "lid" }
[physics]
reynolds = 100
reference_velocity = 1.0
[boundary.lid]
type = "wall"
U = [1.0, 0.0]
[boundary.walls]
type = "wall"
[solver]
max_iter = 2000
tol = 1e-6
""")
    out = tmp_path / "out"
    main(["sweep", str(case), "-o", str(out), "-q", "--param", "physics.reynolds",
          "--values", "10", "100", "--jobs", "2"])
    with open(out / "balayage.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert [float(r["physics.reynolds"]) for r in rows] == [10.0, 100.0]
    assert all(r["converged"] == "True" for r in rows)
    # frottement sur le couvercle : décroît avec Re (Cd = force / (½ U² L))
    cd = [abs(float(r["Cd_lid"])) for r in rows]
    assert np.all(np.isfinite(cd)) and cd[0] > 5 * cd[1]
    assert (out / "balayage.png").exists()


def test_reynolds_sweep_on_case_given_in_nu(tmp_path):
    """Audit : sur un cas donné en nu (prioritaire), balayer physics.reynolds calculait
    tous les points avec le même ν (résultats identiques, sans message)."""
    cfg = _small_cavity()
    cfg["solver"]["max_iter"] = 300
    rows = run_sweep(cfg, "physics.reynolds", [10, 100], out_dir=tmp_path, verbose=False,
                     plot=False, continuation=False)
    cd = [abs(r["Cd_lid"]) for r in rows]
    assert cd[0] > 5 * cd[1]                        # frottement ∝ 1/Re environ
    c = {**cfg, "physics": {"reynolds": 100.0, "reference_velocity": 1.0}}
    rows = run_sweep(c, "physics.nu", [0.1, 0.01], out_dir=tmp_path / "nu", verbose=False,
                     plot=False, continuation=False)
    assert abs(rows[0]["Cd_lid"]) > 5 * abs(rows[1]["Cd_lid"])


def _small_cavity():
    return {"mesh": {"type": "rectangle", "x0": 0.0, "x1": 1.0, "y0": 0.0, "y1": 1.0,
                     "nx": 12, "ny": 12,
                     "names": {"left": "walls", "right": "walls", "bottom": "walls",
                               "top": "lid"}},
            "physics": {"nu": 0.01, "reference_velocity": 1.0},
            "boundary": {"lid": {"type": "wall", "U": [1.0, 0.0]}, "walls": {"type": "wall"}},
            "solver": {"max_iter": 2000, "tol": 1e-6},
            "output": {"plots": False, "vtk": False}}


def test_parallel_sweep_equals_serial(tmp_path):
    """jobs = 2 (processus « spawn ») : sans continuation, résultats identiques au bit près
    au calcul séquentiel ; avec continuation (2 blocs contigus), identiques au calcul
    séquentiel de chaque bloc et égaux au balayage séquentiel à la tolérance près."""
    cfg, vals = _small_cavity(), "0.01, 0.02, 0.04"
    seen = []
    ser = run_sweep(cfg, "physics.nu", vals, out_dir=tmp_path / "s", continuation=False,
                    verbose=False, plot=False)
    par = run_sweep(cfg, "physics.nu", vals, out_dir=tmp_path / "p", continuation=False,
                    verbose=False, plot=False, jobs=2, on_point=seen.append)
    assert par == ser                                  # dans l'ordre des valeurs
    assert sorted(r["physics.nu"] for r in seen) == [0.01, 0.02, 0.04]
    assert (tmp_path / "p" / "balayage.csv").exists()
    cser = run_sweep(cfg, "physics.nu", vals, out_dir=tmp_path / "cs", verbose=False,
                     plot=False)
    cpar = run_sweep(cfg, "physics.nu", vals, out_dir=tmp_path / "cp", verbose=True,
                     plot=False, jobs=2)
    # blocs [0.01, 0.02] et [0.04] : ce dernier part de l'état initial
    assert cpar[:2] == cser[:2] and cpar[2] == ser[2]
    assert all(r["converged"] for r in cpar)
    assert cpar[2]["Cd_lid"] == pytest.approx(cser[2]["Cd_lid"], rel=1e-4)
    assert (tmp_path / "cp" / "nu_p0_04" / "journal.txt").is_file()   # journal par point


def test_parallel_sweep_stop(tmp_path):
    """Arrêt demandé : plus aucun point lancé, les points en cours s'arrêtent."""
    rows = run_sweep(_small_cavity(), "physics.nu", "0.01, 0.02, 0.04, 0.08",
                     out_dir=tmp_path, continuation=False, verbose=False, plot=False, jobs=2,
                     should_stop=lambda: True)
    assert len(rows) <= 2 and not any(r["converged"] for r in rows)


def test_pitching_moment_sign_nose_up_positive(tmp_path):
    """Convention aéronautique : autour du bord d'attaque, la portance (appliquée vers le
    quart de corde) pique le profil : C_m,BA ≈ −C_l/4 < 0."""
    cfg = _ogrid_case({"type": "naca", "code": "0012", "chord": 1.0, "name": "airfoil"},
                      500, 48, 24, 5e-3, radius=10.0)
    cfg["physics"]["angle_of_attack"] = 4.0
    cfg["output"]["moment_center"] = [0.0, 0.0]
    cfg["solver"].update(tol=1e-5, algorithm="coupled")
    s = run_case(cfg, out_dir=tmp_path, verbose=False, plot=False)
    cl, cm = s["airfoil"]["Cl"], s["airfoil"]["Cm"]
    assert cl > 0.1
    assert cm == pytest.approx(-cl / 4, rel=0.15)


@pytest.mark.parametrize("text,expected", [
    ("0 à 4", "« 0 à 4 » n'est pas un nombre. Écrire début:fin:pas"),
    ("0:4", "trois nombres attendus, début:fin:pas"),
    ("0:10:-1", "le pas -1 ne mène pas de 0 à 10"),
    ("0:4:0", "pas nul"),
    ("", "Aucune valeur"),
])
def test_parse_values_errors_explained(text, expected):
    """Audit U10 : « could not convert string to float » en anglais."""
    with pytest.raises(ValueError, match=expected):
        parse_values(text)


def test_parse_values_decimal_comma():
    assert parse_values("0,5; 1; 1,5") == [0.5, 1.0, 1.5]
    assert parse_values("0:1:0,5") == [0.0, 0.5, 1.0]


_CAV = """
[mesh]
type = "rectangle"
x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 1.0
nx = 8
ny = 8
names = { left = "walls", right = "walls", bottom = "walls", top = "lid" }
[physics]
reynolds = 100
[boundary.lid]
type = "wall"
U = [1.0, 0.0]
[boundary.walls]
type = "wall"
[solver]
max_iter = 15
"""


def test_sweep_run_outside_case_directory(tmp_path, monkeypatch):
    """C11 (audit 2) : lancé depuis un autre dossier que celui du cas, avec une sortie
    relative, le 2e point repartait d'un checkpoint introuvable (« Fichier de reprise
    introuvable ») ; même chose pour le démarrage multigrille."""
    (tmp_path / "cas").mkdir()
    (tmp_path / "cas" / "cav.toml").write_text(_CAV)
    (tmp_path / "ailleurs").mkdir()
    monkeypatch.chdir(tmp_path / "ailleurs")
    rc = main(["sweep", "../cas/cav.toml", "-o", "sw", "-q", "--no-plot", "--param",
               "physics.reynolds", "--values", "50", "100", "200"])
    assert rc in (0, 1)                                 # 1 = non convergé (15 itérations)
    with open("sw/balayage.csv", encoding="utf-8") as fh:
        assert len(list(csv.DictReader(fh))) == 3
    rc = main(["run2d", "../cas/cav.toml", "-o", "fmg", "-q", "--no-plot", "--set",
               "solver.fmg_levels=1"])
    assert rc in (0, 1) and (tmp_path / "ailleurs" / "fmg" / "fmg_niveau1.npz").is_file()


@pytest.mark.parametrize("key, expected", [
    ("physics.reynold", "vouliez-vous dire « reynolds »"),
    ("solver.max_iterr", "vouliez-vous dire « max_iter »"),
    ("physique.reynolds", "[physics]"),
    ("solver.cfl", "sans effet avec le solveur incompressible"),
])
def test_sweep_refuses_unknown_or_useless_key(tmp_path, key, expected):
    """C21 (audit 2 approfondi) : une clé mal orthographiée était balayée en silence
    (points identiques étiquetés de valeurs différentes)."""
    from microrans.tomlio import loads
    with pytest.raises(ValueError, match="Balayage") as exc:
        run_sweep(loads(_CAV), key, [1, 2], out_dir=tmp_path, verbose=False, plot=False)
    assert expected in str(exc.value)
    assert not (tmp_path / "balayage.csv").exists()


def test_mesh_resolution_sweep(tmp_path):
    """C20 (audit 2 approfondi) : --param mesh.nx --values 8 16 (étude de convergence en
    maillage) → « Erreur interne (TypeError) » ; les valeurs lues sont des réels."""
    from microrans.tomlio import loads
    rows = run_sweep(loads(_CAV), "mesh.nx", [6, 10], out_dir=tmp_path, verbose=False,
                     plot=False)
    assert [r["mesh.nx"] for r in rows] == [6.0, 10.0]
    assert rows[0]["Cd_lid"] != rows[1]["Cd_lid"]
