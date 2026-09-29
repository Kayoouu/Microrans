"""Sauvegarde / reprise : reprise exacte sur le même maillage, interpolation sinon."""
import csv
import json

import numpy as np
import pytest

from microrans.cli import main
from microrans.fv2d import Settings, Solver2D
from microrans.fv2d.restart import load_checkpoint, save_checkpoint
from microrans.mesh2d import rectangle_mesh


def _cavity(n=16, scheme="backward", model="laminar", energy=None, adjust=False):
    m = rectangle_mesh(0, 1, 0, 1, n, n, names={"left": "walls", "right": "walls",
                                                "bottom": "walls", "top": "lid"})
    bc = {"lid": {"type": "wall", "U": [1.0, 0.0]}, "walls": {"type": "wall"}}
    if energy:
        bc["lid"]["T"], bc["walls"]["T"] = 1.0, 0.0
    return Solver2D(m, 0.01 if model == "laminar" else 1e-4, bc, model=model, energy=energy,
                    settings=Settings(time_scheme=scheme, adjust_dt=adjust, max_co=0.8),
                    reference_velocity=1.0)


@pytest.mark.parametrize("scheme,model,energy,adjust", [
    ("backward", "sa", {"Pr": 0.7}, True),      # BDF2 à pas variable + turbulence + thermique
    ("rk3", "laminar", None, True),             # explicite, pas adaptatif
    ("ab2", "laminar", None, False),            # historique du second membre
])
def test_transient_restart_is_bitwise_exact(tmp_path, scheme, model, energy, adjust):
    """Arrêt au milieu d'un calcul puis reprise : résultat identique au bit près à un calcul
    sans interruption (niveaux de temps précédents et Δt sauvegardés)."""
    f = tmp_path / "c.npz"

    def cb(s, n):
        if n == 5:
            save_checkpoint(s, f)

    a = _cavity(scheme=scheme, model=model, energy=energy, adjust=adjust)
    a.run_transient(0.02, 0.3, callback=cb)
    b = _cavity(scheme=scheme, model=model, energy=energy, adjust=adjust)
    info = load_checkpoint(b, f)
    assert info["mode"] == "exact" and 0 < b.time < 0.3
    b.run_transient(0.02, 0.3)
    assert b.time == a.time
    assert np.array_equal(a.U, b.U) and np.array_equal(a.p, b.p)
    for k in a.state:
        assert np.array_equal(a.state[k], b.state[k])
    if energy:
        assert np.array_equal(a.T, b.T)


def test_steady_restart_exact_and_history_continues(tmp_path):
    a = _cavity(24)
    a.run_steady(max_iter=200, tol=1e-12)
    b = _cavity(24)
    b.run_steady(max_iter=80, tol=1e-12)
    save_checkpoint(b, tmp_path / "c.npz")
    c = _cavity(24)
    load_checkpoint(c, tmp_path / "c.npz")
    c.run_steady(max_iter=a.iterations - 80, tol=1e-12)
    assert np.array_equal(a.U, c.U)
    assert c.iterations_total == a.iterations_total
    assert [h["iteration"] for h in c.history] == list(range(1, a.iterations + 1))


def test_interpolation_onto_finer_mesh_speeds_up_convergence(tmp_path):
    """Démarrage d'un maillage fin depuis un calcul grossier (mapFields) : même solution,
    moins d'itérations (mesuré : 128² depuis 32², 508 itérations au lieu de 1135)."""
    g = _cavity(16)
    g.run_steady(max_iter=2000, tol=1e-6)
    save_checkpoint(g, tmp_path / "c.npz")
    f0 = _cavity(40)
    f0.run_steady(max_iter=3000, tol=1e-6)
    f1 = _cavity(40)
    assert load_checkpoint(f1, tmp_path / "c.npz")["mode"] == "interpolé"
    f1.run_steady(max_iter=3000, tol=1e-6)
    assert f1.iterations < 0.85 * f0.iterations          # 129 contre 167 ici (40²)
    assert np.abs(f1.U - f0.U).max() < 1e-4


def test_other_model_keeps_freestream_for_missing_variables(tmp_path):
    lam = _cavity(12)
    lam.run_steady(max_iter=50)
    save_checkpoint(lam, tmp_path / "c.npz")
    sa = _cavity(12, model="sa")
    info = load_checkpoint(sa, tmp_path / "c.npz")
    assert info["ignored"] == ["nu_tilde"] or info["ignored"] == list(sa.state)
    assert np.array_equal(sa.U, lam.U)


CASE = """
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
nu = 0.01
reference_velocity = 1.0
[boundary.lid]
type = "wall"
U = [1.0, 0.0]
[boundary.walls]
type = "wall"
[solver]
mode = "%s"
max_iter = 40
tol = 1e-12
dt = 0.05
t_end = 0.5
[output]
plots = false
"""


def test_cli_continue_steady_and_transient(tmp_path):
    for mode in ("steady", "transient"):
        case = tmp_path / f"{mode}.toml"
        case.write_text(CASE % mode)
        out = tmp_path / mode
        main(["run2d", str(case), "-o", str(out), "-q"])
        assert (out / "checkpoint.npz").exists()
        if mode == "steady":
            main(["run2d", str(case), "-o", str(out), "-q", "--continue"])
            s = json.loads((out / "summary.json").read_text())
            assert s["iterations"] == 80 and s["restart"]["mode"] == "exact"
        else:
            main(["run2d", str(case), "-o", str(out), "-q", "--continue",
                  "--set", "solver.t_end=1.0"])
            with open(out / "history.csv", encoding="utf-8") as fh:
                t = [float(r["time"]) for r in csv.DictReader(fh)]
            assert len(t) == 20 and t[-1] == pytest.approx(1.0) and np.all(np.diff(t) > 0)
