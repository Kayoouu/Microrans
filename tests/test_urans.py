"""URANS : cohérence avec le RANS, bilan moyen, périodicité, interface CLI."""
import json

import numpy as np
import pytest

from microrans import channel_grid, get_model, solve_steady, solve_unsteady
from microrans.cases import PulsatingForcing, run_pulsating_channel
from microrans.cli import main
from microrans.models import TURBULENT_MODELS


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_urans_without_forcing_stays_on_rans_solution(model):
    re_tau = 395.0
    nu = 1.0 / re_tau
    g = channel_grid(96, re_tau, 0.5)
    m = get_model(model, g, nu)
    steady = solve_steady(m, g, nu)
    assert steady.converged
    res = solve_unsteady(m, g, nu, PulsatingForcing(1.0, 0.0, 1.0), steady,
                         t_end=0.5, dt=0.01)
    assert np.max(np.abs(res.final.U - steady.U)) / np.max(steady.U) < 1e-7
    for k in m.variables:
        rel = np.max(np.abs(res.final.state[k] - steady.state[k])) / np.max(steady.state[k])
        assert rel < 1e-6, k


@pytest.mark.parametrize("model", ["sa", "sst"])
def test_pulsating_channel_periodic_regime(model):
    r = run_pulsating_channel(model, re_tau=180.0, n_cells=64, y1_plus=0.5,
                              omega_plus=0.02, amplitude=5.0, steps_per_period=64,
                              t_transient=60.0, n_average=3)
    s = r.summary
    # Bilan de quantité de mouvement moyenné sur une période : ⟨τ_w⟩ = ⟨f⟩ h = 1.
    assert s["mean_tau_wall"] == pytest.approx(1.0, abs=2e-3)
    assert s["periodicity_error"] < 1e-3
    # Cœur en écoulement bouchon : amplitude de U_b ≈ A/ω à haute fréquence (à ~15 % près).
    assert s["Ub_amplitude"] == pytest.approx(5.0 / r.forcing.omega, rel=0.15)
    assert s["max_inner_iterations"] < 30
    assert r.phase_U.shape == (8, r.grid.n)


def test_cli_rans_and_urans(tmp_path):
    out = tmp_path / "rans"
    assert main(["rans", "-m", "sa", "kw", "--n-cells", "96", "--y1plus", "0.5",
                 "-o", str(out), "-q"]) == 0
    for m in ("sa", "kw"):
        s = json.loads((out / m / "summary.json").read_text(encoding="utf-8"))
        assert s["converged"]
        assert (out / m / "profiles.csv").exists()
        assert (out / m / f"rans_{m}.png").exists()
    assert (out / "comparison.png").exists()
    assert (out / "summary.csv").exists()

    out = tmp_path / "urans"
    assert main(["urans", "-m", "sa", "--re-tau", "180", "--n-cells", "64", "--y1plus", "0.5",
                 "--omega-plus", "0.02", "--amplitude", "5", "--periods", "4",
                 "--average", "2", "--steps-per-period", "32", "-o", str(out), "-q"]) == 0
    for f in ("summary.json", "history.csv", "phase_profiles.csv", "harmonic.csv",
              "final_profiles.csv", "urans_sa.png"):
        assert (out / "sa" / f).exists(), f


def test_cli_rejects_unknown_model(tmp_path, capsys):
    assert main(["rans", "-m", "foo", "-o", str(tmp_path), "--no-plot"]) == 2
    assert "Modèle inconnu" in capsys.readouterr().err


def test_cli_verify():
    assert main(["verify"]) == 0
