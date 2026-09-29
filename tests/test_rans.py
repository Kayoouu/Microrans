"""Canal plan turbulent stationnaire : convergence, bilans et propriétés connues des modèles."""
import numpy as np
import pytest

from microrans.cases import run_rans_channel
from microrans.models import TURBULENT_MODELS
from microrans.numerics import ddy


@pytest.fixture(scope="module")
def solutions():
    return {m: run_rans_channel(m, re_tau=395.0, n_cells=128, y1_plus=0.5)
            for m in TURBULENT_MODELS}


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_converges_to_machine_precision(solutions, model):
    s = solutions[model].summary
    assert s["converged"], s
    assert s["final_residual"] < 1e-10


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_wall_shear_matches_momentum_balance(solutions, model):
    # Écoulement établi : τ_w = f·h = 1 exactement (à l'erreur de la dérivée pariétale près).
    s = solutions[model].summary
    assert s["tau_wall_bottom"] == pytest.approx(1.0, abs=1e-4)
    assert s["tau_wall_top"] == pytest.approx(1.0, abs=1e-4)


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_symmetry_and_positivity(solutions, model):
    sol = solutions[model].solution
    assert np.allclose(sol.U, sol.U[::-1], rtol=1e-8, atol=1e-10)
    for name, v in sol.state.items():
        assert np.all(v[1:-1] > 0.0), name
    assert np.all(sol.nut >= 0.0)
    assert sol.nut[0] == 0.0 and sol.nut[-1] == 0.0


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_bulk_velocity_in_turbulent_range(solutions, model):
    # Garde-fou grossier (PAS une validation) : un écoulement relaminarisé donnerait
    # U_b+ = Re_τ/3 ≈ 132 ; la corrélation de Dean donne 17.2. Les modèles sont à ±10 %.
    s = solutions[model].summary
    assert abs(s["Ub_plus_error_vs_dean_pct"]) < 12.0, s


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_viscous_sublayer(solutions, model):
    sol = solutions[model].solution
    yp = sol.grid.wall_distance / sol.nu
    mask = (yp > 0) & (yp < 2.0)
    assert np.allclose(sol.U[mask], yp[mask], rtol=0.03)


@pytest.fixture(scope="module")
def sa_high_re():
    return run_rans_channel("sa", re_tau=2000.0, n_cells=192, y1_plus=0.5).solution


def test_sa_nu_tilde_is_kappa_y_near_wall(sa_high_re):
    # Propriété de conception de SA : ν̃ = κ u_τ y dans la couche interne (y/h ≪ 1).
    sol = sa_high_re
    yp = sol.grid.wall_distance / sol.nu
    mask = (yp > 1.0) & (yp < 30.0)
    ratio = sol.state["nu_tilde"][mask] / (0.41 * yp[mask] * sol.nu)
    assert np.all(np.abs(ratio - 1.0) < 0.02)


def test_sa_log_law_slope(sa_high_re):
    sol = sa_high_re
    yp = sol.grid.wall_distance / sol.nu
    xi = yp * ddy(sol.grid, sol.U) * sol.nu  # y+ dU+/dy+ = 1/κ dans la zone log
    mask = (yp > 50) & (yp < 100) & (sol.grid.y < 1.0)
    assert np.all(np.abs(xi[mask] * 0.41 - 1.0) < 0.03)


@pytest.mark.parametrize("model", ["ke", "kw", "sst"])
def test_two_equation_k_follows_local_stress(model):
    # Équilibre production = dissipation : k = τ(y)/√C_μ, avec τ(y) = u_τ² (1 − y/h).
    r = run_rans_channel(model, re_tau=2000.0, n_cells=192, y1_plus=0.5)
    sol = r.solution
    yp = sol.grid.wall_distance / sol.nu
    i = np.argmin(np.abs(yp[: sol.grid.n // 2] - 200.0))
    expected = (1.0 - sol.grid.y[i]) / np.sqrt(0.09)
    assert sol.state["k"][i] == pytest.approx(expected, rel=0.05)


@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_grid_convergence(model):
    coarse = run_rans_channel(model, re_tau=395.0, n_cells=192, y1_plus=0.2).summary
    fine = run_rans_channel(model, re_tau=395.0, n_cells=384, y1_plus=0.1).summary
    assert coarse["Ub_plus"] == pytest.approx(fine["Ub_plus"], rel=0.01)


def test_sa_ft2_variant_runs():
    r = run_rans_channel("sa", re_tau=395.0, n_cells=128, y1_plus=0.5,
                         model_options={"ft2": True})
    assert r.summary["converged"]
    assert "f_t2" in r.solution.model.label


@pytest.mark.parametrize("re_tau", [180.0, 5200.0])
@pytest.mark.parametrize("model", TURBULENT_MODELS)
def test_robust_over_reynolds_range(model, re_tau):
    r = run_rans_channel(model, re_tau=re_tau, n_cells=192, y1_plus=0.5)
    assert r.summary["converged"], r.summary
