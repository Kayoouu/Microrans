"""Cas de calcul prêts à l'emploi : canal plan RANS et canal pulsé URANS."""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from .grid import Grid, channel_grid
from .models import canonical_name, get_model
from .reference import dean_bulk_velocity_plus, stokes_oscillation
from .solver import Solution, UnsteadyResult, solve_steady, solve_unsteady


@dataclass(frozen=True)
class PulsatingForcing:
    """f(t) = mean + amplitude · sin(ω t)  (f = −(1/ρ) ∂p/∂x, unités u_τ²/h)."""
    mean: float = 1.0
    amplitude: float = 0.0
    omega: float = 0.0

    def __call__(self, t):
        return self.mean + self.amplitude * np.sin(self.omega * t)


# ----------------------------------------------------------------------------- RANS
@dataclass
class RansResult:
    model_name: str
    re_tau: float
    solution: Solution
    wall_time: float

    @property
    def summary(self) -> dict:
        s = self.solution
        u_tau = s.u_tau
        tb, tt = s.wall_shear()
        ub_plus = s.bulk_velocity / u_tau
        re_tau_eff = u_tau * s.grid.h / s.nu
        re_b = s.bulk_velocity * 2.0 * s.grid.h / s.nu
        ub_dean = dean_bulk_velocity_plus(re_tau_eff)
        return {
            "model": self.model_name,
            "label": s.model.label,
            "re_tau_nominal": self.re_tau,
            "re_tau": re_tau_eff,
            "converged": bool(s.converged),
            "iterations": int(s.iterations),
            "final_residual": float(max(s.residuals[-1].values())) if s.residuals else None,
            "wall_time_s": round(self.wall_time, 3),
            "n_nodes": int(s.grid.n),
            "y1_plus": float(s.grid.first_cell_height * u_tau / s.nu),
            "tau_wall_bottom": tb,
            "tau_wall_top": tt,
            "u_tau": u_tau,
            "Ub_plus": ub_plus,
            "Uc_plus": s.centerline_velocity / u_tau,
            "Re_b": re_b,
            "Cf": 2.0 / ub_plus ** 2,
            "Ub_plus_dean": ub_dean,
            "Cf_dean": 2.0 / ub_dean ** 2,
            "Ub_plus_error_vs_dean_pct": 100.0 * (ub_plus / ub_dean - 1.0),
        }


def run_rans_channel(model: str = "sa", re_tau: float = 395.0, n_cells: int = 192,
                     y1_plus: float = 0.2, grid: Grid | None = None,
                     model_options: dict | None = None, **solver_kwargs) -> RansResult:
    """Canal plan turbulent établi, gradient de pression constant (u_τ nominale = 1)."""
    name = canonical_name(model)
    nu = 1.0 / re_tau
    grid = grid if grid is not None else channel_grid(n_cells, re_tau, y1_plus)
    m = get_model(name, grid, nu, **(model_options or {}))
    t0 = time.perf_counter()
    sol = solve_steady(m, grid, nu, forcing=1.0, **solver_kwargs)
    return RansResult(name, re_tau, sol, time.perf_counter() - t0)


# ---------------------------------------------------------------------------- URANS
@dataclass
class UransResult:
    model_name: str
    re_tau: float
    forcing: PulsatingForcing
    period: float
    steps_per_period: int
    history: UnsteadyResult
    steady: Solution
    phase: np.ndarray                       # phases ωt ∈ [0, 2π)
    phase_U: np.ndarray                     # (n_phases, n_nodes)
    phase_nut: np.ndarray
    harmonic_U: np.ndarray                  # 1er harmonique complexe de U(y)
    wall_time: float
    extra: dict = field(default_factory=dict)

    @property
    def grid(self) -> Grid:
        return self.steady.grid

    @property
    def nu(self) -> float:
        return self.steady.nu

    def laminar_stokes_harmonic(self) -> np.ndarray:
        """1er harmonique de la réponse laminaire au même forçage (référence Stokes)."""
        f = self.forcing
        return -1.0j * f.amplitude * stokes_oscillation(self.grid.y, f.omega, self.nu, self.grid.h)

    @property
    def summary(self) -> dict:
        h = self.history
        spp = self.steps_per_period
        tau = 0.5 * (h.tau_bottom + h.tau_top)
        last = slice(-spp, None)
        prev = slice(-2 * spp, -spp)
        omega = self.forcing.omega
        t_last = h.times[last]
        tau_h = 2.0 / spp * np.sum(tau[last] * np.exp(-1.0j * omega * t_last))
        ub_h = 2.0 / spp * np.sum(h.bulk_velocity[last] * np.exp(-1.0j * omega * t_last))
        # Forçage A sin(ωt) = Re[−iA e^{iωt}] : phase de référence −90°.
        ref = -1.0j
        periodicity = float(np.max(np.abs(tau[last] - tau[prev])) / max(np.max(np.abs(tau[last])), 1e-300))
        steady_tb, steady_tt = self.steady.wall_shear()
        return {
            "model": self.model_name,
            "label": self.steady.model.label,
            "re_tau_nominal": self.re_tau,
            "omega": omega,
            "omega_plus": omega / self.re_tau,
            "period": self.period,
            "period_plus": self.period * self.re_tau,
            "forcing_amplitude": self.forcing.amplitude,
            "stokes_thickness_plus": float(np.sqrt(2.0 * self.re_tau / omega)),
            "steps_per_period": spp,
            "n_periods": int(round((h.times[-1] - h.times[0]) / self.period)),
            "mean_inner_iterations": float(np.mean(h.inner_iterations[1:])),
            "max_inner_iterations": int(np.max(h.inner_iterations[1:])),
            "wall_time_s": round(self.wall_time, 3),
            "steady_tau_wall": 0.5 * (steady_tb + steady_tt),
            "steady_Ub": self.steady.bulk_velocity,
            "mean_tau_wall": float(np.mean(tau[last])),
            "mean_Ub": float(np.mean(h.bulk_velocity[last])),
            "tau_wall_amplitude": float(abs(tau_h)),
            "tau_wall_phase_deg": float(np.degrees(np.angle(tau_h / ref))),
            "Ub_amplitude": float(abs(ub_h)),
            "Ub_phase_deg": float(np.degrees(np.angle(ub_h / ref))),
            "periodicity_error": periodicity,
        }


def run_pulsating_channel(model: str = "sa", re_tau: float = 395.0, n_cells: int = 192,
                          y1_plus: float = 0.2, omega_plus: float = 0.01,
                          amplitude: float = 10.0, n_periods: int | None = None,
                          t_transient: float = 80.0, steps_per_period: int = 64,
                          n_average: int = 5, n_phases: int = 8, scheme: str = "sdirk2",
                          max_inner: int = 30, inner_tol: float = 1e-6, relax: float = 1.0,
                          model_options: dict | None = None, verbose: bool = False,
                          steady_kwargs: dict | None = None) -> UransResult:
    """Canal pulsé : f(t) = 1 + A sin(ωt), ω = ω⁺ Re_τ (unités h/u_τ).

    1) RANS stationnaire à f = 1 (condition initiale), 2) URANS sur n_periods,
    3) moyenne de phase et 1er harmonique sur les n_average dernières périodes.

    Par défaut n_periods = ceil(t_transient / T) + n_average : le débit moyen relaxe
    avec une constante de temps de l'ordre de 6 h/u_τ à Re_τ = 395, il faut donc
    ~80 h/u_τ avant d'atteindre le régime périodique (⟨τ_w⟩ → 1 à ~1e-4 près).
    """
    if steps_per_period % n_phases:                       # M19 : noms de la ligne de commande
        lo = steps_per_period // n_phases * n_phases
        raise ValueError(f"--steps-per-period = {steps_per_period} : doit être un multiple de "
                         f"{n_phases} (nombre de phases de la moyenne de phase), par exemple "
                         f"{max(lo, n_phases)} ou {lo + n_phases}.")
    omega = omega_plus * re_tau
    period = 2.0 * np.pi / omega
    if n_periods is None:
        n_periods = int(np.ceil(t_transient / period)) + n_average
    if n_average > n_periods:
        raise ValueError("n_average doit être <= n_periods.")
    name = canonical_name(model)
    nu = 1.0 / re_tau
    grid = channel_grid(n_cells, re_tau, y1_plus)
    m = get_model(name, grid, nu, **(model_options or {}))
    t0 = time.perf_counter()
    steady = solve_steady(m, grid, nu, forcing=1.0, **(steady_kwargs or {}))
    if not steady.converged:
        raise RuntimeError(f"Le RANS initial ({name}) n'a pas convergé.")

    dt = period / steps_per_period
    forcing = PulsatingForcing(1.0, amplitude, omega)
    t_end = n_periods * period
    store_from = (n_periods - n_average) * period + 0.5 * dt
    hist = solve_unsteady(m, grid, nu, forcing, steady, t_end=t_end, dt=dt, scheme=scheme,
                          max_inner=max_inner, inner_tol=inner_tol, relax=relax,
                          store_from=store_from, verbose=verbose,
                          log_every=max(steps_per_period, 1))

    n_store = n_average * steps_per_period
    U = hist.stored_U[-n_store:].reshape(n_average, steps_per_period, grid.n)
    nut = hist.stored_nut[-n_store:].reshape(n_average, steps_per_period, grid.n)
    t_store = hist.stored_times[-n_store:]
    phase_U_all = U.mean(axis=0)
    phase_nut_all = nut.mean(axis=0)
    stride = steps_per_period // n_phases
    # Échantillons du stockage : t = t_start + (j+1) dt → phase ωt modulo 2π
    phases = np.mod(omega * t_store[:steps_per_period], 2.0 * np.pi)
    harmonic = 2.0 / steps_per_period * np.sum(
        phase_U_all * np.exp(-1.0j * phases)[:, None], axis=0)
    idx = np.arange(stride - 1, steps_per_period, stride)
    idx = idx[np.argsort(phases[idx])]
    return UransResult(
        model_name=name, re_tau=re_tau, forcing=forcing, period=period,
        steps_per_period=steps_per_period, history=hist, steady=steady,
        phase=phases[idx], phase_U=phase_U_all[idx], phase_nut=phase_nut_all[idx],
        harmonic_U=harmonic, wall_time=time.perf_counter() - t0)
