"""Intégration en temps : RANS stationnaire (pseudo-temps) et URANS (temps physique).

Équation de quantité de mouvement du canal (écoulement établi, U = U(y, t)) :
    ∂U/∂t = f(t) + ∂/∂y[(ν + ν_t) ∂U/∂y],   f = −(1/ρ) ∂p/∂x
Unités : longueurs en h, vitesses en u_τ nominale (f moyen = 1), temps en h/u_τ, ν = 1/Re_τ.

Résolution ségréguée : U puis les variables de turbulence, chacune par un système
tridiagonal implicite. En URANS, des sous-itérations (Picard) convergent le couplage
non linéaire à chaque pas de temps.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .flow import FlowField, bulk_velocity, centerline_velocity, initial_guess, wall_shear
from .grid import Grid
from .models.base import TurbulenceModel
from .numerics import solve_transport


class ImplicitStep:
    """Un pas implicite  a0·φ^{n+1} + hist = S_exp − D_imp·φ^{n+1} + ∂/∂y(Γ ∂φ/∂y).

    Euler implicite : a0 = 1/Δt,  hist = −φ^n/Δt
    BDF2            : a0 = 3/(2Δt), hist = (−4φ^n + φ^{n−1})/(2Δt)
    """

    def __init__(self, grid: Grid, a0: float, history: dict[str, np.ndarray]):
        self.grid = grid
        self.a0 = a0
        self.history = history

    def solve(self, name, gamma, source, sink, wall_values=(0.0, 0.0), model=None):
        if model is not None:
            wall_values = model.wall_values(name)
        rhs = np.asarray(source, dtype=float) - self.history[name]
        diag = self.a0 + np.asarray(sink, dtype=float)
        phi = solve_transport(self.grid, gamma, diag, rhs, wall_values)
        if model is not None and name in model.floors:
            phi[1:-1] = np.maximum(phi[1:-1], model.floors[name])
        return phi

    @classmethod
    def euler(cls, grid, dt, current):
        if np.isinf(dt):
            return cls(grid, 0.0, {k: np.zeros_like(v) for k, v in current.items()})
        return cls(grid, 1.0 / dt, {k: -v / dt for k, v in current.items()})

    @classmethod
    def bdf2(cls, grid, dt, current, previous):
        return cls(grid, 1.5 / dt,
                   {k: (-2.0 * current[k] + 0.5 * previous[k]) / dt for k in current})


@dataclass
class Solution:
    """État d'une solution (stationnaire ou instantanée)."""
    model: TurbulenceModel
    grid: Grid
    nu: float
    U: np.ndarray
    state: dict[str, np.ndarray]
    time: float = 0.0
    iterations: int = 0
    converged: bool = True
    residuals: list[dict[str, float]] = field(default_factory=list)

    @property
    def flow(self) -> FlowField:
        return FlowField.from_velocity(self.grid, self.U)

    @property
    def nut(self) -> np.ndarray:
        return self.model.eddy_viscosity(self.state, self.flow)

    def wall_shear(self):
        return wall_shear(self.grid, self.nu, self.U)

    @property
    def u_tau(self) -> float:
        tb, tt = self.wall_shear()
        return float(np.sqrt(max(0.5 * (tb + tt), 0.0)))

    @property
    def bulk_velocity(self) -> float:
        return bulk_velocity(self.grid, self.U)

    @property
    def centerline_velocity(self) -> float:
        return centerline_velocity(self.grid, self.U)


def _relative_change(new: dict, old: dict) -> dict[str, float]:
    out = {}
    for k in new:
        scale = max(np.max(np.abs(new[k])), 1e-300)
        out[k] = float(np.max(np.abs(new[k] - old[k])) / scale)
    return out


def _check_finite(fields: dict, where: str):
    for k, v in fields.items():
        if not np.all(np.isfinite(v)):
            raise FloatingPointError(f"Valeurs non finies dans '{k}' ({where}).")


def initial_solution(model: TurbulenceModel, grid: Grid, nu: float,
                     u_tau: float = 1.0) -> Solution:
    """Champ initial turbulent (profil de Cess) adapté au modèle."""
    U0, nut0 = initial_guess(grid, nu, u_tau)
    state = model.initial_state(FlowField.from_velocity(grid, U0), nut0)
    return Solution(model, grid, nu, U0, state, converged=False)


def solve_steady(model: TurbulenceModel, grid: Grid, nu: float, forcing: float = 1.0,
                 initial: Solution | None = None, dt: float = 5.0,
                 max_iter: int = 20000, tol: float = 1e-10, relax: float = 0.5,
                 verbose: bool = False, log_every: int = 200) -> Solution:
    """RANS stationnaire par marche en pseudo-temps (Euler implicite, Δτ global).

    - `dt` : pas de pseudo-temps (unités h/u_τ) ; `np.inf` = itérations de Picard pures.
    - `relax` : sous-relaxation des variables de turbulence. Le couplage ségrégué
      U ↔ ν_t se comporte comme x ↦ a/x (pente −1, oscillation de période 2) ;
      relax = 0.5 ramène la pente à ~0. Testé de Re_τ = 180 à 5200 pour les 4 modèles.
    - critère d'arrêt : variation relative max par itération de chaque variable < tol.
    """
    sol = initial if initial is not None else initial_solution(
        model, grid, nu, np.sqrt(max(forcing, 1e-300) * grid.h))
    U, state = sol.U.copy(), {k: v.copy() for k, v in sol.state.items()}
    source = np.full(grid.n, float(forcing))
    residuals = []
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        old = {"U": U, **state}
        flow = FlowField.from_velocity(grid, U)
        nut = model.eddy_viscosity(state, flow)
        step = ImplicitStep.euler(grid, dt, old)
        U = step.solve("U", nu + nut, source, 0.0)
        flow = FlowField.from_velocity(grid, U)
        state_new = model.update(state, flow, step)
        state = {k: state[k] + relax * (state_new[k] - state[k]) for k in state_new}
        new = {"U": U, **state}
        _check_finite(new, f"itération {it}")
        res = _relative_change(new, old)
        residuals.append(res)
        if verbose and (it % log_every == 0 or it == 1):
            txt = "  ".join(f"{k}={v:.2e}" for k, v in res.items())
            print(f"  [{model.name}] it {it:6d}  {txt}")
        if max(res.values()) < tol:
            converged = True
            break
    if verbose:
        status = "convergé" if converged else "NON convergé"
        print(f"  [{model.name}] {status} en {it} itérations "
              f"(résidu max {max(residuals[-1].values()):.2e})")
    return Solution(model, grid, nu, U, state, iterations=it, converged=converged,
                    residuals=residuals)


@dataclass
class UnsteadyResult:
    """Historique d'un calcul URANS."""
    final: Solution
    times: np.ndarray
    forcing: np.ndarray
    tau_bottom: np.ndarray
    tau_top: np.ndarray
    bulk_velocity: np.ndarray
    centerline_velocity: np.ndarray
    inner_iterations: np.ndarray
    stored_times: np.ndarray
    stored_U: np.ndarray
    stored_nut: np.ndarray
    stored_extra: dict[str, np.ndarray] = field(default_factory=dict)


def solve_unsteady(model: TurbulenceModel, grid: Grid, nu: float,
                   forcing: Callable[[float], float], initial: Solution,
                   t_end: float, dt: float, scheme: str = "bdf2",
                   max_inner: int = 30, inner_tol: float = 1e-8, relax: float = 1.0,
                   store_from: float | None = None, verbose: bool = False,
                   log_every: int = 1000) -> UnsteadyResult:
    """URANS : intégration en temps physique (Euler implicite ou BDF2) avec sous-itérations.

    Les profils U, ν_t (et les variables de turbulence) sont stockés à chaque pas pour
    t > store_from (moyennes de phase, analyse harmonique).
    """
    if scheme not in ("bdf2", "euler"):
        raise ValueError("scheme doit être 'bdf2' ou 'euler'.")
    t0 = initial.time
    n_steps = int(round((t_end - t0) / dt))
    if n_steps < 1:
        raise ValueError("t_end doit être > temps initial + dt.")
    linear = len(model.variables) == 0

    current = {"U": initial.U.copy(), **{k: v.copy() for k, v in initial.state.items()}}
    previous = None
    rec = {k: np.empty(n_steps + 1) for k in ("t", "f", "tb", "tt", "ub", "uc", "inner")}

    def record(i, t, f, U, n_inner):
        tb, tt = wall_shear(grid, nu, U)
        rec["t"][i], rec["f"][i], rec["tb"][i], rec["tt"][i] = t, f, tb, tt
        rec["ub"][i] = bulk_velocity(grid, U)
        rec["uc"][i] = centerline_velocity(grid, U)
        rec["inner"][i] = n_inner

    record(0, t0, forcing(t0), current["U"], 0)
    stored_t, stored_u, stored_nut = [], [], []
    stored_extra = {k: [] for k in model.variables}

    for n in range(1, n_steps + 1):
        t = t0 + n * dt
        f = float(forcing(t))
        if scheme == "euler" or previous is None:
            step = ImplicitStep.euler(grid, dt, current)
        else:
            step = ImplicitStep.bdf2(grid, dt, current, previous)
        source = np.full(grid.n, f)
        if previous is None:
            U = current["U"]
            state = {k: current[k] for k in model.variables}
        else:
            # Prédicteur : extrapolation linéaire en temps (réduit les sous-itérations).
            U = 2.0 * current["U"] - previous["U"]
            state = {}
            for k in model.variables:
                guess = 2.0 * current[k] - previous[k]
                state[k] = np.where(guess > 0.0, guess, current[k])
        for inner in range(1, max_inner + 1):
            nut = model.eddy_viscosity(state, FlowField.from_velocity(grid, U))
            U_new = step.solve("U", nu + nut, source, 0.0)
            flow = FlowField.from_velocity(grid, U_new)
            state_new = model.update(state, flow, step)
            state_new = {k: state[k] + relax * (state_new[k] - state[k]) for k in state_new}
            change = _relative_change({"U": U_new, **state_new}, {"U": U, **state})
            U, state = U_new, state_new
            if linear or max(change.values()) < inner_tol:
                break
        new = {"U": U, **state}
        _check_finite(new, f"t = {t:.4g}")
        previous, current = current, new
        record(n, t, f, U, inner)
        if store_from is not None and t > store_from - 1e-12 * max(1.0, abs(t)):
            stored_t.append(t)
            stored_u.append(U.copy())
            stored_nut.append(model.eddy_viscosity(state, FlowField.from_velocity(grid, U)))
            for k in model.variables:
                stored_extra[k].append(state[k].copy())
        if verbose and (n % log_every == 0 or n == n_steps):
            print(f"  [{model.name}] pas {n}/{n_steps}  t={t:.4f}  "
                  f"τ_w={0.5 * (rec['tb'][n] + rec['tt'][n]):.4f}  "
                  f"U_b={rec['ub'][n]:.4f}  sous-it={inner}")

    final = Solution(model, grid, nu, current["U"],
                     {k: current[k] for k in model.variables}, time=t0 + n_steps * dt)
    empty = np.empty((0, grid.n))
    return UnsteadyResult(
        final=final, times=rec["t"], forcing=rec["f"], tau_bottom=rec["tb"],
        tau_top=rec["tt"], bulk_velocity=rec["ub"], centerline_velocity=rec["uc"],
        inner_iterations=rec["inner"].astype(int),
        stored_times=np.array(stored_t),
        stored_U=np.array(stored_u) if stored_u else empty,
        stored_nut=np.array(stored_nut) if stored_nut else empty,
        stored_extra={k: np.array(v) for k, v in stored_extra.items() if v})
