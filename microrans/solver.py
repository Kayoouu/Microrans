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
from .numerics import diffusion, solve_transport


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


class EvalStep:
    """« Pas » d'évaluation explicite : calcule F(φ) = S − D·φ + d/dy(Γ dφ/dy) pour chaque
    variable que le modèle résoudrait, sans rien résoudre (φ inchangé : évaluation de
    Jacobi, chaque équation voit l'état de départ). Utilisé par les schémas explicites."""

    def __init__(self, grid: Grid, state: dict[str, np.ndarray]):
        self.grid, self.state, self.F, self.gamma, self.sink = grid, state, {}, {}, {}

    def solve(self, name, gamma, source, sink, wall_values=(0.0, 0.0), model=None):
        phi = self.state[name]
        self.gamma[name] = np.broadcast_to(np.asarray(gamma, dtype=float), phi.shape)
        self.sink[name] = np.broadcast_to(np.asarray(sink, dtype=float), phi.shape)
        f = (np.asarray(source, dtype=float) - np.asarray(sink, dtype=float) * phi
             + diffusion(self.grid, np.asarray(gamma, dtype=float), phi))
        f[0] = f[-1] = 0.0
        self.F[name] = f
        return phi


# Schémas en temps du solveur 1D.
#   implicites : euler, bdf2 (multipas) ; cn (ESDIRK, 1er étage explicite, A-stable) ;
#                sdirk2, sdirk3 (Alexander 1977, L-stables, raides-précis)
#   explicites : rk1..rk4, ab2 (limités par Δt ≲ Δy²/ν_eff : voir `explicit_dt_limit`)
_G2 = 1.0 - 1.0 / np.sqrt(2.0)
_G3 = 0.435866521508459
_T3 = 0.5 * (1.0 + _G3)
_B31 = -(6 * _G3 ** 2 - 16 * _G3 + 1) / 4
_B32 = (6 * _G3 ** 2 - 20 * _G3 + 5) / 4
DIRK_TABLES = {   # (A triangulaire inférieure avec diagonale, c) ; b = dernière ligne
    "cn": ([[0.0, 0.0], [0.5, 0.5]], [0.0, 1.0]),
    "sdirk2": ([[_G2, 0.0], [1 - _G2, _G2]], [_G2, 1.0]),
    "sdirk3": ([[_G3, 0.0, 0.0], [_T3 - _G3, _G3, 0.0], [_B31, _B32, _G3]], [_G3, _T3, 1.0]),
}
ERK_TABLES = {
    "rk1": ([[0.0]], [1.0], [0.0]),
    "rk2": ([[0.0, 0.0], [1.0, 0.0]], [0.5, 0.5], [0.0, 1.0]),
    "rk3": ([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.25, 0.25, 0.0]], [1 / 6, 1 / 6, 2 / 3],
            [0.0, 1.0, 0.5]),
    "rk4": ([[0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0], [0.0, 0.5, 0.0, 0.0],
             [0.0, 0.0, 1.0, 0.0]], [1 / 6, 1 / 3, 1 / 3, 1 / 6], [0.0, 0.5, 0.5, 1.0]),
}
TIME_SCHEMES_1D = {
    "euler": "Euler implicite (ordre 1, L-stable)",
    "bdf2": "BDF2 (ordre 2, L-stable, multipas)",
    "cn": "Crank-Nicolson (ordre 2, A-stable non L-stable)",
    "sdirk2": "SDIRK 2 étages (ordre 2, L-stable)",
    "sdirk3": "SDIRK 3 étages (ordre 3, L-stable)",
    "rk1": "Euler explicite (ordre 1)",
    "rk2": "RK2 Heun (ordre 2, explicite)",
    "rk3": "RK3 SSP (ordre 3, explicite)",
    "rk4": "RK4 classique (ordre 4, explicite)",
    "ab2": "Adams-Bashforth 2 (ordre 2, explicite)",
}
# demi-longueur de l'intervalle de stabilité réel des schémas explicites
_REAL_STABILITY = {"rk1": 2.0, "rk2": 2.0, "rk3": 2.5127, "rk4": 2.7853, "ab2": 1.0}


def explicit_dt_limit(grid: Grid, gammas, sinks=None, scheme: str = "rk4") -> float:
    """Δt max de stabilité d'un schéma explicite (borne de Gershgorin du spectre).

    Pour chaque équation : λ_max ≤ max_i [2(Γ_w/Δy_w + Γ_e/Δy_e)/vol_i + D_i], avec Γ la
    diffusivité (U : ν + ν_t ; SA : (ν + ν̃)/σ...) et D le puits linéarisé (ex. 2βω pour
    ω, très raide près de la paroi en k-ω). En proche paroi (Δy ~ y⁺/Re_τ) cette limite
    est minuscule : c'est pourquoi le RANS/URANS résolu à la paroi utilise l'implicite.
    """
    if isinstance(gammas, np.ndarray):
        gammas = [gammas]
    sinks = list(sinks) if sinks is not None else [None] * len(gammas)
    lam = 0.0
    for gam, snk in zip(gammas, sinks):
        g = 0.5 * (gam[:-1] + gam[1:]) / grid.dy
        rate = 2.0 * (g[:-1] + g[1:]) / grid.volumes[1:-1]
        if snk is not None:
            rate = rate + np.maximum(snk[1:-1], 0.0)
        lam = max(lam, float(np.max(rate)))
    return float(_REAL_STABILITY[scheme] / lam)


def model_diffusivities(model, grid, nu, U, state):
    """(diffusivités, puits) de toutes les équations transportées (U et turbulence)."""
    flow = FlowField.from_velocity(grid, U)
    ev = EvalStep(grid, dict(state))
    model.update(dict(state), flow, ev)
    return ([nu + model.eddy_viscosity(state, flow), *ev.gamma.values()],
            [None, *ev.sink.values()])


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


def _explicit_rhs(model, grid, nu, f, Y):
    """F(Y, t) de tous les champs (U et turbulence), évaluation explicite."""
    U = Y["U"]
    state = {k: Y[k] for k in model.variables}
    flow = FlowField.from_velocity(grid, U)
    nut = model.eddy_viscosity(state, flow)
    FU = f + diffusion(grid, nu + nut, U)
    FU[0] = FU[-1] = 0.0
    ev = EvalStep(grid, state)
    model.update(state, flow, ev)
    return {"U": FU, **ev.F}


def _clip(model, Y):
    """Valeurs pariétales et planchers après un étage explicite."""
    Y["U"][0] = Y["U"][-1] = 0.0
    for k in model.variables:
        Y[k][0], Y[k][-1] = model.wall_values(k)
        if k in model.floors:
            Y[k][1:-1] = np.maximum(Y[k][1:-1], model.floors[k])
    return Y


def solve_unsteady(model: TurbulenceModel, grid: Grid, nu: float,
                   forcing: Callable[[float], float], initial: Solution,
                   t_end: float, dt: float, scheme: str = "bdf2",
                   max_inner: int = 30, inner_tol: float = 1e-8, relax: float = 1.0,
                   store_from: float | None = None, verbose: bool = False,
                   log_every: int = 1000) -> UnsteadyResult:
    """URANS : intégration en temps physique, schéma au choix (voir TIME_SCHEMES_1D).

    Implicites : sous-itérations de Picard à chaque pas (ou étage) jusqu'à inner_tol.
    Explicites : aucune sous-itération, mais Δt limité par `explicit_dt_limit`.
    Les profils U, ν_t (et les variables de turbulence) sont stockés à chaque pas pour
    t > store_from (moyennes de phase, analyse harmonique).
    """
    if scheme not in TIME_SCHEMES_1D:
        raise ValueError(f"scheme doit être l'un de {list(TIME_SCHEMES_1D)}.")
    t0 = initial.time
    n_steps = int(round((t_end - t0) / dt))
    if n_steps < 1:
        raise ValueError("t_end doit être > temps initial + dt.")
    linear = len(model.variables) == 0

    current = {"U": initial.U.copy(), **{k: v.copy() for k, v in initial.state.items()}}
    previous = None
    F_prev = None
    rec = {k: np.empty(n_steps + 1) for k in ("t", "f", "tb", "tt", "ub", "uc", "inner")}

    def record(i, t, f, U, n_inner):
        tb, tt = wall_shear(grid, nu, U)
        rec["t"][i], rec["f"][i], rec["tb"][i], rec["tt"][i] = t, f, tb, tt
        rec["ub"][i] = bulk_velocity(grid, U)
        rec["uc"][i] = centerline_velocity(grid, U)
        rec["inner"][i] = n_inner

    def implicit_solve(step, f, U, state):
        source = np.full(grid.n, f)
        inner = 0
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
        return U, state, inner

    if scheme in ERK_TABLES or scheme == "ab2":
        dt_lim = explicit_dt_limit(grid, *model_diffusivities(
            model, grid, nu, initial.U, initial.state), scheme=scheme)
        if dt > dt_lim:
            import warnings
            warnings.warn(f"{scheme} : Δt = {dt:.3g} > limite de stabilité explicite "
                          f"≈ {dt_lim:.3g} (diffusion près des parois) : divergence probable.")

    record(0, t0, forcing(t0), current["U"], 0)
    stored_t, stored_u, stored_nut = [], [], []
    stored_extra = {k: [] for k in model.variables}

    for n in range(1, n_steps + 1):
        t = t0 + n * dt
        tn = t - dt
        f = float(forcing(t))
        if scheme in ("euler", "bdf2"):
            if scheme == "euler" or previous is None:
                step = ImplicitStep.euler(grid, dt, current)
            else:
                step = ImplicitStep.bdf2(grid, dt, current, previous)
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
            U, state, inner = implicit_solve(step, f, U, state)
            new = {"U": U, **state}
        elif scheme in DIRK_TABLES:
            A, c = DIRK_TABLES[scheme]
            Fs, inner = [], 0
            Y = current
            for i in range(len(c)):
                ti = tn + c[i] * dt
                if A[i][i] == 0.0:                   # étage explicite (ESDIRK)
                    Fs.append(_explicit_rhs(model, grid, nu, float(forcing(ti)), current))
                    continue
                a0 = 1.0 / (A[i][i] * dt)
                known = {k: current[k] + dt * sum(A[i][j] * Fs[j][k] for j in range(i))
                         for k in current}
                step = ImplicitStep(grid, a0, {k: -a0 * v for k, v in known.items()})
                U, state, it_i = implicit_solve(step, float(forcing(ti)), Y["U"].copy(),
                                                {k: Y[k].copy() for k in model.variables})
                inner += it_i
                Y = {"U": U, **state}
                Fs.append({k: a0 * (Y[k] - known[k]) for k in Y})
            new = Y                                  # schémas raides-précis : yⁿ⁺¹ = dernier étage
        elif scheme in ERK_TABLES:
            A, b, c = ERK_TABLES[scheme]
            Fs = []
            for i in range(len(b)):
                Yi = current if i == 0 else _clip(model, {
                    k: current[k] + dt * sum(A[i][j] * Fs[j][k] for j in range(i) if A[i][j])
                    for k in current})
                Fs.append(_explicit_rhs(model, grid, nu, float(forcing(tn + c[i] * dt)), Yi))
            new = _clip(model, {k: current[k] + dt * sum(bj * Fj[k] for bj, Fj in zip(b, Fs))
                                for k in current})
            inner = len(b)
        else:                                        # ab2
            Fn = _explicit_rhs(model, grid, nu, float(forcing(tn)), current)
            if F_prev is None:
                new = {k: current[k] + dt * Fn[k] for k in current}
            else:
                new = {k: current[k] + dt * (1.5 * Fn[k] - 0.5 * F_prev[k]) for k in current}
            new = _clip(model, new)
            F_prev = Fn
            inner = 1
        U = new["U"]
        state = {k: new[k] for k in model.variables}
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
