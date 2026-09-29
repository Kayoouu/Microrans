"""Cas de vérification du code (solutions exactes) — « est-ce que les équations sont bien résolues ? ».

La *validation* des modèles (« sont-ce les bonnes équations ? ») demande des données
DNS/expérimentales et n'est pas couverte ici au-delà de la corrélation de Dean.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid import tanh_grid
from .models import get_model
from .numerics import solve_transport
from .reference import laminar_poiseuille, womersley_channel
from .solver import Solution, solve_steady, solve_unsteady


@dataclass
class OrderStudy:
    name: str
    sizes: list
    errors: list
    expected_order: float

    @property
    def orders(self) -> list[float]:
        e = self.errors
        return [float(np.log(e[i] / e[i + 1]) / np.log(self.sizes[i + 1] / self.sizes[i]))
                for i in range(len(e) - 1)]

    @property
    def passed(self) -> bool:
        return abs(self.orders[-1] - self.expected_order) < 0.25

    def report(self) -> str:
        errs = "  ".join(f"{e:.2e}" for e in self.errors)
        ords = "  ".join(f"{o:.2f}" for o in self.orders)
        status = "OK " if self.passed else "ÉCHEC"
        return (f"[{status}] {self.name}\n        erreurs : {errs}\n"
                f"        ordres  : {ords}  (attendu {self.expected_order:g})")


def diffusion_mms(sizes=(16, 32, 64, 128), gamma_stretch: float = 2.0) -> OrderStudy:
    """Solution manufacturée : −d/dy((1+y²) du/dy) = s, u = sin(πy/2) sur [0, 2]."""
    errors = []
    for n in sizes:
        g = tanh_grid(n, gamma_stretch)
        y = g.y
        u_ex = np.sin(0.5 * np.pi * y)
        s = -(np.pi * y * np.cos(0.5 * np.pi * y)
              - (1 + y ** 2) * (np.pi ** 2 / 4) * np.sin(0.5 * np.pi * y))
        u = solve_transport(g, 1.0 + y ** 2, 0.0, s, (0.0, 0.0))
        errors.append(float(np.max(np.abs(u - u_ex))))
    return OrderStudy("Diffusion à coefficient variable (MMS), maillage étiré", list(sizes),
                      errors, 2.0)


def poiseuille_error(n: int = 64, nu: float = 0.01) -> float:
    """Erreur relative max du Poiseuille laminaire (le schéma est exact pour un polynôme de degré 2)."""
    g = tanh_grid(n, 2.0)
    m = get_model("laminar", g, nu)
    U0 = np.zeros(g.n)
    sol = solve_steady(m, g, nu, forcing=1.0, initial=Solution(m, g, nu, U0, {}),
                       dt=np.inf, tol=1e-14, max_iter=5)
    ex = laminar_poiseuille(g.y, 1.0, nu)
    return float(np.max(np.abs(sol.U - ex)) / np.max(ex))


def _womersley_error(n_cells, steps, scheme, nu=0.02, amp=5.0, omega=2 * np.pi):
    g = tanh_grid(n_cells, 2.0)
    m = get_model("laminar", g, nu)
    period = 2 * np.pi / omega
    U0 = womersley_channel(g.y, 0.0, 1.0, amp, omega, nu)
    res = solve_unsteady(m, g, nu, lambda t: 1.0 + amp * np.sin(omega * t),
                         Solution(m, g, nu, U0, {}), t_end=period, dt=period / steps,
                         scheme=scheme)
    ex = womersley_channel(g.y, period, 1.0, amp, omega, nu)
    return float(np.max(np.abs(res.final.U - ex)) / np.max(np.abs(ex)))


def womersley_time_order(scheme: str = "bdf2", steps=(20, 40, 80)) -> OrderStudy:
    errors = [_womersley_error(256, s, scheme) for s in steps]
    expected = 2.0 if scheme == "bdf2" else 1.0
    return OrderStudy(f"Womersley laminaire, convergence en temps ({scheme})", list(steps),
                      errors, expected)


def womersley_space_order(sizes=(16, 32, 64)) -> OrderStudy:
    errors = [_womersley_error(n, 2000, "bdf2") for n in sizes]
    return OrderStudy("Womersley laminaire, convergence en espace (BDF2, Δt petit)",
                      list(sizes), errors, 2.0)


def run_all(verbose: bool = True) -> bool:
    ok = True
    studies = [diffusion_mms(), womersley_time_order("euler"), womersley_time_order("bdf2"),
               womersley_space_order()]
    for s in studies:
        ok &= s.passed
        if verbose:
            print(s.report())
    err = poiseuille_error()
    passed = err < 1e-10
    ok &= passed
    if verbose:
        print(f"[{'OK ' if passed else 'ÉCHEC'}] Poiseuille laminaire stationnaire : "
              f"erreur relative {err:.1e} (schéma exact pour un profil parabolique)")
    return ok
