"""Solutions de référence de la dynamique des gaz (gaz parfait, γ constant), pour valider
le solveur compressible (compressible.py).

- Problème de Riemann exact (tube à choc) : méthode de Toro, « Riemann Solvers and
  Numerical Methods for Fluid Dynamics », 3e éd., Springer 2009, chap. 4 (fonction de
  pression f_L + f_R + Δu = 0 résolue par Newton, puis échantillonnage de la solution
  autosemblable en ξ = x/t).
- Choc oblique attaché (relation θ-β-M, choc faible), choc droit : NACA Report 1135 (Ames
  Research Staff, « Equations, tables and charts for compressible flow », 1953), éq. 148
  à 150 ; Anderson, « Modern Compressible Flow », 3e éd., § 4.3.
- Relations isentropiques (conditions génératrices).
- Couche limite laminaire de Blasius : C_f = 0.664 / √Re_x (Schlichting & Gersten,
  « Boundary-Layer Theory », 8e éd., § 7.1 ; incompressible, bonne approximation pour
  M ≲ 0.3 à paroi adiabatique).
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------- tube à choc
def _pressure_function(p, rho, pk, ck, gamma):
    """f_K(p) et f_K'(p) de Toro (éq. 4.6, 4.7, 4.37)."""
    A = 2.0 / ((gamma + 1.0) * rho)
    B = (gamma - 1.0) / (gamma + 1.0) * pk
    if p > pk:                                             # choc
        q = np.sqrt(A / (p + B))
        return (p - pk) * q, q * (1.0 - 0.5 * (p - pk) / (B + p))
    # détente
    r = (p / pk) ** ((gamma - 1.0) / (2.0 * gamma))
    return (2.0 * ck / (gamma - 1.0) * (r - 1.0),
            1.0 / (rho * ck) * (p / pk) ** (-(gamma + 1.0) / (2.0 * gamma)))


def riemann_star(left, right, gamma: float = 1.4, tol: float = 1e-12):
    """État étoile (p*, u*) du problème de Riemann ; left/right = (ρ, u, p)."""
    rl, ul, pl = (float(v) for v in left)
    rr, ur, pr = (float(v) for v in right)
    cl, cr = np.sqrt(gamma * pl / rl), np.sqrt(gamma * pr / rr)
    if 2.0 * (cl + cr) / (gamma - 1.0) <= ur - ul:
        raise ValueError("Problème de Riemann : vide créé (condition de positivité violée).")
    # départ : approximation « two-rarefaction » (Toro éq. 4.46), robuste
    z = (gamma - 1.0) / (2.0 * gamma)
    p = ((cl + cr - 0.5 * (gamma - 1.0) * (ur - ul))
         / (cl / pl ** z + cr / pr ** z)) ** (1.0 / z)
    p = max(p, 1e-12 * (pl + pr))
    for _ in range(100):
        fl, dfl = _pressure_function(p, rl, pl, cl, gamma)
        fr, dfr = _pressure_function(p, rr, pr, cr, gamma)
        p_new = p - (fl + fr + ur - ul) / (dfl + dfr)
        p_new = max(p_new, 1e-12 * (pl + pr))
        if abs(p_new - p) < tol * 0.5 * (p_new + p):
            p = p_new
            break
        p = p_new
    fl, _ = _pressure_function(p, rl, pl, cl, gamma)
    fr, _ = _pressure_function(p, rr, pr, cr, gamma)
    return p, 0.5 * (ul + ur) + 0.5 * (fr - fl)


def riemann_exact(left, right, xi, gamma: float = 1.4):
    """Solution exacte (ρ, u, p) en ξ = (x − x0)/t (tableau), et dictionnaire des vitesses
    des ondes (tête/queue de détente, contact, choc) et de l'état étoile."""
    rl, ul, pl = (float(v) for v in left)
    rr, ur, pr = (float(v) for v in right)
    g = gamma
    cl, cr = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
    ps, us = riemann_star(left, right, g)
    xi = np.asarray(xi, float)
    rho, u, p = np.empty_like(xi), np.empty_like(xi), np.empty_like(xi)
    gm, gp = (g - 1.0) / (g + 1.0), 2.0 / (g + 1.0)
    waves = {"p_star": ps, "u_star": us, "contact": us}
    # --- côté gauche
    if ps > pl:                                            # choc gauche
        rls = rl * (ps / pl + gm) / (gm * ps / pl + 1.0)
        sl = ul - cl * np.sqrt((g + 1) / (2 * g) * ps / pl + (g - 1) / (2 * g))
        waves.update(left_shock=sl, rho_star_left=rls)
        L = xi < sl
        S = (xi >= sl) & (xi < us)
        rho[S], u[S], p[S] = rls, us, ps
    else:                                                  # détente gauche
        rls = rl * (ps / pl) ** (1.0 / g)
        cls = cl * (ps / pl) ** ((g - 1) / (2 * g))
        head, tail = ul - cl, us - cls
        waves.update(left_head=head, left_tail=tail, rho_star_left=rls)
        L = xi < head
        F = (xi >= head) & (xi < tail)
        c = gp * (cl + 0.5 * (g - 1) * (ul - xi[F]))
        u[F] = gp * (cl + 0.5 * (g - 1) * ul + xi[F])
        rho[F] = rl * (c / cl) ** (2.0 / (g - 1))
        p[F] = pl * (c / cl) ** (2.0 * g / (g - 1))
        S = (xi >= tail) & (xi < us)
        rho[S], u[S], p[S] = rls, us, ps
    rho[L], u[L], p[L] = rl, ul, pl
    # --- côté droit
    if ps > pr:                                            # choc droit
        rrs = rr * (ps / pr + gm) / (gm * ps / pr + 1.0)
        sr = ur + cr * np.sqrt((g + 1) / (2 * g) * ps / pr + (g - 1) / (2 * g))
        waves.update(right_shock=sr, rho_star_right=rrs)
        S = (xi >= us) & (xi < sr)
        rho[S], u[S], p[S] = rrs, us, ps
        R = xi >= sr
    else:                                                  # détente droite
        rrs = rr * (ps / pr) ** (1.0 / g)
        crs = cr * (ps / pr) ** ((g - 1) / (2 * g))
        head, tail = ur + cr, us + crs
        waves.update(right_head=head, right_tail=tail, rho_star_right=rrs)
        S = (xi >= us) & (xi < tail)
        rho[S], u[S], p[S] = rrs, us, ps
        F = (xi >= tail) & (xi < head)
        c = gp * (cr - 0.5 * (g - 1) * (ur - xi[F]))
        u[F] = gp * (-cr + 0.5 * (g - 1) * ur + xi[F])
        rho[F] = rr * (c / cr) ** (2.0 / (g - 1))
        p[F] = pr * (c / cr) ** (2.0 * g / (g - 1))
        R = xi >= head
    rho[R], u[R], p[R] = rr, ur, pr
    return rho, u, p, waves


# ----------------------------------------------------------------------- chocs
def normal_shock(M1: float, gamma: float = 1.4) -> dict:
    """Choc droit : rapports p2/p1, ρ2/ρ1, T2/T1, p02/p01 et Mach aval (NACA 1135)."""
    g = gamma
    m2 = M1 * M1
    if M1 <= 1.0:
        raise ValueError("Choc droit : M1 > 1 requis.")
    pr = 1.0 + 2.0 * g / (g + 1.0) * (m2 - 1.0)
    rr = (g + 1.0) * m2 / ((g - 1.0) * m2 + 2.0)
    M2 = np.sqrt((1.0 + 0.5 * (g - 1.0) * m2) / (g * m2 - 0.5 * (g - 1.0)))
    p0r = (rr ** (g / (g - 1.0))) * (pr ** (-1.0 / (g - 1.0)))
    return {"p_ratio": pr, "rho_ratio": rr, "T_ratio": pr / rr, "M2": float(M2),
            "p0_ratio": float(p0r)}


def oblique_shock(M1: float, theta_deg: float, gamma: float = 1.4) -> dict:
    """Choc oblique attaché, solution faible : angle β (degrés), p2/p1, ρ2/ρ1, T2/T1, M2."""
    from scipy.optimize import brentq
    g = gamma
    th = np.radians(theta_deg)

    def f(beta):
        return (2.0 / np.tan(beta) * (M1 ** 2 * np.sin(beta) ** 2 - 1.0)
                / (M1 ** 2 * (g + np.cos(2.0 * beta)) + 2.0) - np.tan(th))
    mu = np.arcsin(1.0 / M1)
    # β de déviation maximale : maximum de θ(β) sur ]μ, π/2[
    bb = np.linspace(mu + 1e-9, 0.5 * np.pi - 1e-9, 20001)
    tmax = bb[np.argmax(f(bb))]
    if f(tmax) < 0.0:
        raise ValueError(f"Choc détaché : θ = {theta_deg}° > déviation maximale à M = {M1}.")
    beta = brentq(f, mu + 1e-12, tmax, xtol=1e-14)
    Mn1 = M1 * np.sin(beta)
    ns = normal_shock(Mn1, g)
    M2 = ns["M2"] / np.sin(beta - th)
    return {"beta_deg": float(np.degrees(beta)), "p_ratio": ns["p_ratio"],
            "rho_ratio": ns["rho_ratio"], "T_ratio": ns["T_ratio"], "M2": float(M2),
            "theta_max_deg": float(np.degrees(np.arctan(f(tmax) + np.tan(th))))}


# ----------------------------------------------------------------------- isentropique
def total_to_static(M, gamma: float = 1.4) -> dict:
    """Rapports T0/T, p0/p, ρ0/ρ à Mach M (écoulement isentropique)."""
    t = 1.0 + 0.5 * (gamma - 1.0) * np.asarray(M, float) ** 2
    return {"T0_T": t, "p0_p": t ** (gamma / (gamma - 1.0)),
            "rho0_rho": t ** (1.0 / (gamma - 1.0))}


def blasius_cf(re_x):
    """Frottement pariétal laminaire de Blasius, C_f = 0.664 / √Re_x."""
    return 0.664 / np.sqrt(np.asarray(re_x, float))
