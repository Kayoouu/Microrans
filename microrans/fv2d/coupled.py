"""Solveur stationnaire COUPLÉ pression-vitesse ([solver] algorithm = "coupled"), à la
manière du « Coupled » de Fluent et de pUCoupledFoam (foam-extend).

À chaque itération, vitesse (u, v) et pression p sont résolues ENSEMBLE dans un seul
système linéaire creux de 3 N inconnues, au lieu d'alterner prédiction de la vitesse et
correction de pression (SIMPLE). Blocs :

    [ A_u   0    G_x ] [u]   [b_u]      A : convection-diffusion (équation de quantité de
    [ 0    A_v   G_y ] [v] = [b_v]          mouvement du solveur SIMPLE, même assemblage)
    [ D_x  D_y  -L   ] [p]   [b_p]      G : gradient de Green-Gauss (p aux faces, CL de p)

La ligne de continuité est le bilan des flux de Rhie-Chow du solveur segmenté :
F_f = (Ū_f + (r_AU ∇p)‾_f)·S_f − (r_AU)_f g_f (p_N − p_P) (+ corrections), avec Ū_f et
p IMPLICITES (D, L) et le terme (r_AU ∇p)‾_f, les corrections non orthogonales et les
forces aux faces calculés avec l'itéré précédent : à convergence, les flux sont
exactement ceux de SIMPLE, donc la solution est la même. Les équations de turbulence,
de température et des scalaires restent résolues séparément après chaque itération
couplée (comme dans Fluent).

Résolution : LU creuse (SuperLU) jusqu'à direct_max cellules, sinon GMRES préconditionné
par une factorisation incomplète (ILU). Sous-relaxation implicite de la quantité de
mouvement (relax_U, défaut 0.9 en couplé) ; pas de relaxation de la pression.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

DIRECT_MAX = 60000


def _coo(rows, cols, vals, shape):
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=shape)


def coupled_step(s, eqs, relax=1.0):
    """Une itération couplée ; `eqs` = équations de quantité de mouvement (sans gradient
    de pression) déjà relaxées par `relax`. Met à jour s.U, s.p, s.F_i, s.F_b.

    Les flux de Rhie-Chow utilisent la diagonale NON relaxée (a_P · relax) : la solution
    convergée ne dépend pas du facteur de sous-relaxation (correction de Majumdar 1988),
    contrairement au SIMPLE(C) segmenté qui utilise la diagonale relaxée."""
    if s.backend.is_gpu:
        raise ValueError("algorithm = coupled : CPU seulement.")
    fvm = s.fvm
    nc, ni, nb = fvm.nc, fvm.ni, fvm.nb
    P, N, Pb = np.asarray(fvm.P), np.asarray(fvm.N), np.asarray(fvm.Pb)
    w, g = fvm.w, fvm.g
    Si, Sb, magSb, V = fvm.Si, fvm.Sb, fvm.magSb, fvm.V
    U, p = s.U, s.p

    diag_u, diag_v = eqs[0][0], eqs[1][0]
    aP = 0.5 * (diag_u + diag_v) * relax
    rAU = V / aP
    rA_f = fvm.interp(rAU)

    rows, cols, vals = [], [], []
    rhs = np.zeros(3 * nc)
    cells = np.arange(nc)

    # --- blocs quantité de mouvement A_c (u : lignes 0..nc-1, v : nc..2nc-1)
    for c, (diag, up, lo, b) in enumerate(eqs):
        o = c * nc
        rows += [o + cells, o + P, o + N]
        cols += [o + cells, o + N, o + P]
        vals += [diag, up, lo]
        rhs[o:o + nc] = b

    # --- blocs G (gradient de pression × V) : Σ p_f S_f, p_f = w p_P + (1 − w) p_N
    pa, pbeta, pg, pd = s.pressure_bc()
    for c in range(2):
        o, q = c * nc, 2 * nc
        Sc = Si[:, c]
        rows += [o + P, o + P, o + N, o + N]
        cols += [q + P, q + N, q + P, q + N]
        vals += [w * Sc, (1 - w) * Sc, -w * Sc, -(1 - w) * Sc]
        if nb:
            rows.append(o + Pb)
            cols.append(q + Pb)
            vals.append(pa * Sb[:, c])
            rhs[o:o + nc] -= np.bincount(Pb, pbeta * Sb[:, c], nc)
        if fvm.axisymmetric:                           # faces latérales du secteur
            rows.append(o + cells)
            cols.append(q + cells)
            vals.append(-V * fvm._side[:, c])

    # --- continuité : Σ flux sortants = 0
    q = 2 * nc
    gradp = fvm.grad(p, s.boundary_p(p))
    # partie explicite des flux internes : (r_AU ∇p)‾·S (+ forces aux faces), non orthogonal
    E_i = np.sum(fvm.interp(rAU[:, None] * gradp) * Si, axis=1)
    if s._has_force():
        fc, ff = s.cell_force(), s.cell_force(face=True)
        E_i = E_i + np.sum((rA_f[:, None] * ff - fvm.interp(rAU[:, None] * fc)) * Si, axis=1)
    nf = 0.0
    if not fvm.orthogonal:
        nf = fvm.nonorth_flux(rA_f, gradp, p, s.settings.nonorth_limit)
    D = rA_f * g
    for c in range(2):
        Sc = Si[:, c]
        rows += [q + P, q + P, q + N, q + N]
        cols += [c * nc + P, c * nc + N, c * nc + P, c * nc + N]
        vals += [w * Sc, (1 - w) * Sc, -w * Sc, -(1 - w) * Sc]
    rows += [q + P, q + P, q + N, q + N]
    cols += [q + P, q + N, q + P, q + N]
    vals += [D, -D, -D, D]
    rhs[q:] -= np.bincount(P, E_i - nf, nc) - np.bincount(N, E_i - nf, nc)
    # frontières : U_b·S (vitesse imposée : constante ; gradient nul : U_P implicite +
    # r_AU ∇p_P·S), pression imposée : − r_AU |S| ∂p/∂n (implicite en p_P) ; glissement :
    # flux nul
    Fb_const = np.zeros(nb)
    if nb:
        kU = np.asarray(s.kindU)
        fixed, zg = kU == 0, kU == 1
        Fb_const += np.where(fixed, np.sum(s.U_fixed * Sb, axis=1), 0.0)
        rAb = rAU[Pb]
        for c in range(2):
            rows.append(q + Pb)
            cols.append(c * nc + Pb)
            vals.append(np.where(zg, Sb[:, c], 0.0))
        Fb_const += np.where(zg, rAb * np.sum(gradp[Pb] * Sb, axis=1), 0.0)
        if s._has_force():
            fb = np.sum(s._boundary_force() * Sb, axis=1)
            Fb_const += np.where(np.asarray(s.kindP) == 1, rAb * fb, 0.0)
        rows.append(q + Pb)
        cols.append(q + Pb)
        vals.append(-rAb * magSb * pg)
        Fb_const -= rAb * magSb * pd
        rhs[q:] -= np.bincount(Pb, Fb_const, nc)
    if not s.has_fixed_p:
        # niveau de pression indéterminé : diagonale de p renforcée dans la 1re cellule
        # (comme le SIMPLE segmenté) ; sans effet à convergence (p₀ ne bouge plus)
        pin = float(np.sum(D[(P == 0) | (N == 0)])) or 1.0
        rows.append(np.array([q]))
        cols.append(np.array([q]))
        vals.append(np.array([pin]))
        rhs[q] += pin * p[0]
    A = _coo(rows, cols, vals, (3 * nc, 3 * nc)).tocsr()
    x0 = np.concatenate([U[:, 0], U[:, 1], p])
    if nc <= DIRECT_MAX:
        x = spla.spsolve(A.tocsc(), rhs, permc_spec="COLAMD")
    else:
        ilu = spla.spilu(A.tocsc(), drop_tol=1e-4, fill_factor=10)
        M = spla.LinearOperator(A.shape, ilu.solve)
        x, info = spla.gmres(A, rhs, x0=x0, M=M, rtol=1e-6, restart=60, maxiter=20)
    r = rhs - A @ x0
    s.residuals_now["continuity_coupled"] = float(np.sum(np.abs(r[q:]))) / max(
        float(np.sum(np.abs(A[q:] @ x0))) + 1e-300, 1e-300)
    s.U = np.column_stack([x[:nc], x[nc:2 * nc]])
    s.p = x[q:]
    # flux conservatifs (même expression que la ligne de continuité)
    Uf = fvm.interp(s.U)
    s.F_i = np.sum(Uf * Si, axis=1) + E_i - D * (s.p[N] - s.p[P]) - nf
    if nb:
        s.F_b = (Fb_const + np.where(zg, np.sum(s.U[Pb] * Sb, axis=1), 0.0)
                 - rAb * magSb * pg * s.p[Pb])
