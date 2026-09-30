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
forces aux faces calculés avec l'itéré précédent (ils se compensent à convergence).
r_AU = V/a_P avec la diagonale NON relaxée : la solution convergée ne dépend pas de
relax_U (Majumdar 1988) ; elle ne diffère de celle de SIMPLEC (diagonale relaxée, r_AtU)
que par la dissipation de Rhie-Chow, c.-à-d. à l'ordre de l'erreur de discrétisation
(0,5 % de la vitesse du couvercle au plus, cavité 24 × 24, dans les coins). Turbulence,
température et scalaires restent résolus séparément après chaque itération couplée
(comme dans Fluent) : c'est eux qui limitent le gain en turbulent et en convection
naturelle.

Résolution linéaire : inconnues entrelacées par cellule, LU creuse (SuperLU, ordre
minimum-degré) jusqu'à DIRECT_MAX cellules, ILU au-delà (non mesuré). La factorisation est
RÉUTILISÉE comme préconditionneur de GMRES aux itérations suivantes tant qu'elle reste
efficace : une factorisation pour 7 à 60 itérations couplées sur les cas d'exemple. Sous-relaxation
implicite de la quantité de mouvement (relax_U, 1 par défaut en couplé) ; pas de
relaxation de la pression. CPU seulement.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import solve_triangular

DIRECT_MAX = 200000     # cellules : LU complète en dessous (160 000 : 1.7 Go), ILU au-dessus
INNER_TOL = 1e-3        # réduction du résidu linéaire demandée à chaque itération couplée
KRYLOV_MAX = 20         # itérations GMRES avant de refactoriser immédiatement
REFACTOR_AFTER = 6      # au-delà, factorisation jugée périmée : refaite à l'itération suivante


def _coo(rows, cols, vals, shape):
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=shape)


def _factor(A, exact, pivot=False):
    """LU (ou ILU) de la matrice couplée entrelacée.

    Par défaut : ordre minimum-degré sur A + Aᵀ et pivots pris SUR LA DIAGONALE, sans
    pivotage : la structure 3 × 3 par cellule est respectée. Le pivotage partiel
    (seuillé) de SuperLU est à proscrire ici : le coefficient de p dans l'équation de
    continuité (~ r_AU |S|/d ~ h²/ν) est petit devant ceux de u, v (~ h) sur les maillages
    fins ou visqueux, le seuil déclenche alors des permutations qui détruisent l'ordre
    (convection naturelle Ra = 10⁶, 96² : 17 millions de termes et 4.2 s au lieu de 2.8
    millions et 0.17 s). Les petites erreurs d'arrondi sont corrigées par GMRES. En cas de
    pivot nul ou d'échec (pivot=True) : pivotage partiel classique (COLAMD)."""
    if pivot:
        opts = dict(permc_spec="COLAMD")
    else:
        opts = dict(permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                    options=dict(SymmetricMode=True))
    A = A.tocsc()
    try:
        if exact:
            return spla.splu(A, **opts).solve
        return spla.spilu(A, drop_tol=1e-5, fill_factor=20, **opts).solve
    except RuntimeError:                               # pivot nul
        if pivot:
            raise
        return _factor(A, exact, pivot=True)


def _gmres(A, b, x0, M, atol, m):
    """GMRES préconditionné à droite, sans redémarrage, au plus m itérations : la norme
    minimisée est celle du vrai résidu ‖b − A x‖₂. Renvoie (x, itérations, convergé)."""
    r = b - A @ x0
    beta = float(np.linalg.norm(r))
    if beta <= atol:
        return x0, 0, True
    Q = np.empty((m + 1, b.size))
    Z = np.empty((m, b.size))
    H = np.zeros((m + 1, m))
    cs, sn, gv = np.zeros(m), np.zeros(m), np.zeros(m + 1)
    Q[0], gv[0] = r / beta, beta
    k, ok = 0, False
    for j in range(m):
        Z[j] = M(Q[j])
        v = A @ Z[j]
        for i in range(j + 1):                      # Gram-Schmidt modifié
            H[i, j] = v @ Q[i]
            v -= H[i, j] * Q[i]
        h = float(np.linalg.norm(v))
        if h > 0.0:
            Q[j + 1] = v / h
        for i in range(j):                          # rotations de Givens précédentes
            H[i, j], H[i + 1, j] = (cs[i] * H[i, j] + sn[i] * H[i + 1, j],
                                    -sn[i] * H[i, j] + cs[i] * H[i + 1, j])
        d = float(np.hypot(H[j, j], h))
        cs[j], sn[j] = H[j, j] / d, h / d
        H[j, j] = d
        gv[j + 1], gv[j] = -sn[j] * gv[j], cs[j] * gv[j]
        k = j + 1
        if abs(gv[k]) <= atol or h == 0.0:
            ok = True
            break
    y = solve_triangular(H[:k, :k], gv[:k])
    return x0 + y @ Z[:k], k, ok


def _solve(s, A, rhs, x0, r0):
    """Résolution du système couplé. La factorisation d'une itération sert de
    préconditionneur GMRES aux suivantes tant qu'elle reste efficace (≤ REFACTOR_AFTER
    itérations) : d'une itération couplée à l'autre la matrice change peu, et une
    descente-remontée coûte ~20 fois moins qu'une factorisation."""
    n = A.shape[0]
    st = getattr(s, "_coupled_lin", None)
    if st is None or st["n"] != n:
        st = s._coupled_lin = {"n": n, "M": None, "factorizations": 0, "krylov": 0}
    atol = INNER_TOL * float(np.linalg.norm(r0))
    if atol == 0.0:
        return x0
    if st["M"] is not None:
        x, k, ok = _gmres(A, rhs, x0, st["M"], atol, KRYLOV_MAX)
        st["krylov"] += k
        if ok:
            if k > REFACTOR_AFTER:
                st["M"] = None
            return x
    exact = n <= 3 * DIRECT_MAX
    st["M"] = _factor(A, exact, st.get("pivot", False))
    st["factorizations"] += 1
    x, k, ok = _gmres(A, rhs, x0, st["M"], atol, KRYLOV_MAX)
    st["krylov"] += k
    if not ok and not st.get("pivot", False):
        # factorisation sans pivotage trop imprécise : pivotage partiel désormais
        st["pivot"] = True
        st["M"] = _factor(A, exact, True)
        st["factorizations"] += 1
        x, k, ok = _gmres(A, rhs, x0, st["M"], atol, KRYLOV_MAX)
        st["krylov"] += k
    if not ok:
        st["M"] = None
    return x


def coupled_step(s, eqs, relax=1.0, a0=0.0):
    """Une itération couplée ; `eqs` = équations de quantité de mouvement (sans gradient
    de pression) déjà relaxées par `relax` (ou avec le terme de pseudo-temps a0·V). Met à
    jour s.U, s.p, s.F_i, s.F_b.

    Les flux de Rhie-Chow utilisent la diagonale PHYSIQUE (a_P · relax − a0·V) : la
    solution convergée ne dépend ni du facteur de sous-relaxation (correction de Majumdar
    1988) ni du pas de pseudo-temps, contrairement au SIMPLE(C) segmenté qui utilise la
    diagonale relaxée."""
    if s.backend.is_gpu:
        raise ValueError("algorithm = coupled : CPU seulement.")
    fvm = s.fvm
    nc, nb = fvm.nc, fvm.nb
    P, N, Pb = np.asarray(fvm.P), np.asarray(fvm.N), np.asarray(fvm.Pb)
    w, g = fvm.w, fvm.g
    Si, Sb, magSb, V = fvm.Si, fvm.Sb, fvm.magSb, fvm.V
    U, p = s.U, s.p

    diag_u, diag_v = eqs[0][0], eqs[1][0]
    aP = 0.5 * (diag_u + diag_v) * relax - a0 * V
    rAU = V / aP
    rA_f = fvm.interp(rAU)

    # inconnues entrelacées (u, v, p de la cellule i : 3i, 3i+1, 3i+2) : les couplages
    # restent locaux et la factorisation LU remplit ~2 fois moins qu'en blocs [u; v; p]
    rows, cols, vals = [], [], []
    rhs = np.zeros(3 * nc)
    cells = np.arange(nc)
    iP, iN = 3 * P, 3 * N

    # --- blocs quantité de mouvement A_c
    for c, (diag, up, lo, b) in enumerate(eqs):
        rows += [3 * cells + c, iP + c, iN + c]
        cols += [3 * cells + c, iN + c, iP + c]
        vals += [diag, up, lo]
        rhs[c::3] = b

    # --- blocs G (gradient de pression × V) : Σ p_f S_f, p_f = w p_P + (1 − w) p_N
    pa, pbeta, pg, pd = s.pressure_bc()
    for c in range(2):
        Sc = Si[:, c]
        rows += [iP + c, iP + c, iN + c, iN + c]
        cols += [iP + 2, iN + 2, iP + 2, iN + 2]
        vals += [w * Sc, (1 - w) * Sc, -w * Sc, -(1 - w) * Sc]
        if nb:
            rows.append(3 * Pb + c)
            cols.append(3 * Pb + 2)
            vals.append(pa * Sb[:, c])
            rhs[c::3] -= np.bincount(Pb, pbeta * Sb[:, c], nc)
        if fvm.axisymmetric:                           # faces latérales du secteur
            rows.append(3 * cells + c)
            cols.append(3 * cells + 2)
            vals.append(-V * fvm._side[:, c])

    # --- continuité : Σ flux sortants = 0
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
        rows += [iP + 2, iP + 2, iN + 2, iN + 2]
        cols += [iP + c, iN + c, iP + c, iN + c]
        vals += [w * Sc, (1 - w) * Sc, -w * Sc, -(1 - w) * Sc]
    rows += [iP + 2, iP + 2, iN + 2, iN + 2]
    cols += [iP + 2, iN + 2, iP + 2, iN + 2]
    vals += [D, -D, -D, D]
    rhs[2::3] -= np.bincount(P, E_i - nf, nc) - np.bincount(N, E_i - nf, nc)
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
            rows.append(3 * Pb + 2)
            cols.append(3 * Pb + c)
            vals.append(np.where(zg, Sb[:, c], 0.0))
        Fb_const += np.where(zg, rAb * np.sum(gradp[Pb] * Sb, axis=1), 0.0)
        if s._has_force():
            fb = np.sum(s._boundary_force() * Sb, axis=1)
            Fb_const += np.where(np.asarray(s.kindP) == 1, rAb * fb, 0.0)
        rows.append(3 * Pb + 2)
        cols.append(3 * Pb + 2)
        vals.append(-rAb * magSb * pg)
        Fb_const -= rAb * magSb * pd
        rhs[2::3] -= np.bincount(Pb, Fb_const, nc)
    if not s.has_fixed_p:
        # niveau de pression indéterminé : diagonale de p renforcée dans la 1re cellule
        # (comme le SIMPLE segmenté) ; sans effet à convergence (p₀ ne bouge plus)
        pin = float(np.sum(D[(P == 0) | (N == 0)])) or 1.0
        rows.append(np.array([2]))
        cols.append(np.array([2]))
        vals.append(np.array([pin]))
        rhs[2] += pin * p[0]
    A = _coo(rows, cols, vals, (3 * nc, 3 * nc)).tocsr()
    x0 = np.column_stack([U[:, 0], U[:, 1], p]).ravel()
    r0 = rhs - A @ x0
    # résidu de continuité de l'itéré courant (flux de Rhie-Chow avec les coefficients à
    # jour), normalisé comme « continuity » : Σ|Σ_f F_f| / (U_ref · Σ|S_f|)
    s.residuals_now["p"] = float(np.sum(np.abs(r0[2::3]))) / (
        s.u_scale * fvm.total_area + 1e-300)
    x = _solve(s, A, rhs, x0, r0)
    X = x.reshape(nc, 3)
    s.U = np.ascontiguousarray(X[:, :2])
    s.p = np.ascontiguousarray(X[:, 2])
    # flux conservatifs (même expression que la ligne de continuité)
    Uf = fvm.interp(s.U)
    s.F_i = np.sum(Uf * Si, axis=1) + E_i - D * (s.p[N] - s.p[P]) - nf
    if nb:
        s.F_b = (Fb_const + np.where(zg, np.sum(s.U[Pb] * Sb, axis=1), 0.0)
                 - rAb * magSb * pg * s.p[Pb])
