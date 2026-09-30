"""Marche implicite en pseudo-temps pour le solveur compressible (stationnaire).

Euler implicite linéarisé à pas de temps local (Blazek § 6.2) :

    (V_i/Δt_i I + J) ΔQ = −R(Q),      V_i/Δt_i = Λ_i / CFL,

- R : résidu COMPLET (ordre 2, flux de Roe/HLLC, visqueux) — la solution convergée est
  celle du schéma explicite ;
- J : jacobienne APPROCHÉE d'ordre 1 : flux de Rusanov (Lax-Friedrichs local)
  F = ½(F_L + F_R) − ½ λ (Q_R − Q_L), λ = max(|u·n| + c) des deux côtés, d'où les blocs
  4×4 ½(A_L + λI)·S et ½(A_R − λI)·S (A = ∂(F·n)/∂Q, Blazek annexe A.9) ; conditions aux
  limites : état fantôme miroir des parois dérivé exactement (G = M Q_P), extrapolation
  exacte, autres conditions à état extérieur figé ; visqueux : diffusion scalaire
  max(4/3, γ/Pr) μ/ρ · S/d sur les équations de quantité de mouvement et d'énergie ;
- résolution approchée : mise à l'échelle par l'inverse des blocs diagonaux
  (D⁻¹A = I + L' + U', L'/U' : blocs strictement inférieurs / supérieurs) puis
  `linear_sweeps` balayages de Gauss-Seidel symétrique par blocs (avant / arrière) —
  l'équivalent « matriciel » du LU-SGS (Yoon & Jameson 1988), les solves triangulaires
  étant faits par SuperLU (scipy.sparse.linalg.spsolve_triangular, code compilé) ;
- sécurité : sous-relaxation locale de ΔQ (variation relative de ρ et ρE ≤ 20 % par
  itération, comme la relaxation par point de SU2) ; CFL multiplié par `cfl_growth` à
  chaque itération jusqu'à `cfl_max`, divisé par 2 si la mise à jour est trop limitée.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve_triangular

from .compressible import _EXTRAP, _NOSLIP, _SLIP

_E = np.diag([0.0, 1.0, 1.0, 1.0])            # équations touchées par la diffusion
_I4 = np.eye(4)


def flux_jacobian(W, nx, ny, g):
    """A = ∂(F·n)/∂Q (n, 4, 4) aux états primitifs W (4, n) ; n non unitaire accepté."""
    r, u, v, p = W
    g1 = g - 1.0
    qn = u * nx + v * ny
    phi2 = 0.5 * g1 * (u * u + v * v)
    H = g / g1 * p / r + 0.5 * (u * u + v * v)
    n = len(r)
    A = np.empty((n, 4, 4))
    A[:, 0, 0] = 0.0
    A[:, 0, 1] = nx
    A[:, 0, 2] = ny
    A[:, 0, 3] = 0.0
    A[:, 1, 0] = phi2 * nx - u * qn
    A[:, 1, 1] = qn - (g - 2.0) * u * nx
    A[:, 1, 2] = u * ny - g1 * v * nx
    A[:, 1, 3] = g1 * nx
    A[:, 2, 0] = phi2 * ny - v * qn
    A[:, 2, 1] = v * nx - g1 * u * ny
    A[:, 2, 2] = qn - (g - 2.0) * v * ny
    A[:, 2, 3] = g1 * ny
    A[:, 3, 0] = qn * (phi2 - H)
    A[:, 3, 1] = H * nx - g1 * u * qn
    A[:, 3, 2] = H * ny - g1 * v * qn
    A[:, 3, 3] = g * qn
    return A


def _pattern(rowb, colb, nc):
    """Structure CSR (éléments) d'une matrice faite de blocs 4×4 (rowb, colb) + diagonale
    unité ; perm : position de chaque coefficient dans [blocs.ravel(), 1.0]."""
    nf = len(rowb)
    a = np.arange(4)
    r = (4 * rowb[:, None, None] + a[None, :, None] + 0 * a[None, None, :]).ravel()
    c = (4 * colb[:, None, None] + 0 * a[None, :, None] + a[None, None, :]).ravel()
    r = np.concatenate([r, np.arange(4 * nc)])
    c = np.concatenate([c, np.arange(4 * nc)])
    src = np.concatenate([np.arange(16 * nf), np.full(4 * nc, 16 * nf)])
    order = np.lexsort((c, r))
    indptr = np.concatenate([[0], np.cumsum(np.bincount(r, minlength=4 * nc))])
    M = sp.csr_matrix((np.ones(len(r)), c[order].astype(np.int32),
                       indptr.astype(np.int32)), shape=(4 * nc, 4 * nc))
    M.sum_duplicates()                       # (blocs en double : maillages dégénérés)
    canonical = M.nnz == len(r)
    return M, src[order], canonical


class ImplicitStepper:
    """Un pas implicite (voir l'en-tête du module) pour un CompressibleSolver2D."""

    def __init__(self, solver):
        s = solver
        self.s = s
        nc = s.nc
        keep = s.P != s.N                        # face périodique d'une cellule sur elle-même
        self.f = np.nonzero(keep)[0]
        P, N = s.P[self.f], s.N[self.f]
        self.fP, self.fN = P, N
        self.low = P > N                         # bloc (P, N) sous la diagonale
        lo, hi = np.minimum(P, N), np.maximum(P, N)
        self.Lm, self.Lperm, cL = _pattern(hi, lo, nc)
        self.Um, self.Uperm, cU = _pattern(lo, hi, nc)
        self.canonical = cL and cU
        nf = len(self.f)
        self.IP = sp.csr_matrix((np.ones(nf), (P, np.arange(nf))), shape=(nc, nf))
        self.IN = sp.csr_matrix((np.ones(nf), (N, np.arange(nf))), shape=(nc, nf))
        self.IB = sp.csr_matrix((np.ones(s.nb), (s.Pb, np.arange(s.nb))), shape=(nc, s.nb))
        # parois : matrices M = ∂G/∂Q_P de l'état fantôme
        nb = s.nb
        M = np.zeros((nb, 4, 4))
        M[:, 0, 0] = M[:, 3, 3] = 1.0
        nx, ny = s.nbx, s.nby
        sl = s.kind == _SLIP
        M[sl, 1, 1] = 1.0 - 2.0 * nx[sl] ** 2
        M[sl, 1, 2] = M[sl, 2, 1] = -2.0 * nx[sl] * ny[sl]
        M[sl, 2, 2] = 1.0 - 2.0 * ny[sl] ** 2
        ns = s.kind == _NOSLIP
        M[ns, 1, 1] = M[ns, 2, 2] = -1.0
        self.Mwall = M
        self.wall = sl | ns
        self.extrap = s.kind == _EXTRAP
        self.limited = 0.0

    # ------------------------------------------------------------------ jacobienne
    def _assemble(self, cfl):
        s = self.s
        g = s.gas.gamma
        W = s.primitive()
        f, P, N = self.f, self.fP, self.fN
        nx, ny, S = s.nx[f], s.ny[f], s.magS[f]
        WP, WN = np.take(W, P, axis=1), np.take(W, N, axis=1)
        AP = flux_jacobian(WP, nx, ny, g)
        AN = flux_jacobian(WN, nx, ny, g)
        c = np.sqrt(g * W[3] / W[0])
        qn = W[1] * 0.0
        lamP = np.abs(WP[1] * nx + WP[2] * ny) + c[P]
        lamN = np.abs(WN[1] * nx + WN[2] * ny) + c[N]
        lam = np.maximum(lamP, lamN)
        del qn
        hS = 0.5 * S
        Bpp = (AP + lam[:, None, None] * _I4) * hS[:, None, None]      # ∂R_P/∂Q_P
        Bpn = (AN - lam[:, None, None] * _I4) * hS[:, None, None]      # ∂R_P/∂Q_N
        if s.gas.viscous:
            T = W[3] / (W[0] * s.gas.R)
            nu = s.gas.mu_of(T) / W[0]
            kv = max(4.0 / 3.0, g / s.gas.Pr)
            vf = kv * 0.5 * (nu[P] + nu[N]) * S / s.dist[f]
            Bpp = Bpp + vf[:, None, None] * _E
            Bpn = Bpn - vf[:, None, None] * _E
        # blocs diagonaux : Σ ∂R_P/∂Q_P (owner) + Σ ∂R_N/∂Q_N = −Bpn (neighbour)
        D = np.asarray(self.IP @ Bpp.reshape(-1, 16) - self.IN @ Bpn.reshape(-1, 16))
        # frontières
        Pb = s.Pb
        Wc = np.take(W, Pb, axis=1)
        nbx, nby, Sb = s.nbx, s.nby, s.magSb
        Ab = flux_jacobian(Wc, nbx, nby, g)
        lb = np.abs(Wc[1] * nbx + Wc[2] * nby) + c[Pb]
        Bb = 0.5 * (Ab + lb[:, None, None] * _I4)
        w = self.wall
        if w.any():
            G = s.ghost(Wc)[:, w]
            AG = flux_jacobian(G, nbx[w], nby[w], g)
            Bb[w] += 0.5 * np.matmul(AG - lb[w, None, None] * _I4, self.Mwall[w])
        e = self.extrap
        if e.any():
            Bb[e] = Ab[e]
        Bb *= Sb[:, None, None]
        if s.gas.viscous and len(s._noslip):
            ns = s._noslip
            T = W[3, Pb[ns]] / (W[0, Pb[ns]] * s.gas.R)
            vb = kv * s.gas.mu_of(T) / W[0, Pb[ns]] * Sb[ns] / s.distb[ns]
            Bb[ns] += vb[:, None, None] * _E
        D += np.asarray(self.IB @ Bb.reshape(-1, 16))
        D = D.reshape(-1, 4, 4)
        D += (s.spectral_radius(s.Q) / cfl)[:, None, None] * _I4
        Dinv = np.linalg.inv(D)
        # blocs hors diagonale mis à l'échelle : ligne P (colonne N) et ligne N (colonne P)
        Opn = np.matmul(Dinv[P], Bpn)
        Onp = np.matmul(Dinv[N], -Bpp)
        low = self.low
        Lb = np.where(low[:, None, None], Opn, Onp)
        Ub = np.where(low[:, None, None], Onp, Opn)
        self.Lm.data[:] = np.append(Lb.ravel(), 1.0)[self.Lperm]
        self.Um.data[:] = np.append(Ub.ravel(), 1.0)[self.Uperm]
        if self.canonical:
            self.Lm.has_canonical_format = True
            self.Um.has_canonical_format = True
        return Dinv

    # ------------------------------------------------------------------ pas
    def step(self, R, cfl, order=2) -> bool:
        """Met à jour s.Q ; renvoie False si la mise à jour a dû être fortement limitée
        (le CFL doit alors baisser)."""
        s = self.s
        Dinv = self._assemble(cfl)
        b = np.einsum("cij,cj->ci", Dinv, -R.T).ravel()
        sweeps = max(int(getattr(s.settings, "linear_sweeps", 2)), 1)
        x = np.zeros_like(b)
        for k in range(sweeps):
            rhs = b if k == 0 else b - (self.Um @ x - x)
            x = spsolve_triangular(self.Lm, rhs, lower=True, unit_diagonal=True)
            rhs = b - (self.Lm @ x - x)
            x = spsolve_triangular(self.Um, rhs, lower=False, unit_diagonal=True)
        dQ = x.reshape(-1, 4).T
        Q = s.Q
        # sous-relaxation locale : |Δρ| ≤ 0.2 ρ, |Δ(ρE)| ≤ 0.2 ρE
        ratio = np.maximum(np.abs(dQ[0]) / Q[0], np.abs(dQ[3]) / Q[3])
        alpha = np.minimum(1.0, 0.2 / np.maximum(ratio, 1e-300))
        for _ in range(6):
            Qn = Q + alpha * dQ
            Wn = s.primitive(Qn)
            bad = (Wn[0] <= 0) | (Wn[3] <= 0) | ~np.isfinite(Wn).all(axis=0)
            if not bad.any():
                break
            alpha[bad] *= 0.25
        else:
            raise FloatingPointError("Pas implicite : état non physique malgré la "
                                     "sous-relaxation (réduire cfl / cfl_max).")
        s.Q = Qn
        self.limited = float(np.mean(alpha < 1.0))
        return self.limited < 0.05
