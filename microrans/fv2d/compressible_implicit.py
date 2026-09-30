"""Marche implicite en pseudo-temps pour le solveur compressible (stationnaire).

Euler implicite linéarisé à pas de temps local (Blazek § 6.2) :

    (V_i/Δt_i I + J) ΔQ = −R(Q),      V_i/Δt_i = Λ_i / CFL,

- R : résidu COMPLET (ordre 2, flux de Roe/HLLC, visqueux) — la solution convergée est
  celle du schéma explicite ;
- J : jacobienne APPROCHÉE d'ordre 1 : F = ½(F_L + F_R) − ½ |Ã| (Q_R − Q_L), d'où les
  blocs 4×4 ½(A_L + |Ã|)·S et ½(A_R − |Ã|)·S (A = ∂(F·n)/∂Q, Blazek annexe A.9) ; |Ã| :
  matrice de dissipation de Roe (moyennes de Roe des cellules, correction d'entropie
  comprise), obtenue colonne par colonne avec la décomposition en ondes du flux de Roe —
  jacobian = "roe" (défaut). jacobian = "rusanov" : |Ã| ≈ λ I, λ = max(|u·n| + c), plus
  simple mais trop dissipatif à bas Mach (ondes de convection ralenties d'un facteur
  |u|/(|u| + c) : 6 fois à M = 0.2). Conditions aux limites : état fantôme miroir des
  parois dérivé exactement (G = M Q_P), extrapolation exacte, autres conditions à état
  extérieur figé ; visqueux : approximation « couche mince » (dérivées selon la normale
  seulement, Blazek § A.10) τ·n ≈ μ/d (Δu + ⅓ n (n·Δu)), q·n ≈ −k ΔT/d, linéarisée
  exactement par rapport à Q via ∂(u, v, T)/∂Q (paroi adiabatique : pas de flux de chaleur,
  isotherme : ΔT = T_w − T_P) ;
- résolution approchée : mise à l'échelle par l'inverse des blocs diagonaux
  (D⁻¹A = I + L' + U', L'/U' : blocs strictement inférieurs / supérieurs), puis
    * linear_solver = "gmres" (défaut) : GMRES (Saad & Schultz 1986) préconditionné par un
      balayage de Gauss-Seidel symétrique par blocs, (I + L')(I + U'), au plus
      `linear_iter` itérations (défaut 20) ou tolérance relative `linear_tol` (0.05) —
      comme le FGMRES + ILU de SU2 ;
    * linear_solver = "sgs" : `linear_sweeps` balayages de Gauss-Seidel symétrique par
      blocs (l'équivalent « matriciel » du LU-SGS, Yoon & Jameson 1988) ; moins cher par
      itération mais inefficace aux grands CFL (mesuré : 2 balayages ne réduisent le résidu
      linéaire que de 0.94 à 0.61 à CFL 1000 sur la plaque plane) ;
  descentes / remontées triangulaires compilées (SuperLU) ;
- sécurité : sous-relaxation locale de ΔQ (variation relative de ρ et ρE ≤ 20 % par
  itération, comme la relaxation par point de SU2) ; CFL multiplié par `cfl_growth` à
  chaque itération jusqu'à `cfl_max`, divisé par 2 si la mise à jour est trop limitée.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve_triangular

try:                                          # descente / remontée compilées de SuperLU
    from scipy.sparse.linalg._dsolve import _superlu
except ImportError:                           # pragma: no cover
    _superlu = None

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


def roe_dissipation(WL, WR, nx, ny, g, eps: float = 0.1):
    """Matrice |Ã| (n, 4, 4) de la dissipation de Roe entre WL et WR (4, n) : la
    colonne k est la dissipation de roe_flux pour ΔQ = e_k, linéarisée à l'état moyen de
    Roe (mêmes valeurs propres, même correction d'entropie de Harten)."""
    rL, uL, vL, pL = WL
    rR, uR, vR, pR = WR
    sL, sR = np.sqrt(rL), np.sqrt(rR)
    iw = 1.0 / (sL + sR)
    HL = g / (g - 1.0) * pL / rL + 0.5 * (uL * uL + vL * vL)
    HR = g / (g - 1.0) * pR / rR + 0.5 * (uR * uR + vR * vR)
    r = sL * sR
    u = (sL * uL + sR * uR) * iw
    v = (sL * vL + sR * vR) * iw
    H = (sL * HL + sR * HR) * iw
    k = 0.5 * (u * u + v * v)
    c2 = np.maximum((g - 1.0) * (H - k), 1e-300)
    c = np.sqrt(c2)
    q = u * nx + v * ny
    l1, l2, l3 = np.abs(q - c), np.abs(q), np.abs(q + c)
    if eps > 0.0:
        d = eps * c
        h1, h3 = np.maximum(d - l1, 0.0), np.maximum(d - l3, 0.0)
        l1 = l1 + 0.5 * h1 * h1 / d
        l3 = l3 + 0.5 * h3 * h3 / d
    g1 = g - 1.0
    n = len(r)
    D = np.empty((n, 4, 4))
    for col in range(4):
        # ΔW = (∂W/∂Q) e_col à l'état de Roe
        if col == 0:
            dr, du, dv, dp = 1.0, -u / r, -v / r, g1 * k
        elif col == 1:
            dr, du, dv, dp = 0.0, 1.0 / r, 0.0, -g1 * u
        elif col == 2:
            dr, du, dv, dp = 0.0, 0.0, 1.0 / r, -g1 * v
        else:
            dr, du, dv, dp = 0.0, 0.0, 0.0, g1
        dq = du * nx + dv * ny
        a1 = l1 * (dp - r * c * dq) / (2.0 * c2)
        a2 = l2 * (dr - dp / c2)
        a3 = l3 * (dp + r * c * dq) / (2.0 * c2)
        a4 = l2 * r
        D[:, 0, col] = a1 + a2 + a3
        D[:, 1, col] = a1 * (u - c * nx) + a2 * u + a3 * (u + c * nx) + a4 * (du - dq * nx)
        D[:, 2, col] = a1 * (v - c * ny) + a2 * v + a3 * (v + c * ny) + a4 * (dv - dq * ny)
        D[:, 3, col] = (a1 * (H - c * q) + a2 * k + a3 * (H + c * q)
                        + a4 * (u * du + v * dv - q * dq))
    return D


def _dvars(W, R, g):
    """∂(u, v, T)/∂Q (n, 3, 4) aux états primitifs W (4, n)."""
    r, u, v, p = W
    E = p / ((g - 1.0) * r) + 0.5 * (u * u + v * v)
    n = len(r)
    T = np.zeros((n, 3, 4))
    T[:, 0, 0], T[:, 0, 1] = -u / r, 1.0 / r
    T[:, 1, 0], T[:, 1, 2] = -v / r, 1.0 / r
    a = (g - 1.0) / (R * r)
    T[:, 2, 0] = a * (u * u + v * v - E)
    T[:, 2, 1] = -a * u
    T[:, 2, 2] = -a * v
    T[:, 2, 3] = a
    return T


def _thin_layer(mu, k, u, v, nx, ny, coef, heat=None):
    """B (n, 4, 3) : flux visqueux « couche mince » = coef · B · Δ(u, v, T),
    coef = S/d ; heat : masque des faces sans flux de chaleur (paroi adiabatique)."""
    n = len(mu)
    B = np.zeros((n, 4, 3))
    B[:, 1, 0] = mu * (1.0 + nx * nx / 3.0)
    B[:, 1, 1] = B[:, 2, 0] = mu * nx * ny / 3.0
    B[:, 2, 1] = mu * (1.0 + ny * ny / 3.0)
    B[:, 3, 0] = u * B[:, 1, 0] + v * B[:, 2, 0]
    B[:, 3, 1] = u * B[:, 1, 1] + v * B[:, 2, 1]
    B[:, 3, 2] = k if heat is None else np.where(heat, 0.0, k)
    return B * coef[:, None, None]


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


class _Triangular:
    """Résolution (I + T) x = b, T strictement triangulaire (CSR, diagonale unité stockée).

    Appel direct de la descente / remontée compilée de SuperLU (celle qu'utilise
    scipy.sparse.linalg.spsolve_triangular, sans ses copies et conversions à chaque appel :
    3 à 10 fois plus rapide ici) ; vérifiée à la construction contre spsolve_triangular,
    qui sert de repli si l'interface interne de SciPy change."""

    def __init__(self, M, lower: bool):
        self.M, self.lower = M, lower
        n = M.shape[0]
        self.n = n
        rows = np.repeat(np.arange(n), np.diff(M.indptr))
        self.dpos = np.nonzero(M.indices == rows)[0]
        self.ind = M.indices.astype(np.intc)
        self.ptr = M.indptr.astype(np.intc)
        self.eye = (np.ones(n), np.arange(n, dtype=np.intc), np.arange(n + 1, dtype=np.intc))
        self.empty = (np.zeros(0), np.zeros(0, dtype=np.intc), np.zeros(n + 1, dtype=np.intc))
        self.fast = _superlu is not None and len(self.dpos) == n
        if self.fast:
            saved = M.data.copy()
            rng = np.random.default_rng(1)
            M.data[:] = 0.1 * rng.standard_normal(M.nnz) / 4.0
            M.data[self.dpos] = 1.0
            b = rng.standard_normal(n)
            try:
                ok = np.allclose(self._fast(b), spsolve_triangular(
                    M, b, lower=lower, unit_diagonal=True), rtol=1e-10, atol=1e-12)
            except Exception:                              # noqa: BLE001
                ok = False
            self.fast = bool(ok)
            M.data[:] = saved

    def _fast(self, b):
        M, n = self.M, self.n
        b = np.array(b, dtype=float, copy=True)
        if self.lower:
            # Mᵀ (CSC) triangulaire supérieure : L = I, U = Mᵀ de diagonale nulle
            data = M.data.copy()
            data[self.dpos] = 0.0
            x, info = _superlu.gstrs("T", n, n, *self.eye, n, M.nnz, data, self.ind,
                                     self.ptr, b)
        else:
            # Mᵀ (CSC) triangulaire inférieure à diagonale unité : L = Mᵀ, U vide
            x, info = _superlu.gstrs("T", n, M.nnz, M.data, self.ind, self.ptr,
                                     n, 0, *self.empty, b)
        if info:
            raise np.linalg.LinAlgError("Résolution triangulaire : matrice singulière.")
        return x

    def solve(self, b):
        if self.fast:
            return self._fast(b)
        return spsolve_triangular(self.M, b, lower=self.lower, unit_diagonal=True)


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
        self.any_low = bool(self.low.any())
        lo, hi = np.minimum(P, N), np.maximum(P, N)
        self.Lm, self.Lperm, cL = _pattern(hi, lo, nc)
        self.Um, self.Uperm, cU = _pattern(lo, hi, nc)
        self.canonical = cL and cU
        nf = len(self.f)
        self._buf = np.empty(16 * nf + 1)
        self._buf[-1] = 1.0
        self.Lsolve = _Triangular(self.Lm, True)
        self.Usolve = _Triangular(self.Um, False)
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
        self.linear_iterations = 0
        st = s.settings
        self.roe = str(getattr(st, "implicit_jacobian", "roe")).lower() == "roe"
        self.method = str(getattr(st, "linear_solver", "gmres")).lower()
        if self.method not in ("gmres", "sgs"):
            raise ValueError(f"linear_solver = '{self.method}' : gmres | sgs.")

    # ------------------------------------------------------------------ jacobienne
    def _assemble(self, cfl):
        s = self.s
        g = s.gas.gamma
        W = s.primitive()
        f, P, N = self.f, self.fP, self.fN
        nx, ny, S = s.nx[f], s.ny[f], s.magS[f]
        WP, WN = np.take(W, P, axis=1), np.take(W, N, axis=1)
        hS = (0.5 * S)[:, None, None]
        if self.roe:
            Dis = roe_dissipation(WP, WN, nx, ny, g, s.settings.entropy_fix)
        else:
            c = np.sqrt(g * W[3] / W[0])
            lam = np.maximum(np.abs(WP[1] * nx + WP[2] * ny) + c[P],
                             np.abs(WN[1] * nx + WN[2] * ny) + c[N])
            Dis = lam[:, None, None] * _I4
        Bpp = flux_jacobian(WP, nx, ny, g)
        Bpp += Dis
        Bpp *= hS                                 # ∂R_P/∂Q_P
        Bpn = flux_jacobian(WN, nx, ny, g)
        Bpn -= Dis
        Bpn *= hS                                 # ∂R_P/∂Q_N
        gas = s.gas
        if gas.viscous:
            # ∂R_P/∂Q_P += (S/d) B T_P ; ∂R_P/∂Q_N −= (S/d) B T_N
            w = s.w[f]
            Tc = W[3] / (W[0] * gas.R)
            mu_f = gas.mu_of(w * Tc[P] + (1 - w) * Tc[N])
            uf = w * WP[1] + (1 - w) * WN[1]
            vf = w * WP[2] + (1 - w) * WN[2]
            Bv = _thin_layer(mu_f, mu_f * gas.cp / gas.Pr, uf, vf, nx, ny, S / s.dist[f])
            Bpp += np.matmul(Bv, _dvars(WP, gas.R, g))
            Bpn -= np.matmul(Bv, _dvars(WN, gas.R, g))
        # blocs diagonaux : Σ ∂R_P/∂Q_P (owner) + Σ ∂R_N/∂Q_N = −Bpn (neighbour)
        D = np.asarray(self.IP @ Bpp.reshape(-1, 16) - self.IN @ Bpn.reshape(-1, 16))
        # frontières
        Pb = s.Pb
        Wc = np.take(W, Pb, axis=1)
        nbx, nby, Sb = s.nbx, s.nby, s.magSb
        Ab = flux_jacobian(Wc, nbx, nby, g)
        G = s.ghost(Wc)
        if self.roe:
            Disb = roe_dissipation(Wc, G, nbx, nby, g, s.settings.entropy_fix)
        else:
            c = np.sqrt(g * W[3] / W[0])
            Disb = (np.abs(Wc[1] * nbx + Wc[2] * nby) + c[Pb])[:, None, None] * _I4
        Bb = 0.5 * (Ab + Disb)
        w = self.wall
        if w.any():
            AG = flux_jacobian(G[:, w], nbx[w], nby[w], g)
            Bb[w] += 0.5 * np.matmul(AG - Disb[w], self.Mwall[w])
        e = self.extrap
        if e.any():
            Bb[e] = Ab[e]
        Bb *= Sb[:, None, None]
        if gas.viscous and len(s._noslip):
            # parois adhérentes : Δ(u, v, T) = (u_w − u_P, …) ; travail u_w·τ nul
            ns = s._noslip
            Wn = Wc[:, ns]
            mu_b = gas.mu_of(Wn[3] / (Wn[0] * gas.R))
            z = np.zeros(len(ns))
            Bv = _thin_layer(mu_b, mu_b * gas.cp / gas.Pr, z, z, nbx[ns], nby[ns],
                             Sb[ns] / s.distb[ns], heat=~np.isfinite(s.bc_Tw[ns]))
            Bb[ns] += np.matmul(Bv, _dvars(Wn, gas.R, g))
        D += np.asarray(self.IB @ Bb.reshape(-1, 16))
        D = D.reshape(-1, 4, 4)
        dt_term = s.spectral_radius(s.Q) / cfl
        for k in range(4):
            D[:, k, k] += dt_term
        Dinv = np.linalg.inv(D)
        # blocs hors diagonale mis à l'échelle : ligne P (colonne N) et ligne N (colonne P)
        Opn = np.matmul(Dinv[P], Bpn)
        Onp = np.matmul(Dinv[N], Bpp)
        Onp *= -1.0
        if self.any_low:
            low = self.low
            Lb, Ub = Onp.copy(), Opn.copy()
            Lb[low], Ub[low] = Opn[low], Onp[low]
        else:
            Lb, Ub = Onp, Opn
        buf = self._buf
        buf[:-1] = Lb.reshape(-1)
        np.take(buf, self.Lperm, out=self.Lm.data)
        buf[:-1] = Ub.reshape(-1)
        np.take(buf, self.Uperm, out=self.Um.data)
        if self.canonical:
            self.Lm.has_canonical_format = True
            self.Um.has_canonical_format = True
        return Dinv

    # ------------------------------------------------------------------ résolution
    def _sgs(self, b, sweeps):
        x = None
        for k in range(sweeps):
            rhs = b if x is None else b - (self.Um @ x - x)
            x = self.Lsolve.solve(rhs)
            rhs = b - (self.Lm @ x - x)
            x = self.Usolve.solve(rhs)
        return x

    def _gmres(self, b, restart, rtol):
        from scipy.sparse.linalg import LinearOperator, gmres
        n = len(b)
        Lm, Um = self.Lm, self.Um
        A = LinearOperator((n, n), matvec=lambda x: Lm @ x + Um @ x - x, dtype=float)
        Mp = LinearOperator((n, n), matvec=lambda r: self.Usolve.solve(self.Lsolve.solve(r)),
                            dtype=float)
        count = [0]

        def cb(_):
            count[0] += 1
        x, _ = gmres(A, b, rtol=rtol, restart=restart, maxiter=1, M=Mp, callback=cb,
                     callback_type="pr_norm")
        self.linear_iterations = count[0]
        return x

    # ------------------------------------------------------------------ pas
    def step(self, R, cfl, order=2) -> bool:
        """Met à jour s.Q ; renvoie False si la mise à jour a dû être fortement limitée
        (le CFL doit alors baisser)."""
        s = self.s
        st = s.settings
        Dinv = self._assemble(cfl)
        b = np.einsum("cij,cj->ci", Dinv, -R.T).ravel()
        if self.method == "gmres":
            x = self._gmres(b, max(int(getattr(st, "linear_iter", 20)), 1),
                            float(getattr(st, "linear_tol", 0.05)))
        else:
            x = self._sgs(b, max(int(getattr(st, "linear_sweeps", 2)), 1))
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
