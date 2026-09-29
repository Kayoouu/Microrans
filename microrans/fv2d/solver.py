"""Solveur incompressible 2D RANS / URANS, volumes finis colocalisés.

Algorithmes (à la manière d'OpenFOAM) :
- stationnaire : SIMPLE ou SIMPLEC (« consistent »), sous-relaxations U, p, turbulence ;
- instationnaire : PIMPLE (boucles externes + corrections PISO), schémas Euler implicite ou
  « backward » (BDF2), correction ddtCorr du flux.
Couplage vitesse-pression : forme OpenFOAM de l'interpolation de Rhie-Chow
(HbyA = H/a_P, flux φ = φ(HbyA) − (1/a_P)_f ∇p·S).

Conditions aux limites physiques par patch (type, à la manière de SU2) :
  wall      paroi adhérente (option U = [ux, uy] pour une paroi mobile)
  inlet     vitesse imposée U (constante ou expressions en x, y), turbulence amont
  outlet    pression imposée p (défaut 0), gradient nul pour U et la turbulence
  symmetry  glissement (vitesse normale nulle)
  farfield  champ lointain : U∞ imposée en entrée, p∞ en sortie (selon le signe de U∞·n)
Les patches périodiques sont gérés par le maillage (pas de condition à donner).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from ..mesh2d.mesh import Mesh2D
from ..models import canonical_name, get_model
from .fvm import FVM, normalized_residual, solve_linear

BC_TYPES = ("wall", "inlet", "outlet", "symmetry", "farfield")


@dataclass
class Settings:
    """Paramètres numériques (valeurs par défaut raisonnables pour débuter)."""
    algorithm: str = "SIMPLEC"            # SIMPLE | SIMPLEC (stationnaire)
    relax_U: float | None = None          # None = auto : 0.9 (maillage ~orthogonal) / 0.7
    relax_p: float = 1.0                  # SIMPLE : ~0.3 ; SIMPLEC : 1
    relax_turb: float = 0.8
    convection_U: str = "linearUpwind"    # upwind | linearUpwind
    convection_turb: str = "upwind"
    max_iter: int = 3000
    tol: float = 1e-5                     # résidus normalisés (OpenFOAM) de tous les champs
    solver_p: str = "direct"              # direct | bicgstab | amg (pyamg si installé)
    solver_U: str = "auto"                # auto = direct si < 40 000 cellules, sinon bicgstab
    solver_turb: str = "auto"
    # instationnaire
    time_scheme: str = "backward"         # euler | backward
    n_outer: int = 2
    n_corr: int = 2
    turbulence_every_outer: bool = False
    # maillages non orthogonaux
    n_nonorth: int = 1                    # boucles de correction non orthogonale (pression)
    nonorth_limit: float = 0.5            # limiteur ψ de la correction (1 = aucun)


def _eval_expr(expr, x, y):
    if isinstance(expr, (int, float)):
        return np.full_like(x, float(expr))
    ns = {"x": x, "y": y, "np": np, "pi": np.pi, "sqrt": np.sqrt, "sin": np.sin, "cos": np.cos,
          "exp": np.exp, "tanh": np.tanh, "abs": np.abs, "minimum": np.minimum,
          "maximum": np.maximum}
    val = eval(str(expr), {"__builtins__": {}}, ns)    # expression fournie par l'utilisateur
    return np.broadcast_to(np.asarray(val, dtype=float), x.shape).copy()


class FlowField2D:
    """Grandeurs du champ moyen utilisées par les modèles de turbulence."""

    def __init__(self, solver: "Solver2D", U, gradU):
        self.solver = solver
        self.U = U
        self.gradU = gradU               # gradU[:, i, j] = ∂u_i/∂x_j
        self._d2 = None

    @property
    def strain(self):
        g = self.gradU
        return np.sqrt(2 * g[:, 0, 0] ** 2 + 2 * g[:, 1, 1] ** 2 + (g[:, 0, 1] + g[:, 1, 0]) ** 2)

    @property
    def vorticity(self):
        g = self.gradU
        return np.abs(g[:, 1, 0] - g[:, 0, 1])

    @property
    def second_derivative_sq(self):
        if self._d2 is None:
            s = self.solver
            tot = np.zeros(s.fvm.nc)
            for i in range(2):
                for j in range(2):
                    comp = self.gradU[:, i, j]
                    gg = s.fvm.grad(comp, comp[s.fvm.Pb])
                    tot += np.sum(gg ** 2, axis=1)
            self._d2 = tot
        return self._d2


class Ops2D:
    """Gradients pour les modèles (valeur pariétale imposée, gradient nul ailleurs)."""

    def __init__(self, solver: "Solver2D"):
        self.s = solver

    def _grad(self, f, name):
        s = self.s
        fb = f[s.fvm.Pb].copy()
        if name is not None and np.any(s.is_wall):
            base = "k" if name == "sqrt_k" else name
            wv = s.model.wall_value(base, s.fvm.dperp[s.is_wall])
            fb[s.is_wall] = np.sqrt(np.maximum(wv, 0.0)) if name == "sqrt_k" else wv
        return s.fvm.grad(f, fb)

    def grad_sq(self, f, name=None):
        return np.sum(self._grad(f, name) ** 2, axis=1)

    def grad_dot(self, f, g, fname=None, gname=None):
        return np.sum(self._grad(f, fname) * self._grad(g, gname), axis=1)


class Step2D:
    """Résolution implicite d'une équation de transport de turbulence en 2D."""

    def __init__(self, solver: "Solver2D", a0: float, hist: dict, relax: float):
        self.s, self.a0, self.hist, self.relax = solver, a0, hist, relax

    def solve(self, name, gamma, source, sink, model=None):
        s = self.s
        fvm = s.fvm
        phi = s.state[name]
        bc = s.scalar_bc(name)
        gam_i = fvm.interp(gamma)
        gam_b = gamma[fvm.Pb]
        grad = None if fvm.orthogonal else fvm.grad(phi, bc[0] * phi[fvm.Pb] + bc[1])
        diag, up, lo, rhs = fvm.assemble(s.F_i, s.F_b, gam_i, gam_b, bc, grad_phi=grad,
                                         scheme=s.settings.convection_turb, bounded=s.steady,
                                         phi=phi, nonorth_limit=s.settings.nonorth_limit)
        V = fvm.V
        rhs += np.asarray(source) * V
        diag += np.asarray(sink) * V
        if self.a0:
            diag += self.a0 * V
            rhs -= self.hist[name] * V
        if self.relax < 1.0:
            rhs += (1.0 - self.relax) / self.relax * diag * phi
            diag = diag / self.relax
        A = fvm.matrix(diag, up, lo)
        s.residuals_now[name] = normalized_residual(A, phi, rhs, s.freestream.get(name, 0.0))
        new = solve_linear(A, rhs, phi, s.settings.solver_turb, rtol=0.1 if s.steady else 1e-4)
        floor = s.model.floors.get(name)
        if floor is not None:
            new = np.maximum(new, floor)
        return new


class Solver2D:
    """Écoulement incompressible 2D (ρ = 1, grandeurs cinématiques)."""

    def __init__(self, mesh: Mesh2D, nu: float, boundaries: dict, model: str = "laminar",
                 model_options: dict | None = None, body_force=(0.0, 0.0),
                 initial_U=(0.0, 0.0), turbulence_inflow: dict | None = None,
                 settings: Settings | None = None, reference_velocity: float | None = None):
        self.mesh = mesh
        self.nu = float(nu)
        self.settings = settings or Settings()
        self.body_force = np.asarray(body_force, dtype=float)
        missing = [p.name for p in mesh.patches if p.name not in boundaries]
        if missing:
            raise ValueError(f"Conditions aux limites manquantes pour les patches {missing}.")
        types = {}
        for name, spec in boundaries.items():
            if spec["type"] not in BC_TYPES:
                raise ValueError(f"Type de condition inconnu '{spec['type']}' ({name}). "
                                 f"Choix : {BC_TYPES}")
            types[name] = {"wall": "wall", "symmetry": "symmetry"}.get(spec["type"], "patch")
        mesh.set_patch_types({k: v for k, v in types.items() if k in mesh.patch_types})
        self.boundaries = boundaries
        self.fvm = FVM(mesh)
        fvm = self.fvm
        for attr in ("solver_U", "solver_turb", "solver_p"):
            if getattr(self.settings, attr) == "auto":
                setattr(self.settings, attr, "direct" if fvm.nc < 40000 else "bicgstab")
        if self.settings.relax_U is None:
            # SIMPLEC à 0.9 peut diverger sur maillage très non orthogonal (triangles, hybride)
            nonorth = mesh.quality()["non_orthogonality_max_deg"]
            self.settings.relax_U = 0.9 if nonorth < 30.0 else 0.7
        self._setup_bc()
        uref = reference_velocity
        if uref is None:
            uref = max([np.max(np.linalg.norm(self.U_fixed[self.kindU == 0], axis=1))
                        if np.any(self.kindU == 0) else 0.0, 1e-30])
            uref = max(uref, float(np.linalg.norm(initial_U)), 1e-12)
        self.U_ref = uref
        self.model_name = canonical_name(model)
        self.model = get_model(self.model_name, mesh, self.nu, **(model_options or {}))
        self.model.ops = Ops2D(self)
        ti = {"intensity": 0.001, "viscosity_ratio": 0.1, **(turbulence_inflow or {})}
        self.freestream = self.model.freestream_values(uref, ti["intensity"],
                                                       ti["viscosity_ratio"])
        # champs
        nc = fvm.nc
        self.U = np.tile(np.asarray(initial_U, dtype=float), (nc, 1))
        self.p = np.zeros(nc)
        self.state = {k: np.full(nc, v) for k, v in self.freestream.items()}
        Ub = self.boundary_U(self.U)
        self.F_i = np.sum(fvm.interp(self.U) * fvm.Si, axis=1)
        self.F_b = np.sum(Ub * fvm.Sb, axis=1)
        self.nut = np.zeros(nc)
        self.steady = True
        self.time = 0.0
        self.history: list[dict] = []
        self.residuals_now: dict = {}
        self.has_fixed_p = bool(np.any(self.kindP == 0))

    # ------------------------------------------------------------------ conditions limites
    def _setup_bc(self):
        m, fvm = self.mesh, self.fvm
        nb = fvm.nb
        Cb = m.face_centers[m.n_internal:]
        self.kindU = np.ones(nb, dtype=int)          # 0 fixe, 1 gradient nul, 2 glissement
        self.U_fixed = np.zeros((nb, 2))
        self.kindP = np.ones(nb, dtype=int)          # 0 fixe, 1 gradient nul
        self.p_fixed = np.zeros(nb)
        self.kindT = np.ones(nb, dtype=int)          # 0 amont, 1 gradient nul, 2 paroi
        self.is_wall = np.zeros(nb, dtype=bool)
        self.patch_slices = {}
        for p in m.patches:
            spec = self.boundaries[p.name]
            sl = slice(p.start - m.n_internal, p.start - m.n_internal + p.size)
            self.patch_slices[p.name] = sl
            x, y = Cb[sl, 0], Cb[sl, 1]
            kind = spec["type"]
            if kind == "wall":
                self.kindU[sl] = 0
                uw = spec.get("U", [0.0, 0.0])
                self.U_fixed[sl] = np.column_stack([_eval_expr(uw[0], x, y), _eval_expr(uw[1], x, y)])
                self.kindT[sl] = 2
                self.is_wall[sl] = True
            elif kind == "inlet":
                self.kindU[sl] = 0
                u = spec["U"]
                self.U_fixed[sl] = np.column_stack([_eval_expr(u[0], x, y), _eval_expr(u[1], x, y)])
                self.kindT[sl] = 0
            elif kind == "outlet":
                self.kindP[sl] = 0
                self.p_fixed[sl] = spec.get("p", 0.0)
            elif kind == "symmetry":
                self.kindU[sl] = 2
            elif kind == "farfield":
                uinf = np.asarray(spec["U"], dtype=float)
                inflow = self.fvm.nb_hat[sl] @ uinf < 0.0
                self.kindU[sl] = np.where(inflow, 0, 1)
                self.U_fixed[sl] = uinf
                self.kindP[sl] = np.where(inflow, 1, 0)
                self.p_fixed[sl] = spec.get("p", 0.0)
                self.kindT[sl] = np.where(inflow, 0, 1)

    def vector_bc(self, c: int, U):
        """Coefficients (α, β, γ, δ) de la composante c de la vitesse."""
        fvm = self.fvm
        nb, dp = fvm.nb, fvm.dperp
        a, b, g, d = np.zeros(nb), np.zeros(nb), np.zeros(nb), np.zeros(nb)
        fx = self.kindU == 0
        b[fx] = self.U_fixed[fx, c]
        g[fx] = -1.0 / dp[fx]
        d[fx] = self.U_fixed[fx, c] / dp[fx]
        zg = self.kindU == 1
        a[zg] = 1.0
        sl = self.kindU == 2
        if np.any(sl):
            n = fvm.nb_hat[sl]
            o = 1 - c
            uo = U[fvm.Pb[sl], o]
            a[sl] = 1.0 - n[:, c] ** 2
            b[sl] = -n[:, c] * n[:, o] * uo
            g[sl] = -n[:, c] ** 2 / dp[sl]
            d[sl] = -n[:, c] * n[:, o] * uo / dp[sl]
        return a, b, g, d

    def boundary_U(self, U):
        Pb = self.fvm.Pb
        out = np.empty((self.fvm.nb, 2))
        for c in range(2):
            a, b, _, _ = self.vector_bc(c, U)
            out[:, c] = a * U[Pb, c] + b
        return out

    def pressure_bc(self):
        nb, dp = self.fvm.nb, self.fvm.dperp
        fx = self.kindP == 0
        a = np.where(fx, 0.0, 1.0)
        b = np.where(fx, self.p_fixed, 0.0)
        g = np.where(fx, -1.0 / dp, 0.0)
        d = np.where(fx, self.p_fixed / dp, 0.0)
        return a, b, g, d

    def boundary_p(self, p):
        a, b, _, _ = self.pressure_bc()
        return a * p[self.fvm.Pb] + b

    def scalar_bc(self, name):
        nb, dp = self.fvm.nb, self.fvm.dperp
        a, b, g, d = np.zeros(nb), np.zeros(nb), np.zeros(nb), np.zeros(nb)
        fx = self.kindT == 0
        val = self.freestream.get(name, 0.0)
        b[fx] = val
        g[fx] = -1.0 / dp[fx]
        d[fx] = val / dp[fx]
        a[self.kindT == 1] = 1.0
        w = self.kindT == 2
        if np.any(w):
            wv = self.model.wall_value(name, dp[w])
            b[w] = wv
            g[w] = -1.0 / dp[w]
            d[w] = wv / dp[w]
        return a, b, g, d

    # ------------------------------------------------------------------ outils
    @property
    def u_scale(self) -> float:
        return max(float(np.abs(self.U).max()), self.U_ref)

    def grad_U(self, U):
        Ub = self.boundary_U(U)
        return np.stack([self.fvm.grad(U[:, 0], Ub[:, 0]), self.fvm.grad(U[:, 1], Ub[:, 1])],
                        axis=1)

    def flow(self):
        return FlowField2D(self, self.U, self.grad_U(self.U))

    def update_nut(self):
        self.nut = self.model.eddy_viscosity(self.state, self.flow()) if self.model.variables \
            else np.zeros(self.fvm.nc)

    def _offdiag_mult(self, up, lo, x):
        fvm = self.fvm
        return (np.bincount(fvm.P, up * x[fvm.N], fvm.nc) + np.bincount(fvm.N, lo * x[fvm.P], fvm.nc))

    # ------------------------------------------------------------------ quantité de mouvement
    def _momentum(self, a0=0.0, hist=None, relax=1.0):
        fvm, U = self.fvm, self.U
        nu_eff = self.nu + self.nut
        gam_i = fvm.interp(nu_eff)
        gam_b = np.where(self.is_wall, self.nu, nu_eff[fvm.Pb])
        gradU = self.grad_U(U)
        gUf = fvm.interp(gradU)
        V = fvm.V
        eqs = []
        for c in range(2):
            bc = self.vector_bc(c, U)
            diag, up, lo, rhs = fvm.assemble(self.F_i, self.F_b, gam_i, gam_b, bc,
                                             grad_phi=gradU[:, c, :],
                                             scheme=self.settings.convection_U,
                                             bounded=self.steady, phi=U[:, c],
                                             nonorth_limit=self.settings.nonorth_limit)
            # terme ∇·(ν_eff (∇U)ᵀ) explicite
            tf_i = gam_i * (gUf[:, 0, c] * fvm.Si[:, 0] + gUf[:, 1, c] * fvm.Si[:, 1])
            gb = gradU[fvm.Pb]
            tf_b = gam_b * (gb[:, 0, c] * fvm.Sb[:, 0] + gb[:, 1, c] * fvm.Sb[:, 1])
            rhs += fvm.sum_faces(tf_i, tf_b)
            rhs += self.body_force[c] * V
            if a0:
                diag = diag + a0 * V
                rhs -= hist[:, c] * V
            if relax < 1.0:
                rhs += (1.0 - relax) / relax * diag * U[:, c]
                diag = diag / relax
            eqs.append((diag, up, lo, rhs))
        return eqs

    # ------------------------------------------------------------------ pression
    def _pressure_correction(self, eqs, ddt_corr=None, relax_p=1.0):
        fvm, s = self.fvm, self.settings
        V = fvm.V
        H = np.column_stack([rhs - self._offdiag_mult(up, lo, self.U[:, c])
                             for c, (diag, up, lo, rhs) in enumerate(eqs)])
        diag_u, diag_v = eqs[0][0], eqs[1][0]
        aP = 0.5 * (diag_u + diag_v)
        rAU = V / aP
        HbyA = np.column_stack([H[:, 0] / diag_u, H[:, 1] / diag_v])
        pb = self.boundary_p(self.p)
        gradp = fvm.grad(self.p, pb)
        HbyA_b = self.boundary_U(HbyA)
        phiHbyA_i = np.sum(fvm.interp(HbyA) * fvm.Si, axis=1)
        phiHbyA_b = np.sum(HbyA_b * fvm.Sb, axis=1)
        rAtU = rAU
        if s.algorithm.upper() == "SIMPLEC" and self.steady:
            H1 = -0.5 * (np.bincount(fvm.P, eqs[0][1], fvm.nc) + np.bincount(fvm.N, eqs[0][2], fvm.nc)
                         + np.bincount(fvm.P, eqs[1][1], fvm.nc) + np.bincount(fvm.N, eqs[1][2], fvm.nc))
            rAtU = V / np.maximum(aP - H1, 1e-3 * aP)
            # (rAtU − rAU) ∂p/∂n |S| avec le MÊME gradient normal (partie non orthogonale
            # comprise) que l'équation de pression, faces internes et frontières.
            drA = rAtU - rAU
            gam_d = fvm.interp(drA)
            sn = gam_d * fvm.g * (self.p[fvm.N] - self.p[fvm.P])
            if not fvm.orthogonal:
                sn = sn + fvm.nonorth_flux(gam_d, gradp, self.p, s.nonorth_limit)
            phiHbyA_i = phiHbyA_i + sn
            pa, pbb, pg, pd = self.pressure_bc()
            phiHbyA_b = phiHbyA_b + drA[fvm.Pb] * fvm.magSb * (pg * self.p[fvm.Pb] + pd)
            HbyA = HbyA - (rAU - rAtU)[:, None] * gradp
        if ddt_corr is not None:
            phiHbyA_i = phiHbyA_i + fvm.interp(rAU) * ddt_corr
        gam_i = fvm.interp(rAtU)
        gam_b = rAtU[fvm.Pb]
        bc = self.pressure_bc()
        a_, b_, g_, d_ = bc
        div_hbya = fvm.sum_faces(phiHbyA_i, phiHbyA_b)
        p_new = self.p
        n_loops = 1 if fvm.orthogonal else 1 + max(int(s.n_nonorth), 0)
        nf_used = 0.0
        for k in range(n_loops):
            # correction non orthogonale explicite, recalculée avec la dernière pression
            gp_k = gradp if k == 0 else fvm.grad(p_new, a_ * p_new[fvm.Pb] + b_)
            if not fvm.orthogonal:
                nf_used = fvm.nonorth_flux(gam_i, gp_k, p_new, s.nonorth_limit)
            diag, up, lo, rhs = fvm.assemble(np.zeros(fvm.ni), np.zeros(fvm.nb), gam_i, gam_b,
                                             bc, grad_phi=gp_k, phi=p_new,
                                             nonorth_limit=s.nonorth_limit)
            rhs -= div_hbya
            if not self.has_fixed_p:
                diag[0] += diag[0]
            A = fvm.matrix(diag, up, lo)
            if k == 0:
                self.residuals_now["p"] = normalized_residual(A, self.p, rhs, self.u_scale ** 2)
            p_new = solve_linear(A, rhs, p_new, s.solver_p, rtol=1e-6)
        # flux conservatifs (cohérents avec la dernière équation résolue)
        self.F_i = phiHbyA_i - gam_i * fvm.g * (p_new[fvm.N] - p_new[fvm.P]) - nf_used
        self.F_b = phiHbyA_b - gam_b * fvm.magSb * (g_ * p_new[fvm.Pb] + d_)
        self.p = self.p + relax_p * (p_new - self.p)
        gradp = fvm.grad(self.p, self.boundary_p(self.p))
        self.U = HbyA - rAtU[:, None] * gradp

    # ------------------------------------------------------------------ turbulence
    def _turbulence(self, a0=0.0, hist=None, relax=1.0):
        if not self.model.variables:
            return
        step = Step2D(self, a0, hist or {}, relax)
        new = self.model.update(self.state, self.flow(), step)
        self.state.update(new)
        self.update_nut()

    # ------------------------------------------------------------------ stationnaire
    def run_steady(self, max_iter=None, tol=None, verbose=False, log_every=50, callback=None):
        s = self.settings
        max_iter = max_iter or s.max_iter
        tol = tol or s.tol
        self.steady = True
        relax_p = s.relax_p if s.algorithm.upper() == "SIMPLEC" else min(s.relax_p, 0.3)
        self.update_nut()
        t0 = time.perf_counter()
        converged = False
        for it in range(1, max_iter + 1):
            self.residuals_now = {}
            eqs = self._momentum(relax=s.relax_U)
            for c, (diag, up, lo, rhs) in enumerate(eqs):
                A = self.fvm.matrix(diag, up, lo)
                gp = self.fvm.grad(self.p, self.boundary_p(self.p))[:, c] * self.fvm.V
                b = rhs - gp
                self.residuals_now["Ux" if c == 0 else "Uy"] = normalized_residual(
                    A, self.U[:, c], b, self.u_scale)
                self.U[:, c] = solve_linear(A, b, self.U[:, c], s.solver_U, rtol=0.1)
            self._pressure_correction(eqs, relax_p=relax_p)
            self._turbulence(relax=s.relax_turb)
            cont = float(np.sum(np.abs(self.fvm.div(self.F_i, self.F_b))) /
                         (self.u_scale * np.sum(self.fvm.mesh.magSf) + 1e-300))
            rec = {"iteration": it, **self.residuals_now, "continuity": cont}
            self.history.append(rec)
            if not all(np.isfinite(v) for v in rec.values()) or not np.all(np.isfinite(self.U)):
                raise FloatingPointError(f"Divergence à l'itération {it}.")
            if verbose and (it % log_every == 0 or it == 1):
                print("  it %5d  " % it + "  ".join(f"{k}={v:.2e}" for k, v in rec.items()
                                                    if k != "iteration"))
            if callback:
                callback(self, it)
            if max(v for k, v in self.residuals_now.items()) < tol:
                converged = True
                break
        self.converged = converged
        self.iterations = it
        self.wall_time = time.perf_counter() - t0
        if verbose:
            print(f"  {'convergé' if converged else 'NON convergé'} en {it} itérations "
                  f"({self.wall_time:.1f} s)")
        return converged

    # ------------------------------------------------------------------ instationnaire
    def run_transient(self, dt: float, t_end: float, verbose=False, log_every=50,
                      callback=None, probes=None):
        """PIMPLE. probes : fonction(solver) -> dict appelée à chaque pas (historique)."""
        s, fvm = self.settings, self.fvm
        self.steady = False
        self.update_nut()
        n_steps = int(round((t_end - self.time) / dt))
        U_old, U_old2 = self.U.copy(), None
        F_old, F_old2 = self.F_i.copy(), None
        st_old, st_old2 = {k: v.copy() for k, v in self.state.items()}, None
        series = []
        t0 = time.perf_counter()
        for n in range(1, n_steps + 1):
            self.time += dt
            if s.time_scheme == "euler" or U_old2 is None:
                a0 = 1.0 / dt
                histU = -U_old / dt
                histT = {k: -v / dt for k, v in st_old.items()}
                corr_old = (F_old - np.sum(fvm.interp(U_old) * fvm.Si, axis=1))
                coef = 1.0 - np.minimum(np.abs(corr_old) / (np.abs(F_old) + 1e-30), 1.0)
                ddt_corr = coef * corr_old / dt
            else:
                a0 = 1.5 / dt
                histU = (-2.0 * U_old + 0.5 * U_old2) / dt
                histT = {k: (-2.0 * st_old[k] + 0.5 * st_old2[k]) / dt for k in st_old}
                c1 = F_old - np.sum(fvm.interp(U_old) * fvm.Si, axis=1)
                c2 = F_old2 - np.sum(fvm.interp(U_old2) * fvm.Si, axis=1)
                coef = 1.0 - np.minimum(np.abs(c1) / (np.abs(F_old) + 1e-30), 1.0)
                ddt_corr = coef * (2.0 * c1 - 0.5 * c2) / dt
            for outer in range(s.n_outer):
                self.residuals_now = {}
                final = outer == s.n_outer - 1
                eqs = self._momentum(a0=a0, hist=histU, relax=1.0)
                for c, (diag, up, lo, rhs) in enumerate(eqs):
                    A = fvm.matrix(diag, up, lo)
                    gp = fvm.grad(self.p, self.boundary_p(self.p))[:, c] * fvm.V
                    b = rhs - gp
                    self.residuals_now["Ux" if c == 0 else "Uy"] = normalized_residual(
                        A, self.U[:, c], b, self.u_scale)
                    self.U[:, c] = solve_linear(A, b, self.U[:, c], s.solver_U, rtol=1e-4)
                for _ in range(s.n_corr):
                    self._pressure_correction(eqs, ddt_corr=ddt_corr)
                if final or s.turbulence_every_outer:
                    self._turbulence(a0=a0, hist=histT)
            if not np.all(np.isfinite(self.U)):
                raise FloatingPointError(f"Divergence à t = {self.time:.4g}.")
            U_old2, U_old = U_old, self.U.copy()
            F_old2, F_old = F_old, self.F_i.copy()
            st_old2, st_old = st_old, {k: v.copy() for k, v in self.state.items()}
            rec = {"time": self.time, **self.residuals_now}
            if probes:
                rec.update(probes(self))
            series.append(rec)
            if verbose and (n % log_every == 0 or n == 1):
                print(f"  t={self.time:.4f}  " + "  ".join(
                    f"{k}={v:.3e}" for k, v in rec.items() if k != "time"))
            if callback:
                callback(self, n)
        self.wall_time = time.perf_counter() - t0
        return series

    # ------------------------------------------------------------------ post-traitement
    def forces(self, patches=None):
        """Efforts (par unité d'envergure, ρ = 1) sur des patches : pression + frottement."""
        fvm = self.fvm
        names = patches or [p.name for p in self.mesh.patches if p.type == "wall"]
        out = {}
        pb = self.boundary_p(self.p)
        for name in names:
            sl = self.patch_slices[name]
            Pb = fvm.Pb[sl]
            Fp = np.sum(pb[sl, None] * fvm.Sb[sl], axis=0)
            du = self.U[Pb] - self.U_fixed[sl]
            Fv = np.sum((self.nu * fvm.magSb[sl] / fvm.dperp[sl])[:, None] * du, axis=0)
            out[name] = {"pressure": Fp, "viscous": Fv, "total": Fp + Fv}
        return out

    def wall_shear(self, patch):
        """(abscisse curviligne implicite) centres de faces, τ_w signé (tangente locale), y⁺."""
        fvm = self.fvm
        sl = self.patch_slices[patch]
        Pb = fvm.Pb[sl]
        n = fvm.nb_hat[sl]
        t = np.column_stack([-n[:, 1], n[:, 0]])
        du = self.U[Pb] - self.U_fixed[sl]
        ut = np.sum(du * t, axis=1)
        tau = self.nu * ut / fvm.dperp[sl]
        yplus = fvm.dperp[sl] * np.sqrt(np.abs(tau)) / self.nu
        xf = self.mesh.face_centers[self.mesh.n_internal:][sl]
        return xf, tau, yplus

    def fields(self) -> dict:
        out = {"U": self.U, "p": self.p, "U_mag": np.linalg.norm(self.U, axis=1),
               "vorticity": self.flow().vorticity}
        if self.model.variables:
            out["nut_over_nu"] = self.nut / self.nu
            out.update(self.state)
        out["wall_distance"] = self.mesh.wall_distance
        return out
