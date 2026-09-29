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
  axis      axe de révolution (axisymétrique) : comme symmetry, surface nulle
Les patches périodiques sont gérés par le maillage (pas de condition à donner).

Axisymétrique (axisymmetric=True ; x = axe, y = rayon, sans rotation propre) : géométrie
pondérée par r (voir fvm.py) ; terme circonférentiel de la contrainte visqueuse
−τ_θθ/r = −2 ν_eff u_r / r² (implicite) ; taux de déformation complété par (u_r / r)².
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from ..mesh2d.mesh import Mesh2D
from ..models import canonical_name, get_model
from .fvm import FVM, normalized_residual

BC_TYPES = ("wall", "inlet", "outlet", "symmetry", "farfield", "axis")

# Schémas en temps. Limites de stabilité des schémas explicites (convection linearUpwind) :
# co_max = Courant convectif (mesuré sur le tourbillon de Taylor-Green advecté, avec marge :
# instabilité observée à 1.04 / 1.31 / 1.47 / 1.66 / 0.73 pour rk1..rk4 / ab2 ; rk1 est de
# plus faiblement instable en convection pure), dn_max = nombre de diffusion (théorie :
# intervalle réel de stabilité / 2 : 2, 2, 2.51, 2.79, 1).
TIME_SCHEMES = {
    "euler": dict(kind="implicit", order=1, label="Euler implicite (PIMPLE)"),
    "backward": dict(kind="implicit", order=2, label="BDF2 / backward (PIMPLE)"),
    "crankNicolson": dict(kind="implicit", order=2, label="Crank-Nicolson (PIMPLE)"),
    "rk1": dict(kind="explicit", order=1, co_max=0.5, dn_max=1.0, stages=1,
                label="Euler explicite + projection"),
    "rk2": dict(kind="explicit", order=2, co_max=1.0, dn_max=1.0, stages=2,
                label="RK2 (Heun, SSP) + projection"),
    "rk3": dict(kind="explicit", order=3, co_max=1.2, dn_max=1.25, stages=3,
                label="RK3 SSP (Shu-Osher) + projection"),
    "rk4": dict(kind="explicit", order=4, co_max=1.4, dn_max=1.39, stages=4,
                label="RK4 classique + projection"),
    "ab2": dict(kind="explicit", order=2, co_max=0.5, dn_max=0.5, stages=1,
                label="Adams-Bashforth 2 + projection"),
}
_ALIASES = {"bdf1": "euler", "implicit_euler": "euler", "bdf2": "backward",
            "cn": "crankNicolson", "cranknicolson": "crankNicolson", "crank-nicolson":
            "crankNicolson", "explicit_euler": "rk1", "heun": "rk2", "ssprk3": "rk3"}

# tableaux de Butcher (A strictement triangulaire inférieure, b, c)
RK_TABLES = {
    "rk1": ([[0.0]], [1.0], [0.0]),
    "rk2": ([[0.0, 0.0], [1.0, 0.0]], [0.5, 0.5], [0.0, 1.0]),
    "rk3": ([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.25, 0.25, 0.0]], [1 / 6, 1 / 6, 2 / 3],
            [0.0, 1.0, 0.5]),
    "rk4": ([[0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0], [0.0, 0.5, 0.0, 0.0],
             [0.0, 0.0, 1.0, 0.0]], [1 / 6, 1 / 3, 1 / 3, 1 / 6], [0.0, 0.5, 0.5, 1.0]),
}


def canonical_time_scheme(name: str) -> str:
    if str(name).lower() == "auto":
        return "auto"
    key = _ALIASES.get(str(name).lower(), str(name))
    if key not in TIME_SCHEMES:
        low = {k.lower(): k for k in TIME_SCHEMES}
        key = low.get(key.lower(), key)
    if key not in TIME_SCHEMES:
        raise ValueError(f"Schéma en temps inconnu '{name}'. Choix : {list(TIME_SCHEMES)}")
    return key


@dataclass
class Settings:
    """Paramètres numériques (valeurs par défaut raisonnables pour débuter)."""
    algorithm: str = "SIMPLEC"            # SIMPLE | SIMPLEC (stationnaire)
    relax_U: float | None = None          # None = auto : 0.9 (maillage ~orthogonal) / 0.7
    relax_p: float = 1.0                  # SIMPLE : ~0.3 ; SIMPLEC : 1
    relax_turb: float = 0.8
    # stationnaire pseudo-transitoire : pas de pseudo-temps global (remplace relax_U/relax_T ;
    # conseillé sur maillages très fins et étirés, voir README)
    pseudo_dt: float | None = None
    pseudo_cfl: float | None = None
    convection_U: str = "linearUpwind"    # upwind | linearUpwind
    convection_turb: str = "upwind"
    convection_T: str = "linearUpwind"
    relax_T: float = 0.9
    max_iter: int = 3000
    tol: float = 1e-5                     # résidus normalisés (OpenFOAM) de tous les champs
    # critère complémentaire (comme les moniteurs de Fluent) : arrêt quand les efforts sur
    # les parois varient de moins de monitor_tol (relatif) sur monitor_window itérations
    monitor_tol: float | None = None
    monitor_window: int = 100
    # solveurs linéaires (voir microrans.linalg) : auto | direct | amg | bicgstab | pyamg
    solver_p: str = "auto"                # auto : LU si < 20 000 cellules, sinon AMG + CG
    solver_U: str = "auto"                # auto : BiCGStab + Jacobi
    solver_turb: str = "auto"
    # instationnaire (voir TIME_SCHEMES)
    time_scheme: str = "auto"             # auto | euler | backward | crankNicolson | rk1..4 | ab2
    adjust_dt: bool = False               # pas de temps adaptatif (comme adjustTimeStep)
    max_co: float = 1.0                   # Courant visé si adjust_dt
    max_dt: float = float("inf")
    cn_theta: float = 0.5                 # Crank-Nicolson : 0.5 (ordre 2), >0.5 plus dissipatif
    ddt_phi_coeff: float | None = 1.0     # correction ddtCorr (None = coefficient OpenFOAM)
    # traitement pariétal : resolved (y⁺ ≲ 1, défaut) | wall_function (loi de Spalding,
    # y⁺ ≈ 1 à 300 ; SA, k-ω, SST — pas le k-ε Launder-Sharma bas-Reynolds)
    wall_treatment: str = "resolved"
    # matériel : cpu (NumPy/SciPy) | gpu (CuPy, expérimental : voir microrans/backend.py)
    backend: str = "cpu"
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


def _active(a0) -> bool:
    """Terme temporel (ou pseudo-temporel) présent : scalaire non nul ou tableau par cellule."""
    return a0 is not None and (not isinstance(a0, (int, float)) or a0 != 0.0)


class FlowField2D:
    """Grandeurs du champ moyen utilisées par les modèles de turbulence."""

    def __init__(self, solver: "Solver2D", U, gradU):
        self.solver = solver
        self.U = U
        self.gradU = gradU               # gradU[:, i, j] = ∂u_i/∂x_j
        self._d2 = None

    def _wall_override(self, val):
        """Lois de paroi : dans les cellules pariétales, le gradient de Green-Gauss (vitesse
        nulle à la paroi, sous-couche non résolue) surestime fortement le cisaillement ; on
        le remplace par celui de la loi de paroi, u_τ²/(ν + κ u_τ y) (comme la production
        modifiée des cellules pariétales d'OpenFOAM)."""
        s = self.solver
        if not s.wall_function:
            return val
        mask, sw = s._wall_cell_shear()
        return s.xp.where(mask, sw, val)

    @property
    def strain(self):
        xp = self.solver.xp
        g = self.gradU
        s2 = 2 * g[:, 0, 0] ** 2 + 2 * g[:, 1, 1] ** 2 + (g[:, 0, 1] + g[:, 1, 0]) ** 2
        fvm = self.solver.fvm
        if fvm.axisymmetric:
            s2 = s2 + 2 * (self.U[:, 1] / fvm.radius) ** 2      # S_θθ = u_r / r
        return self._wall_override(xp.sqrt(s2))

    @property
    def vorticity(self):
        xp = self.solver.xp
        g = self.gradU
        return self._wall_override(xp.abs(g[:, 1, 0] - g[:, 0, 1]))

    @property
    def second_derivative_sq(self):
        xp = self.solver.xp
        if self._d2 is None:
            s = self.solver
            tot = xp.zeros(s.fvm.nc)
            for i in range(2):
                for j in range(2):
                    comp = self.gradU[:, i, j]
                    gg = s.fvm.grad(comp, comp[s.fvm.Pb])
                    tot += xp.sum(gg ** 2, axis=1)
            self._d2 = tot
        return self._d2


class Ops2D:
    """Gradients pour les modèles (valeur pariétale imposée, gradient nul ailleurs)."""

    def __init__(self, solver: "Solver2D"):
        self.s = solver

    def _grad(self, f, name):
        xp = self.s.xp
        s = self.s
        fb = f[s.fvm.Pb].copy()
        if name is not None and xp.any(s.is_wall) and not (
                s.wall_function and name in ("k", "sqrt_k")):
            base = "k" if name == "sqrt_k" else name
            wv = s.model.wall_value(base, s.fvm.dperp[s.is_wall])
            fb[s.is_wall] = xp.sqrt(xp.maximum(wv, 0.0)) if name == "sqrt_k" else wv
        return s.fvm.grad(f, fb)

    def grad_sq(self, f, name=None):
        xp = self.s.xp
        return xp.sum(self._grad(f, name) ** 2, axis=1)

    def grad_dot(self, f, g, fname=None, gname=None):
        xp = self.s.xp
        return xp.sum(self._grad(f, fname) * self._grad(g, gname), axis=1)


class Step2D:
    """Résolution implicite d'une équation de transport de turbulence en 2D."""

    def __init__(self, solver: "Solver2D", a0: float, hist: dict, relax: float):
        self.s, self.a0, self.hist, self.relax = solver, a0, hist, relax

    def solve(self, name, gamma, source, sink, model=None):
        xp = self.s.xp
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
        rhs += xp.asarray(source) * V
        diag += xp.asarray(sink) * V
        if _active(self.a0):
            diag += self.a0 * V
            rhs -= self.hist[name] * V
        if self.relax < 1.0:
            rhs += (1.0 - self.relax) / self.relax * diag * phi
            diag = diag / self.relax
        if s.wall_function and name == "omega":
            # valeur imposée dans les cellules pariétales (ligne de matrice remplacée)
            mask, val = s._omega_wall_cells()
            up = xp.where(mask[fvm.P], 0.0, up)
            lo = xp.where(mask[fvm.N], 0.0, lo)
            rhs = xp.where(mask, diag * val, rhs)
        A = fvm.matrix(diag, up, lo)
        s.residuals_now[name] = normalized_residual(A, phi, rhs, s.freestream.get(name, 0.0))
        new = s.fvm.lin.solve(A, rhs, phi, s.settings.solver_turb,
                              rtol=0.1 if s.steady else 1e-4, tag=name)
        floor = s.model.floors.get(name)
        if floor is not None:
            new = xp.maximum(new, floor)
        return new


class Solver2D:
    """Écoulement incompressible 2D (ρ = 1, grandeurs cinématiques)."""

    def __init__(self, mesh: Mesh2D, nu: float, boundaries: dict, model: str = "laminar",
                 model_options: dict | None = None, body_force=(0.0, 0.0),
                 initial_U=(0.0, 0.0), turbulence_inflow: dict | None = None,
                 settings: Settings | None = None, reference_velocity: float | None = None,
                 energy: dict | None = None, axisymmetric: bool = False):
        self.mesh = mesh
        self.axisymmetric = bool(axisymmetric)
        # thermique (optionnelle) : Pr, Pr_t, beta (Boussinesq), T_ref, gravity, T0, source
        self.energy = None if energy is None else {
            "Pr": 0.71, "Pr_t": 0.85, "beta": 0.0, "T_ref": 0.0, "gravity": (0.0, -9.81),
            **energy}
        self.nu = float(nu)
        self.settings = settings or Settings()
        # force volumique : constante (fx, fy) ou fonction du temps t -> (fx, fy)
        self.body_force = body_force if callable(body_force) else np.asarray(body_force, float)
        self._t_eval = 0.0
        missing = [p.name for p in mesh.patches if p.name not in boundaries]
        if missing:
            raise ValueError(f"Conditions aux limites manquantes pour les patches {missing}.")
        types = {}
        for name, spec in boundaries.items():
            if spec["type"] not in BC_TYPES:
                raise ValueError(f"Type de condition inconnu '{spec['type']}' ({name}). "
                                 f"Choix : {BC_TYPES}")
            types[name] = {"wall": "wall", "symmetry": "symmetry",
                           "axis": "symmetry"}.get(spec["type"], "patch")
        mesh.set_patch_types({k: v for k, v in types.items() if k in mesh.patch_types})
        self.boundaries = boundaries
        if self.axisymmetric:
            radial = [abs(float(np.asarray(body_force, float)[1]))
                      if not callable(body_force) else 0.0]
            if self.energy is not None and self.energy.get("beta"):
                radial.append(abs(float(self.energy["gravity"][1])))
            if max(radial) > 0.0:
                raise ValueError("Axisymétrique : la force volumique et la gravité doivent être "
                                 "portées par l'axe x (composante radiale y nulle).")
        self.fvm = FVM(mesh, self.settings.backend, self.axisymmetric)
        fvm = self.fvm
        self.backend, self.xp = fvm.backend, fvm.xp
        xp = self.xp
        if self.settings.relax_U is None:
            # SIMPLEC à 0.9 peut diverger sur maillage très non orthogonal (triangles, hybride)
            nonorth = mesh.quality()["non_orthogonality_max_deg"]
            self.settings.relax_U = 0.9 if nonorth < 30.0 else 0.7
        kindU_h, U_fixed_h, kindP_h = self._setup_bc()
        uref = reference_velocity
        if uref is None:
            uref = max([float(np.max(np.linalg.norm(U_fixed_h[kindU_h == 0], axis=1)))
                        if np.any(kindU_h == 0) else 0.0, 1e-30])
            uref = max(uref, float(np.linalg.norm(initial_U)), 1e-12)
        self.U_ref = uref
        self.model_name = canonical_name(model)
        self.model = get_model(self.model_name, mesh, self.nu, **(model_options or {}))
        self.model.ops = Ops2D(self)
        self.model.d = self.backend.asarray(self.model.d)
        self.model.wall_nodes = self.backend.asarray(self.model.wall_nodes)
        self.wall_function = self.settings.wall_treatment == "wall_function"
        if self.settings.wall_treatment not in ("resolved", "wall_function"):
            raise ValueError("wall_treatment : 'resolved' ou 'wall_function'.")
        if self.wall_function and self.model_name == "ke":
            raise ValueError("Lois de paroi incompatibles avec le k-ε Launder-Sharma "
                             "(modèle bas-Reynolds, y⁺ ≈ 1 obligatoire) : utiliser SST, k-ω "
                             "ou SA, ou wall_treatment = 'resolved'.")
        self.nu_wall = xp.full(fvm.nb, self.nu)     # viscosité effective aux faces de paroi
        self.u_tau_wall = xp.zeros(fvm.nb)
        ti = {"intensity": 0.001, "viscosity_ratio": 0.1, **(turbulence_inflow or {})}
        self.freestream = self.model.freestream_values(uref, ti["intensity"],
                                                       ti["viscosity_ratio"])
        # champs
        nc = fvm.nc
        self.U = xp.tile(xp.asarray(np.asarray(initial_U, dtype=float)), (nc, 1))
        self.p = xp.zeros(nc)
        self.state = {k: xp.full(nc, float(v)) for k, v in self.freestream.items()}
        Ub = self.boundary_U(self.U)
        self.F_i = xp.sum(fvm.interp(self.U) * fvm.Si, axis=1)
        self.F_b = xp.sum(Ub * fvm.Sb, axis=1)
        self.nut = xp.zeros(nc)
        if self.energy is not None:
            self.T = xp.full(nc, float(self.energy.get("T0", self.energy["T_ref"])))
        self.steady = True
        self.time = 0.0
        self.dt = 0.0
        self.history: list[dict] = []
        self.residuals_now: dict = {}
        # reprise (fv2d/restart.py) : compteur global d'itérations, niveaux de temps
        # précédents, historique d'efforts du calcul repris
        self.iterations_total = 0
        self._restart_hist = None
        self.series_restart: list[dict] = []
        self.series: list[dict] = []
        self.averager = None               # moyennes temporelles (sampling.TimeAverage)
        self.has_fixed_p = bool(np.any(kindP_h == 0))

    # ------------------------------------------------------------------ conditions limites
    def _setup_bc(self):
        """Types et valeurs de conditions aux limites par face (préparés sur CPU, puis
        transférés sur le matériel de calcul). Retourne les copies CPU utiles."""
        m, fvm = self.mesh, self.fvm
        nb = fvm.nb
        Cb = m.face_centers[m.n_internal:]
        nhat = m.nf[m.n_internal:]
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
            elif kind in ("symmetry", "axis"):
                self.kindU[sl] = 2
            elif kind == "farfield":
                uinf = np.asarray(spec["U"], dtype=float)
                inflow = nhat[sl] @ uinf < 0.0
                self.kindU[sl] = np.where(inflow, 0, 1)
                self.U_fixed[sl] = uinf
                self.kindP[sl] = np.where(inflow, 1, 0)
                self.p_fixed[sl] = spec.get("p", 0.0)
                self.kindT[sl] = np.where(inflow, 0, 1)
        # température : 0 valeur imposée, 1 gradient (flux) imposé
        self.kindTemp = np.ones(nb, dtype=int)
        self.T_fixed = np.zeros(nb)
        self.T_grad = np.zeros(nb)
        if self.energy is not None:
            Tref = float(self.energy["T_ref"])
            alpha = self.nu / float(self.energy["Pr"])
            for p in m.patches:
                spec, sl = self.boundaries[p.name], self.patch_slices[p.name]
                kind = spec["type"]
                if kind == "wall" and "T" in spec:
                    self.kindTemp[sl] = 0
                    self.T_fixed[sl] = spec["T"]
                elif kind == "wall" and "q" in spec:
                    # flux q entrant dans le fluide (unités cinématiques q/(ρ c_p)) :
                    # q = −(−α∇T)·n_sortante = α ∂T/∂n  ⇒  ∂T/∂n = q/α
                    self.T_grad[sl] = float(spec["q"]) / alpha
                elif kind == "inlet":
                    self.kindTemp[sl] = 0
                    self.T_fixed[sl] = spec.get("T", Tref)
                elif kind == "farfield":
                    inflow = nhat[sl] @ np.asarray(spec["U"], dtype=float) < 0.0
                    self.kindTemp[sl] = np.where(inflow, 0, 1)
                    self.T_fixed[sl] = spec.get("T", Tref)
        host = (self.kindU.copy(), self.U_fixed.copy(), self.kindP.copy())
        A = self.backend.asarray
        self.kindTemp, self.T_fixed = A(self.kindTemp), A(self.T_fixed)
        self.T_grad = A(self.T_grad)
        self.kindU, self.U_fixed, self.kindP = A(self.kindU), A(self.U_fixed), A(self.kindP)
        self.p_fixed, self.kindT, self.is_wall = A(self.p_fixed), A(self.kindT), A(self.is_wall)
        return host

    def vector_bc(self, c: int, U):
        """Coefficients (α, β, γ, δ) de la composante c de la vitesse."""
        xp = self.xp
        fvm = self.fvm
        nb, dp = fvm.nb, fvm.dperp
        a, b, g, d = xp.zeros(nb), xp.zeros(nb), xp.zeros(nb), xp.zeros(nb)
        fx = self.kindU == 0
        b[fx] = self.U_fixed[fx, c]
        g[fx] = -1.0 / dp[fx]
        d[fx] = self.U_fixed[fx, c] / dp[fx]
        zg = self.kindU == 1
        a[zg] = 1.0
        sl = self.kindU == 2
        if xp.any(sl):
            n = fvm.nb_hat[sl]
            o = 1 - c
            uo = U[fvm.Pb[sl], o]
            a[sl] = 1.0 - n[:, c] ** 2
            b[sl] = -n[:, c] * n[:, o] * uo
            g[sl] = -n[:, c] ** 2 / dp[sl]
            d[sl] = -n[:, c] * n[:, o] * uo / dp[sl]
        return a, b, g, d

    def boundary_U(self, U):
        xp = self.xp
        Pb = self.fvm.Pb
        out = xp.empty((self.fvm.nb, 2))
        for c in range(2):
            a, b, _, _ = self.vector_bc(c, U)
            out[:, c] = a * U[Pb, c] + b
        return out

    def pressure_bc(self):
        """Pression imposée (sorties) ou gradient normal imposé ailleurs. Avec une force
        volumique, ce gradient vaut f·n (« fixedFluxPressure » d'OpenFOAM) : l'équilibre
        hydrostatique est respecté jusqu'à la paroi, sinon des courants parasites
        apparaissent dans les cellules pariétales."""
        xp = self.xp
        dp = self.fvm.dperp
        fx = self.kindP == 0
        G = self._boundary_force_normal() if self._has_force() else 0.0
        a = xp.where(fx, 0.0, 1.0)
        b = xp.where(fx, self.p_fixed, dp * G)
        g = xp.where(fx, -1.0 / dp, 0.0)
        d = xp.where(fx, self.p_fixed / dp, G)
        return a, b, g, d

    def _boundary_force(self):
        """Force volumique aux faces frontières (flottabilité avec la température de paroi)."""
        xp = self.xp
        f = xp.zeros((self.fvm.nb, 2)) + self.body_force_at(self._t_eval)[None, :]
        e = self.energy
        if e is not None and e["beta"]:
            g = xp.asarray(np.asarray(e["gravity"], dtype=float))
            Tb = self.boundary_T(self.T)
            f = f - float(e["beta"]) * (Tb - float(e["T_ref"]))[:, None] * g[None, :]
        return f

    def _boundary_force_normal(self):
        return self.xp.sum(self._boundary_force() * self.fvm.nb_hat, axis=1)

    def boundary_p(self, p):
        a, b, _, _ = self.pressure_bc()
        return a * p[self.fvm.Pb] + b

    def scalar_bc(self, name):
        xp = self.xp
        nb, dp = self.fvm.nb, self.fvm.dperp
        a, b, g, d = xp.zeros(nb), xp.zeros(nb), xp.zeros(nb), xp.zeros(nb)
        fx = self.kindT == 0
        val = self.freestream.get(name, 0.0)
        b[fx] = val
        g[fx] = -1.0 / dp[fx]
        d[fx] = val / dp[fx]
        a[self.kindT == 1] = 1.0
        w = self.kindT == 2
        if self.wall_function and name == "k":
            a[w] = 1.0                        # kqRWallFunction : gradient nul
            w = xp.zeros_like(w)
        if xp.any(w):
            wv = self.model.wall_value(name, dp[w])
            b[w] = wv
            g[w] = -1.0 / dp[w]
            d[w] = wv / dp[w]
        return a, b, g, d

    # ------------------------------------------------------------------ outils
    @property
    def u_scale(self) -> float:
        xp = self.xp
        return max(float(xp.abs(self.U).max()), self.U_ref)

    def grad_U(self, U):
        xp = self.xp
        Ub = self.boundary_U(U)
        return xp.stack([self.fvm.grad(U[:, 0], Ub[:, 0]), self.fvm.grad(U[:, 1], Ub[:, 1])],
                        axis=1)

    def flow(self):
        return FlowField2D(self, self.U, self.grad_U(self.U))

    def update_nut(self):
        xp = self.xp
        self.nut = self.model.eddy_viscosity(self.state, self.flow()) if self.model.variables \
            else xp.zeros(self.fvm.nc)

    def _offdiag_mult(self, up, lo, x):
        xp = self.xp
        fvm = self.fvm
        return (xp.bincount(fvm.P, up * x[fvm.N], fvm.nc) + xp.bincount(fvm.N, lo * x[fvm.P], fvm.nc))

    # ------------------------------------------------------------------ quantité de mouvement
    def _has_force(self):
        e = self.energy
        if e is not None and e["beta"]:
            return True
        return callable(self.body_force) or bool(np.any(np.asarray(self.body_force) != 0))

    def body_force_at(self, t):
        xp = self.xp
        bf = self.body_force
        return xp.asarray(bf(t) if callable(bf) else bf, dtype=float)

    # ------------------------------------------------------------------ lois de paroi
    KAPPA, B_LOG = 0.41, 5.2

    def _spalding_u_tau(self, Ut, y):
        """u_τ tel que (y⁺, u⁺) vérifie la loi de Spalding (1961), valable de la sous-couche
        visqueuse à la zone logarithmique (Newton vectorisé)."""
        xp = self.xp
        k, E = self.KAPPA, np.exp(-self.KAPPA * self.B_LOG)
        nu = self.nu
        ut = xp.maximum(xp.sqrt(nu * Ut / y), 1e-12 * self.U_ref)       # sous-couche
        for _ in range(30):
            up = Ut / ut
            ku = xp.minimum(k * up, 50.0)
            ex = xp.exp(ku)
            g = up + E * (ex - 1 - ku - ku ** 2 / 2 - ku ** 3 / 6)
            dg = 1 + E * k * (ex - 1 - ku - ku ** 2 / 2)
            f = y * ut / nu - g
            df = y / nu + dg * Ut / ut ** 2
            ut = xp.maximum(ut - f / df, 0.5 * ut)
        return ut

    def _update_wall_function(self):
        if not self.wall_function:
            return
        xp, fvm = self.xp, self.fvm
        w = self.is_wall
        n = fvm.nb_hat
        du = self.U[fvm.Pb] - self.U_fixed
        ut = xp.abs(du[:, 0] * n[:, 1] - du[:, 1] * n[:, 0])            # |U_t| relative
        y = fvm.dperp
        u_tau = self._spalding_u_tau(xp.maximum(ut, 1e-30), y)
        nu_w = xp.where(ut > 1e-12 * self.U_ref, u_tau ** 2 * y / xp.maximum(ut, 1e-30),
                        self.nu)
        self.nu_wall = xp.where(w, xp.maximum(nu_w, self.nu), self.nu)
        self.u_tau_wall = xp.where(w, u_tau, 0.0)

    def _wall_cell_shear(self):
        xp, fvm = self.xp, self.fvm
        w = self.is_wall
        ut = self.u_tau_wall
        sw = ut ** 2 / (self.nu + self.KAPPA * ut * fvm.dperp)
        cnt = xp.bincount(fvm.Pb, weights=xp.where(w, 1.0, 0.0), minlength=fvm.nc)
        val = xp.bincount(fvm.Pb, weights=xp.where(w, sw, 0.0), minlength=fvm.nc)
        return cnt > 0, val / xp.maximum(cnt, 1.0)

    def _omega_wall_cells(self):
        """ω imposé dans les cellules pariétales (omegaWallFunction d'OpenFOAM) :
        ω = √(ω_vis² + ω_log²), ω_vis = 6ν/(β₁y²), ω_log = u_τ/(√β* κ y)."""
        xp, fvm = self.xp, self.fvm
        w = self.is_wall
        y = fvm.dperp
        beta1 = getattr(self.model, "beta1", getattr(self.model, "beta", 0.075))
        bstar = getattr(self.model, "beta_star", 0.09)
        om = xp.sqrt((6 * self.nu / (beta1 * y ** 2)) ** 2
                     + (self.u_tau_wall / (np.sqrt(bstar) * self.KAPPA * y)) ** 2)
        cnt = xp.bincount(fvm.Pb, weights=xp.where(w, 1.0, 0.0), minlength=fvm.nc)
        val = xp.bincount(fvm.Pb, weights=xp.where(w, om, 0.0), minlength=fvm.nc)
        mask = cnt > 0
        return mask, val / xp.maximum(cnt, 1.0)

    def cell_force(self, face=False):
        """Force volumique par unité de masse aux cellules (ou aux faces internes) :
        force imposée + flottabilité de Boussinesq −β(T − T_ref) g."""
        xp = self.xp
        n = self.fvm.ni if face else self.fvm.nc
        f = xp.zeros((n, 2)) + self.body_force_at(self._t_eval)[None, :]
        e = self.energy
        if e is not None and e["beta"]:
            T = self.fvm.interp(self.T) if face else self.T
            g = xp.asarray(np.asarray(e["gravity"], dtype=float))
            f = f - float(e["beta"]) * (T - float(e["T_ref"]))[:, None] * g[None, :]
        return f

    def temperature_bc(self):
        xp = self.xp
        dp = self.fvm.dperp
        fx = self.kindTemp == 0
        a = xp.where(fx, 0.0, 1.0)
        b = xp.where(fx, self.T_fixed, dp * self.T_grad)
        g = xp.where(fx, -1.0 / dp, 0.0)
        d = xp.where(fx, self.T_fixed / dp, self.T_grad)
        return a, b, g, d

    def boundary_T(self, T):
        a, b, _, _ = self.temperature_bc()
        return a * T[self.fvm.Pb] + b

    def _energy_eq(self, a0=0.0, hist=None, relax=1.0):
        """Équation de la température : ∂T/∂t + ∇·(UT) = ∇·((ν/Pr + ν_t/Pr_t)∇T) + Q."""
        if self.energy is None:
            return
        xp = self.xp
        fvm, s, e = self.fvm, self.settings, self.energy
        alpha = self.nu / float(e["Pr"])
        a_eff = alpha + self.nut / float(e["Pr_t"])
        gam_i = fvm.interp(a_eff)
        gam_b = xp.where(self.is_wall, alpha, a_eff[fvm.Pb])
        bc = self.temperature_bc()
        grad = fvm.grad(self.T, bc[0] * self.T[fvm.Pb] + bc[1])
        diag, up, lo, rhs = fvm.assemble(self.F_i, self.F_b, gam_i, gam_b, bc, grad_phi=grad,
                                         scheme=s.convection_T, bounded=self.steady,
                                         phi=self.T, nonorth_limit=s.nonorth_limit)
        if e.get("source"):
            rhs = rhs + float(e["source"]) * fvm.V
        if _active(a0):
            diag = diag + a0 * fvm.V
            rhs = rhs - hist * fvm.V
        if relax < 1.0:
            rhs = rhs + (1.0 - relax) / relax * diag * self.T
            diag = diag / relax
        A = fvm.matrix(diag, up, lo)
        dT = max(float(xp.max(xp.abs(self.T_fixed))) if bool(xp.any(self.kindTemp == 0))
                 else 0.0, 1.0)
        self.residuals_now["T"] = normalized_residual(A, self.T, rhs, dT)
        self.T = fvm.lin.solve(A, rhs, self.T, s.solver_U, rtol=0.1 if self.steady else 1e-6,
                               tag="T")

    def _momentum(self, a0=0.0, hist=None, relax=1.0, theta=1.0, explicit=None):
        """Équations de U (composantes x, y) : (diag, upper, lower, rhs) sans gradient de p.

        theta < 1 (Crank-Nicolson) : opérateur spatial pondéré par θ, plus la partie
        explicite `explicit` (= (1−θ)·R(Uⁿ)·V, tableau (nc, 2)).
        """
        xp = self.xp
        fvm, U = self.fvm, self.U
        force = self.cell_force()
        self._update_wall_function()
        nu_eff = self.nu + self.nut
        gam_i = fvm.interp(nu_eff)
        gam_b = xp.where(self.is_wall, self.nu_wall, nu_eff[fvm.Pb])
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
            rhs += force[:, c] * V
            if c == 1 and fvm.axisymmetric:
                # −τ_θθ/r = −2 ν_eff u_r/r² (les autres termes de ∇·τ sont donnés exactement
                # par la divergence pondérée par r des deux termes ci-dessus)
                diag = diag + 2.0 * nu_eff * V / fvm.radius ** 2
            if theta != 1.0:
                diag, up, lo, rhs = theta * diag, theta * up, theta * lo, theta * rhs
            if explicit is not None:
                rhs = rhs + explicit[:, c]
            if _active(a0):
                diag = diag + a0 * V
                rhs -= hist[:, c] * V
            if relax < 1.0:
                rhs += (1.0 - relax) / relax * diag * U[:, c]
                diag = diag / relax
            eqs.append((diag, up, lo, rhs))
        return eqs

    # ------------------------------------------------------------------ pression
    def _pressure_correction(self, eqs, ddt_corr=None, relax_p=1.0, rtol=1e-6):
        xp = self.xp
        fvm, s = self.fvm, self.settings
        V = fvm.V
        H = xp.column_stack([rhs - self._offdiag_mult(up, lo, self.U[:, c])
                             for c, (diag, up, lo, rhs) in enumerate(eqs)])
        diag_u, diag_v = eqs[0][0], eqs[1][0]
        aP = 0.5 * (diag_u + diag_v)
        rAU = V / aP
        # diagonale commune a_P (moyenne des composantes, cmptAv d'OpenFOAM) ; l'écart propre
        # à chaque composante (symétrie, terme circonférentiel) passe dans H, pour que
        # U = HbyA − rAU ∇p vérifie exactement l'équation de quantité de mouvement
        HbyA = xp.column_stack([(H[:, c] - (eqs[c][0] - aP) * self.U[:, c]) / aP
                                for c in range(2)])
        pb = self.boundary_p(self.p)
        gradp = fvm.grad(self.p, pb)
        HbyA_b = self.boundary_U(HbyA)
        phiHbyA_i = xp.sum(fvm.interp(HbyA) * fvm.Si, axis=1)
        phiHbyA_b = xp.sum(HbyA_b * fvm.Sb, axis=1)
        if self._has_force():
            # force évaluée AUX FACES (comme le gradient de pression) au lieu de
            # l'interpolation de la force des cellules : l'équilibre hydrostatique
            # (flottabilité ⇄ ∇p) est représenté sans courants parasites.
            fc, ff = self.cell_force(), self.cell_force(face=True)
            phiHbyA_i = phiHbyA_i + xp.sum(
                (fvm.interp(rAU)[:, None] * ff - fvm.interp(rAU[:, None] * fc)) * fvm.Si, axis=1)
            # faces à gradient de p imposé (= f·n) : flux de force compensé → débit inchangé
            fb = xp.sum(self._boundary_force() * fvm.Sb, axis=1)
            phiHbyA_b = phiHbyA_b + xp.where(self.kindP == 1, rAU[fvm.Pb] * fb, 0.0)
        rAtU = rAU
        if s.algorithm.upper() == "SIMPLEC" and self.steady:
            H1 = -0.5 * (xp.bincount(fvm.P, eqs[0][1], fvm.nc) + xp.bincount(fvm.N, eqs[0][2], fvm.nc)
                         + xp.bincount(fvm.P, eqs[1][1], fvm.nc) + xp.bincount(fvm.N, eqs[1][2], fvm.nc))
            rAtU = V / xp.maximum(aP - H1, 1e-3 * aP)
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
            diag, up, lo, rhs = fvm.assemble(xp.zeros(fvm.ni), xp.zeros(fvm.nb), gam_i, gam_b,
                                             bc, grad_phi=gp_k, phi=p_new,
                                             nonorth_limit=s.nonorth_limit)
            rhs -= div_hbya
            if not self.has_fixed_p:
                diag[0] += diag[0]
            A = fvm.matrix(diag, up, lo)
            if k == 0:
                self.residuals_now["p"] = normalized_residual(A, self.p, rhs, self.u_scale ** 2)
            p_new = fvm.lin.solve(A, rhs, p_new, s.solver_p, rtol=rtol, symmetric=True,
                                  tag="p")
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
    def run_steady(self, max_iter=None, tol=None, verbose=False, log_every=50, callback=None,
                   probes=None):
        xp = self.xp
        s = self.settings
        max_iter = max_iter or s.max_iter
        tol = tol or s.tol
        self.steady = True
        relax_p = s.relax_p if s.algorithm.upper() == "SIMPLEC" else min(s.relax_p, 0.3)
        self.update_nut()
        t0 = time.perf_counter()
        converged = False
        forces = []
        for it in range(1, max_iter + 1):
            self.residuals_now = {}
            if s.pseudo_dt or s.pseudo_cfl:
                # Pseudo-transitoire. La sous-relaxation implicite équivaut à un pas de
                # pseudo-temps local ∝ V/a_P, donc ∝ Δy²/ν près des parois : convergence en
                # O(N²) itérations sur les maillages fins et étirés. pseudo_cfl : pas local
                # fondé sur le seul Courant CONVECTIF (pas local de SU2 / pseudo-transient
                # de Fluent) ; pseudo_dt : pas global.
                if s.pseudo_cfl:
                    _, _, conv, _ = self.courant(1.0)
                    a0 = conv / s.pseudo_cfl + 1e-6 * self.U_ref / max(
                        float(np.sqrt(np.mean(self.mesh.cell_volumes))), 1e-300)
                    a0v = a0[:, None]
                else:
                    a0 = a0v = 1.0 / s.pseudo_dt
                eqs = self._momentum(a0=a0, hist=-a0v * self.U)
                a0T = {k: -a0 * v for k, v in self.state.items()}
            else:
                eqs = self._momentum(relax=s.relax_U)
            for c, (diag, up, lo, rhs) in enumerate(eqs):
                A = self.fvm.matrix(diag, up, lo)
                gp = self.fvm.grad(self.p, self.boundary_p(self.p))[:, c] * self.fvm.V
                b = rhs - gp
                self.residuals_now["Ux" if c == 0 else "Uy"] = normalized_residual(
                    A, self.U[:, c], b, self.u_scale)
                self.U[:, c] = self.fvm.lin.solve(A, b, self.U[:, c], s.solver_U, rtol=0.1,
                                                  tag="U")
            self._pressure_correction(eqs, relax_p=relax_p, rtol=0.01)
            if s.pseudo_dt or s.pseudo_cfl:
                self._turbulence(a0=a0, hist=a0T, relax=s.relax_turb)
                if self.energy is not None:
                    self._energy_eq(a0=a0, hist=-a0 * self.T)
            else:
                self._turbulence(relax=s.relax_turb)
                self._energy_eq(relax=s.relax_T)
            cont = float(xp.sum(xp.abs(self.fvm.div(self.F_i, self.F_b))) /
                         (self.u_scale * self.fvm.total_area + 1e-300))
            rec = {"iteration": self.iterations_total + it, **self.residuals_now,
                   "continuity": cont}
            self.history.append(rec)
            if not all(xp.isfinite(v) for v in rec.values()) or not xp.all(xp.isfinite(self.U)):
                raise FloatingPointError(f"Divergence à l'itération {it}.")
            if probes:
                rec.update(probes(self))            # sondes (NaN possible hors domaine)
            if verbose and (it % log_every == 0 or it == 1):
                print("  it %5d  " % it + "  ".join(f"{k}={v:.2e}" for k, v in rec.items()
                                                    if k != "iteration"))
            if s.monitor_tol:
                forces.append(xp.concatenate([f["total"] for f in self.forces().values()]))
                w = s.monitor_window
                if len(forces) > w:
                    ref = xp.maximum(xp.abs(forces[-1]), 1e-12 * self.u_scale ** 2)
                    if xp.max(xp.abs(forces[-1] - forces[-1 - w]) / ref) < s.monitor_tol:
                        converged = "forces"
            if callback and callback(self, it):
                break                                  # arrêt demandé (interface graphique)
            if max(v for k, v in self.residuals_now.items()) < tol:
                converged = True
            if converged:
                break
        self.converged = converged
        self.iterations = it
        self.iterations_total += it
        self.wall_time = time.perf_counter() - t0
        if verbose:
            how = " (efforts stabilisés)" if converged == "forces" else ""
            print(f"  {'convergé' if converged else 'NON convergé'}{how} en {it} itérations "
                  f"({self.wall_time:.1f} s)")
        return bool(converged)

    # ------------------------------------------------------------------ instationnaire
    def courant(self, dt: float | None = None):
        """(Co, Dn, taux_conv, taux_diff) : nombre de Courant convectif (définition
        OpenFOAM : 0.5·Σ|φ_f|/V·Δt) et nombre de diffusion Dn = Σ ν_f|S|²/(d·S)/V·Δt
        (maxima sur les cellules) ; les taux sont par unité de Δt."""
        xp = self.xp
        fvm = self.fvm
        dt = self.dt if dt is None else dt
        nu_eff = self.nu + self.nut
        nu_f = fvm.interp(nu_eff) * fvm.g
        nu_b = xp.where(self.kindU != 1, xp.where(self.is_wall, self.nu, nu_eff[fvm.Pb])
                        * fvm.magSb / fvm.dperp, 0.0)
        diff = (xp.bincount(fvm.P, nu_f, fvm.nc) + xp.bincount(fvm.N, nu_f, fvm.nc)
                + xp.bincount(fvm.Pb, nu_b, fvm.nc)) / fvm.V
        aF = xp.abs(self.F_i)
        conv = 0.5 * (xp.bincount(fvm.P, aF, fvm.nc) + xp.bincount(fvm.N, aF, fvm.nc)
                      + xp.bincount(fvm.Pb, xp.abs(self.F_b), fvm.nc)) / fvm.V
        return float(conv.max() * dt), float(diff.max() * dt), conv, diff

    def _stable_dt(self, info, dt_old):
        """Pas de temps respectant max_co (et la limite de diffusion des schémas explicites)."""
        xp = self.xp
        s = self.settings
        _, _, conv, diff = self.courant(1.0)
        if info["kind"] == "explicit":
            # combinaison des limites convective et diffusive du schéma
            rate = xp.max(conv / min(s.max_co, info["co_max"]) + diff / info["dn_max"])
        else:
            rate = conv.max() / s.max_co
        dt = 0.9 / rate if rate > 0 else s.max_dt
        if dt_old is not None:
            dt = min(dt, 1.2 * dt_old)                  # croissance lissée (comme OpenFOAM)
        return min(dt, s.max_dt)

    def auto_time_scheme(self, dt: float) -> str:
        """Choix automatique (études `microrans schemes` et cylindre Re = 100, README) :
        - RK3 explicite (3 à 5× moins cher que PIMPLE à précision égale en convection
          dominante) si le Δt demandé respecte sa stabilité et que la diffusion ne limite
          pas le pas ;
        - sinon Crank-Nicolson en laminaire (au même Δt, amplitude de portance du cylindre
          0.339 contre 0.389 en BDF2 ; valeur convergée ≈ 0.31-0.32) ;
        - BDF2 (L-stable) avec un modèle de turbulence (termes sources raides)."""
        xp = self.xp
        info = TIME_SCHEMES["rk3"]
        _, _, conv, diff = self.courant(1.0)
        dt_conv = info["co_max"] / max(conv.max(), 1e-300)
        dt_diff = info["dn_max"] / max(diff.max(), 1e-300)
        explicit_ok = xp.max(conv * dt / info["co_max"] + diff * dt / info["dn_max"]) <= 1.0
        if self.model.variables:
            return "backward"
        if dt_diff < 0.5 * dt_conv or not explicit_ok:
            return "crankNicolson"
        return "rk3"

    def run_transient(self, dt: float, t_end: float, verbose=False, log_every=50,
                      callback=None, probes=None):
        """Intégration en temps physique jusqu'à t_end.

        Schémas (Settings.time_scheme) : voir TIME_SCHEMES. Implicites (PIMPLE) : euler,
        backward (BDF2 à pas variable), crankNicolson. Explicites à projection : rk1, rk2,
        rk3 (SSP), rk4, ab2. Si Settings.adjust_dt, Δt suit max_co (et la limite de
        diffusion des schémas explicites). probes : fonction(solver) -> dict par pas.
        """
        xp = self.xp
        s = self.settings
        self.steady = False
        self.update_nut()
        name = canonical_time_scheme(s.time_scheme)
        if name == "auto":
            name = self.auto_time_scheme(dt)
            s.time_scheme = name
        info = TIME_SCHEMES[name]
        restart = self._restart_hist
        self._restart_hist = None
        if restart is not None and s.adjust_dt and self.dt > 0:
            self.dt = self._stable_dt(info, self.dt)     # suite exacte d'un calcul repris
        else:
            self.dt = float(dt)
            if s.adjust_dt:
                self.dt = min(self.dt, self._stable_dt(info, None))
        if not s.adjust_dt and info["kind"] == "explicit":
            co, dn, _, _ = self.courant(self.dt)
            if co / info["co_max"] + dn / info["dn_max"] > 1.0:
                import warnings
                warnings.warn(f"Δt = {self.dt:.3g} au-delà de la limite de stabilité de "
                              f"{name} (Co = {co:.2f}, Dn = {dn:.2f}, limites "
                              f"{info['co_max']:.2f} / {info['dn_max']:.2f}) : risque de "
                              f"divergence. Réduire Δt ou activer adjust_dt.")
        series = self.series = []
        t0 = time.perf_counter()
        n = 0
        step = self._pimple_step if info["kind"] == "implicit" else self._explicit_step
        self._hist = restart or {"U": [], "F": [], "state": [], "dt": [], "R": [], "T": []}
        eps = 1e-9 * max(abs(t_end), 1.0)
        while self.time < t_end - eps:
            n += 1
            if s.adjust_dt and n > 1:
                self.dt = self._stable_dt(info, self.dt)
            if s.adjust_dt or n == 1:
                # dernier pas raccourci / allongé pour tomber exactement sur t_end
                remaining = t_end - self.time
                if self.dt > remaining or remaining - self.dt < 0.05 * self.dt:
                    self.dt = remaining
            step(name, info)
            if not xp.all(xp.isfinite(self.U)):
                raise FloatingPointError(f"Divergence à t = {self.time:.4g}.")
            co, dn, _, _ = self.courant()
            rec = {"time": self.time, "dt": self.dt, "Co": co, **self.residuals_now}
            if info["kind"] == "explicit":
                rec["Dn"] = dn
            if probes:
                rec.update(probes(self))
            series.append(rec)
            if verbose and (n % log_every == 0 or n == 1):
                print(f"  t={self.time:.4f}  " + "  ".join(
                    f"{k}={v:.3e}" for k, v in rec.items() if k != "time"))
            if callback and callback(self, n):
                break
        self.wall_time = time.perf_counter() - t0
        self.time_steps = n
        return series

    # ---------------------------------------------------------------- PIMPLE (implicite)
    def _push_history(self):
        h = self._hist
        h["U"].insert(0, self.U.copy())
        h["F"].insert(0, self.F_i.copy())
        h["state"].insert(0, {k: v.copy() for k, v in self.state.items()})
        h["T"].insert(0, self.T.copy() if self.energy is not None else None)
        h["dt"].insert(0, self.dt)
        for key in ("U", "F", "state", "dt", "T"):
            del h[key][3:]

    def _pimple_step(self, name, info):
        xp = self.xp
        s, fvm = self.settings, self.fvm
        dt = self.dt
        if not self._hist["U"]:
            self._push_history()                      # niveau n
        h = self._hist
        U_old, F_old, st_old = h["U"][0], h["F"][0], h["state"][0]
        second = name == "backward" and len(h["U"]) >= 2
        theta, explicit = 1.0, None
        if name == "crankNicolson":
            # partie explicite (1−θ)·R(Uⁿ, tⁿ) évaluée avec les champs de l'instant n
            theta = s.cn_theta
            self._t_eval = self.time
            R0 = self._spatial_operator()
            explicit = (1.0 - theta) * R0
        if second:
            w = dt / h["dt"][1]                       # BDF2 à pas variable
            a0 = (1.0 + 2.0 * w) / ((1.0 + w) * dt)
            b1 = -(1.0 + w) / dt
            b2 = w * w / ((1.0 + w) * dt)
            U_old2, F_old2, st_old2 = h["U"][1], h["F"][1], h["state"][1]
            histU = b1 * U_old + b2 * U_old2
            histT = {k: b1 * st_old[k] + b2 * st_old2[k] for k in st_old}
            c1 = F_old - xp.sum(fvm.interp(U_old) * fvm.Si, axis=1)
            c2 = F_old2 - xp.sum(fvm.interp(U_old2) * fvm.Si, axis=1)
            ddt_corr = -self._ddt_coef(c1, F_old) * (b1 * c1 + b2 * c2)
        else:
            a0 = 1.0 / dt
            histU = -U_old / dt
            histT = {k: -v / dt for k, v in st_old.items()}
            c1 = F_old - xp.sum(fvm.interp(U_old) * fvm.Si, axis=1)
            ddt_corr = self._ddt_coef(c1, F_old) * c1 / dt
        if name == "crankNicolson" and len(h["U"]) >= 2:
            # turbulence : BDF2 (L-stable) même avec Crank-Nicolson sur la quantité de mvt
            w = dt / h["dt"][1]
            aT = (1.0 + 2.0 * w) / ((1.0 + w) * dt)
            b1, b2 = -(1.0 + w) / dt, w * w / ((1.0 + w) * dt)
            histT = {k: b1 * st_old[k] + b2 * h["state"][1][k] for k in st_old}
        else:
            aT = a0
        self.time += dt
        self._t_eval = self.time
        for outer in range(s.n_outer):
            self.residuals_now = {}
            final = outer == s.n_outer - 1
            eqs = self._momentum(a0=a0, hist=histU, relax=1.0, theta=theta, explicit=explicit)
            gp = fvm.grad(self.p, self.boundary_p(self.p)) * fvm.V[:, None]
            for c, (diag, up, lo, rhs) in enumerate(eqs):
                A = fvm.matrix(diag, up, lo)
                b = rhs - gp[:, c]
                self.residuals_now["Ux" if c == 0 else "Uy"] = normalized_residual(
                    A, self.U[:, c], b, self.u_scale)
                self.U[:, c] = fvm.lin.solve(A, b, self.U[:, c], s.solver_U, rtol=1e-4,
                                             tag="U")
            for corr in range(s.n_corr):
                last = final and corr == s.n_corr - 1
                self._pressure_correction(eqs, ddt_corr=ddt_corr, rtol=1e-6 if last else 0.01)
            if final or s.turbulence_every_outer:
                self._turbulence(a0=aT, hist=histT)
            if self.energy is not None:
                self._energy_eq(a0=aT, hist=self._T_hist(aT, name))
        self._push_history()

    def _ddt_coef(self, c1, F_old):
        """Coefficient de la correction ddtCorr du flux. 1 (défaut) : formulation cohérente
        de Tuković, Perić & Jasak (2018) — les flux aux faces gardent leur propre dérivée
        temporelle, résultat indépendant de Δt, ordre du schéma conservé. None : coefficient
        adaptatif historique d'OpenFOAM (1 − |c|/|φ|), qui rend la solution dépendante de Δt."""
        xp = self.xp
        k = self.settings.ddt_phi_coeff
        if k is None:
            return 1.0 - xp.minimum(xp.abs(c1) / (xp.abs(F_old) + 1e-30), 1.0)
        return k

    def _T_hist(self, a0, name):
        """Termes d'histoire de T (Euler ou BDF2 à pas variable, comme la turbulence)."""
        h = self._hist
        T1 = h["T"][0]
        if name in ("backward", "crankNicolson") and len(h["T"]) >= 2:
            w = self.dt / h["dt"][1]
            return -(1.0 + w) / self.dt * T1 + w * w / ((1.0 + w) * self.dt) * h["T"][1]
        return -T1 / self.dt

    def _spatial_operator(self):
        """R(U)·V = b − A·U : convection + diffusion + forces volumiques (sans pression),
        avec les schémas et conditions aux limites de l'équation implicite."""
        xp = self.xp
        U = self.U
        R = xp.empty_like(U)
        for c, (diag, up, lo, rhs) in enumerate(self._momentum()):
            R[:, c] = rhs - diag * U[:, c] - self._offdiag_mult(up, lo, U[:, c])
        return R

    # ---------------------------------------------------------------- projection (explicite)
    def _poisson(self):
        """Laplacien géométrique de pression (Γ = 1) : constant, factorisé une seule fois."""
        xp = self.xp
        if getattr(self, "_poisson_cache", None) is None:
            fvm = self.fvm
            bc = self.pressure_bc()
            diag, up, lo, rhs_bc = fvm.assemble(xp.zeros(fvm.ni), xp.zeros(fvm.nb),
                                                xp.ones(fvm.ni), xp.ones(fvm.nb), bc)
            if not self.has_fixed_p:
                diag = diag.copy()
                diag[0] += diag[0]
            L = fvm.matrix(diag, up, lo)
            if fvm.nc < 100000 and not self.backend.is_gpu:
                from scipy.sparse.linalg import splu
                lu = splu(L.tocsc(), permc_spec="MMD_AT_PLUS_A")
                solve = lambda b, x0: lu.solve(b)                       # noqa: E731
            else:
                from ..linalg import fcg
                M = fvm.lin.amg.setup(L)
                solve = lambda b, x0: fcg(L, b, x0, M, 1e-9, 500)[0]     # noqa: E731
            self._poisson_cache = (L, rhs_bc, solve)
        return self._poisson_cache

    def _project(self, Ustar, tau):
        """U* → champ à divergence (des flux) nulle : −∇·(τ∇p) = −∇·φ*."""
        xp = self.xp
        fvm = self.fvm
        L, rhs_bc, solve = self._poisson()
        Ub = self.boundary_U(Ustar)
        Fi = xp.sum(fvm.interp(Ustar) * fvm.Si, axis=1)
        Fb = xp.sum(Ub * fvm.Sb, axis=1)
        rhs = rhs_bc - fvm.sum_faces(Fi, Fb) / tau
        nf = 0.0
        if not fvm.orthogonal:
            # correction non orthogonale explicite avec la dernière pression connue
            gp = fvm.grad(self.p, self.boundary_p(self.p))
            nf = fvm.nonorth_flux(xp.ones(fvm.ni), gp, self.p, self.settings.nonorth_limit)
            rhs = rhs + xp.bincount(fvm.P, nf, fvm.nc) - xp.bincount(fvm.N, nf, fvm.nc)
        p = solve(rhs, self.p)
        _, _, g_, d_ = self.pressure_bc()
        self.F_i = Fi - tau * (fvm.g * (p[fvm.N] - p[fvm.P]) + nf)
        self.F_b = Fb - tau * fvm.magSb * (g_ * p[fvm.Pb] + d_)
        self.p = p
        self.U = Ustar - tau * fvm.grad(p, self.boundary_p(p))

    def _explicit_step(self, name, info):
        """Runge-Kutta explicite à projection à chaque étage (Sanderse & Koren 2012), ou
        Adams-Bashforth 2. La pression est le multiplicateur de Lagrange de chaque étage."""
        dt, V = self.dt, self.fvm.V[:, None]
        t_n = self.time
        U0, st_old = self.U.copy(), {k: v.copy() for k, v in self.state.items()}
        T0 = self.T.copy() if self.energy is not None else None
        if name == "ab2":
            self._t_eval = t_n
            R = self._spatial_operator() / V
            if self._hist["R"]:
                w = dt / self._hist["dt"][0]                 # AB2 à pas variable
                incr = (1.0 + 0.5 * w) * R - 0.5 * w * self._hist["R"][0]
            else:
                incr = R                                     # 1er pas : Euler explicite
            self._hist["R"] = [R]
            self._hist["dt"] = [dt]
            self._project(U0 + dt * incr, dt)
        else:
            A, b, c = RK_TABLES[name]
            F0 = (self.F_i.copy(), self.F_b.copy())
            Rs = []
            for i in range(len(b)):
                if i > 0:
                    self.F_i, self.F_b = F0
                    Ustar = U0 + dt * sum(A[i][j] * Rs[j] for j in range(i) if A[i][j])
                    self._project(Ustar, c[i] * dt)
                self._t_eval = t_n + c[i] * dt
                Rs.append(self._spatial_operator() / V)
            self.F_i, self.F_b = F0
            self._project(U0 + dt * sum(bj * Rj for bj, Rj in zip(b, Rs) if bj), dt)
        self.time = t_n + dt
        self.residuals_now = {}
        if self.model.variables:
            # turbulence : Euler implicite (termes sources raides) après le pas de U
            self._turbulence(a0=1.0 / dt, hist={k: -v / dt for k, v in st_old.items()})
        if self.energy is not None:
            self._energy_eq(a0=1.0 / dt, hist=-T0 / dt)

    # ------------------------------------------------------------------ matériel
    def to_cpu(self):
        """Rapatrie champs et opérateurs sur CPU (post-traitement après un calcul GPU)."""
        if not self.backend.is_gpu:
            return self
        h = self.backend.to_host
        for name in ("U", "p", "F_i", "F_b", "nut", "kindU", "U_fixed", "kindP", "p_fixed",
                     "kindT", "is_wall", "kindTemp", "T_fixed", "T_grad", "nu_wall",
                     "u_tau_wall"):
            setattr(self, name, h(getattr(self, name)))
        self.state = {k: h(v) for k, v in self.state.items()}
        if self.energy is not None:
            self.T = h(self.T)
        self.model.d = h(self.model.d)
        self.model.wall_nodes = h(self.model.wall_nodes)
        self.settings.backend = "cpu"
        self.fvm = FVM(self.mesh, "cpu", self.axisymmetric)
        self.backend, self.xp = self.fvm.backend, self.fvm.xp
        self._poisson_cache = None
        return self

    # ------------------------------------------------------------------ post-traitement
    def forces(self, patches=None):
        """Efforts (par unité d'envergure, ρ = 1) sur des patches : pression + frottement."""
        xp = self.xp
        fvm = self.fvm
        names = patches or [p.name for p in self.mesh.patches if p.type == "wall"]
        out = {}
        pb = self.boundary_p(self.p)
        for name in names:
            sl = self.patch_slices[name]
            Pb = fvm.Pb[sl]
            Fp = xp.sum(pb[sl, None] * fvm.Sb[sl], axis=0)
            du = self.U[Pb] - self.U_fixed[sl]
            Fv = xp.sum((self.nu_wall[sl] * fvm.magSb[sl] / fvm.dperp[sl])[:, None] * du, axis=0)
            out[name] = {"pressure": Fp, "viscous": Fv, "total": Fp + Fv}
        return out

    def moment(self, patch, center=(0.0, 0.0)):
        """Moment (par unité d'envergure, ρ = 1, sens trigonométrique : positif = cabrer
        pour un écoulement selon +x) des efforts de pression et de frottement sur `patch`
        autour de `center`."""
        xp = self.xp
        fvm = self.fvm
        sl = self.patch_slices[patch]
        pb = self.boundary_p(self.p)[sl]
        du = self.U[fvm.Pb[sl]] - self.U_fixed[sl]
        dF = pb[:, None] * fvm.Sb[sl] + (self.nu_wall[sl] * fvm.magSb[sl]
                                         / fvm.dperp[sl])[:, None] * du
        r = self.backend.asarray(self.mesh.face_centers[self.mesh.n_internal:][
            self.patch_slices[patch]] - np.asarray(center, float))
        return float(xp.sum(r[:, 0] * dF[:, 1] - r[:, 1] * dF[:, 0]))

    def wall_shear(self, patch):
        """(abscisse curviligne implicite) centres de faces, τ_w signé (tangente locale), y⁺."""
        xp = self.xp
        fvm = self.fvm
        sl = self.patch_slices[patch]
        Pb = fvm.Pb[sl]
        n = fvm.nb_hat[sl]
        t = xp.column_stack([-n[:, 1], n[:, 0]])
        du = self.U[Pb] - self.U_fixed[sl]
        ut = xp.sum(du * t, axis=1)
        tau = self.nu_wall[sl] * ut / fvm.dperp[sl]
        yplus = fvm.dperp[sl] * xp.sqrt(xp.abs(tau)) / self.nu
        xf = self.mesh.face_centers[self.mesh.n_internal:][sl]
        return xf, tau, yplus

    def mean_fields(self) -> dict:
        """Moyennes et écarts-types temporels (si [output] average_from est donné)."""
        return self.averager.fields() if self.averager is not None else {}

    def fields(self) -> dict:
        xp = self.xp
        out = {"U": self.U, "p": self.p, "U_mag": xp.linalg.norm(self.U, axis=1),
               "vorticity": self.flow().vorticity}
        out.update(self.mean_fields())
        if self.model.variables:
            out["nut_over_nu"] = self.nut / self.nu
            out.update(self.state)
        if self.energy is not None:
            out["T"] = self.T
        out["wall_distance"] = self.mesh.wall_distance
        return out

    def bulk_temperature(self, x):
        """Température de mélange Σ u_x T V / Σ u_x V de la tranche de cellules dont
        l'étendue en x contient chaque abscisse de `x` (conduites orientées selon x)."""
        xp = self.xp
        m = self.mesh
        nodes_x = np.where(m.cell_nodes >= 0, m.points[np.maximum(m.cell_nodes, 0), 0], np.nan)
        lo, hi = np.nanmin(nodes_x, axis=1), np.nanmax(nodes_x, axis=1)
        w = self.backend.to_host(self.U[:, 0] * self.fvm.V)
        T = self.backend.to_host(self.T)
        out = np.empty(len(x))
        for i, xi in enumerate(np.asarray(x, float)):
            sel = (lo <= xi) & (xi <= hi)
            out[i] = np.sum(w[sel] * T[sel]) / np.sum(w[sel]) if np.any(sel) else np.nan
        return xp.asarray(out) if self.backend.is_gpu else out

    def wall_heat_flux(self, patch):
        """(T_paroi, flux q entrant dans le fluide, en unités cinématiques q/(ρ c_p))."""
        fvm = self.fvm
        sl = self.patch_slices[patch]
        alpha = self.nu / float(self.energy["Pr"])
        Tb = self.boundary_T(self.T)[sl]
        return Tb, alpha * (Tb - self.T[fvm.Pb[sl]]) / fvm.dperp[sl]
