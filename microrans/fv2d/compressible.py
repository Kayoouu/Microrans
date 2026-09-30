"""Solveur compressible 2D en densité (Euler / Navier-Stokes laminaire), volumes finis
centrés sur les cellules, maillages non structurés polygonaux (mêmes maillages que le
solveur incompressible : Mesh2D, périodicité comprise).

Variables conservatives Q = (ρ, ρu, ρv, ρE), gaz parfait p = ρ R T, E = p/((γ−1)ρ) + ½|u|².
Grandeurs dimensionnelles SI (Pa, K, kg/m³, m/s).

Choix numériques (et pourquoi)
- Flux convectifs : solveur de Riemann approché de Roe (Roe 1981, J. Comput. Phys. 43 ;
  forme vectorielle « dissipation |A|ΔQ » de Blazek, « Computational Fluid Dynamics:
  Principles and Applications », 3e éd., 2015, éq. 4.89-4.91) avec correction d'entropie
  de Harten sur les ondes acoustiques (|λ| < δ → (λ² + δ²)/(2δ), δ = ε c̃, ε = 0.1 par
  défaut), ou HLLC (Toro, Spruce & Speares 1994 ; Toro 2009 § 10.4 ; vitesses d'onde
  d'Einfeldt / Roe, Batten et al. 1997), positif et robuste pour les chocs forts et les
  détentes vers le vide. Roe est le choix par défaut : contact et couches de cisaillement
  exacts (dissipation minimale en couche limite) ; HLLC en option.
- Ordre 2 : reconstruction MUSCL linéaire des variables PRIMITIVES (ρ, u, v, p) aux faces,
  φ_f = φ_P + ψ_P ∇φ_P·(x_f − x_P), gradient de Green-Gauss (même opérateur que fvm.grad,
  assemblé ici en matrice creuse pour traiter les 4-5 variables d'un coup), limiteur :
    * venkatakrishnan (défaut) : Venkatakrishnan 1995 (AIAA J. 33(1)), lissé donc
      favorable à la convergence stationnaire, avec le seuil ε² = (K Δφ_global)² de Wang
      (2000, AIAA 2000-0911 ; limiteur VENKATAKRISHNAN_WANG de SU2) : ε indépendant des
      unités et de la taille des mailles (grandeurs SI dimensionnelles ici) ;
    * barth_jespersen : Barth & Jespersen 1989 (min-max strict, fvm.limit_grad) ; très
      robuste mais non différentiable : les résidus stationnaires stagnent souvent ;
    * none.
  États reconstruits non physiques (ρ ≤ 0 ou p ≤ 0) : retour à l'ordre 1 sur la face.
  Option « gel » du limiteur après N itérations (limiter_freeze, comme LIMITER_ITER de SU2)
  pour éviter le cyclage limite des résidus stationnaires.
- Flux visqueux (Navier-Stokes laminaire) : gradient aux faces = moyenne des gradients de
  cellule + correction selon PN (« edge-normal correction », Blazek § 5.4 ; Weiss et al.
  1999), τ = μ(∇u + ∇uᵀ − ⅔ ∇·u I), q = −k ∇T, k = μ c_p / Pr ; μ constante ou loi de
  Sutherland (μ_ref = 1.716e-5 Pa·s, T_ref = 273.15 K, S = 110.4 K : White, « Viscous
  Fluid Flow », 3e éd., éq. 1-36).
- Temps :
    * instationnaire : Runge-Kutta SSP d'ordre 3 (Shu & Osher 1988), pas de temps GLOBAL
      Δt = CFL · min_i V_i / (Λc_i + C_v Λv_i), Λc_i = Σ_f (|u·n| + c)_f S_f,
      Λv_i = max(4/3, γ/Pr) (μ/ρ)_i Σ_f S_f² / V_i (rayons spectraux, Blazek § 6.1.4) ;
    * stationnaire : même schéma avec pas de temps LOCAL (Δt_i avec le même CFL), ou
      schéma multi-étapes à 5 étapes optimisé pour l'amortissement (coefficients de van
      Leer, Tai & Powell 1989, AIAA 89-1933, 2e ordre upwind) — « rk5 » — ou implicite
      (Euler implicite linéarisé, jacobienne d'ordre 1 de Rusanov assemblée en matrice
      creuse par blocs, résolue par GMRES + ILU ; CFL croissant) — « implicit ».
- Conditions aux limites (par état fantôme, puis le même solveur de Riemann à la face) :
    farfield           invariants de Riemann 1D normaux (Blazek § 8.4 / Jameson 1983) ;
                       entrée / sortie, subsonique / supersonique selon la vitesse normale
    supersonic_inlet   état imposé (écoulement amont par défaut, ou mach / p / T / angle)
    inlet              entrée subsonique : pression totale p0, température totale T0,
                       direction ; invariant sortant u_n + 2c/(γ−1) extrapolé (Blazek § 8.5)
    outlet             sortie subsonique : pression statique p ; entropie et invariant
                       sortant extrapolés (sortie supersonique locale : extrapolation)
    supersonic_outlet  extrapolation de tout l'état
    slip_wall          paroi glissante (Euler) : état miroir (u_n → −u_n)
    symmetry           comme slip_wall (et flux visqueux nul)
    wall               paroi adhérente adiabatique ; T = T_w : isotherme ; U = paroi
                       mobile ; en calcul non visqueux, « wall » = slip_wall
  Faces périodiques : faces internes du maillage (rien à faire).

Performance : NumPy vectorisé (aucune boucle Python sur les cellules ou les faces), CPU.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from ..mesh2d.mesh import Mesh2D
from .fvm import FVM

COMP_BC_TYPES = ("farfield", "supersonic_inlet", "inlet", "outlet", "supersonic_outlet",
                 "slip_wall", "symmetry", "wall")
_ALIASES = {"euler_wall": "slip_wall", "inviscid_wall": "slip_wall",
            "subsonic_inlet": "inlet", "total_pressure_inlet": "inlet",
            "pressure_inlet": "inlet", "subsonic_outlet": "outlet",
            "pressure_outlet": "outlet", "extrapolation": "supersonic_outlet",
            "no_slip_wall": "wall", "isothermal_wall": "wall", "adiabatic_wall": "wall"}
FLUXES = ("roe", "hllc")
LIMITERS = ("venkatakrishnan", "barth_jespersen", "none")
STEADY_SCHEMES = ("rk3", "rk5", "implicit")

# codes internes des conditions aux limites
_SLIP, _NOSLIP, _FAR, _SUPIN, _INLET, _OUTLET, _EXTRAP = range(7)
# van Leer, Tai & Powell (1989) : 5 étapes, discrétisation upwind d'ordre 2
_RK5 = (0.0695, 0.1602, 0.2898, 0.5060, 1.0)


# ============================================================================ gaz
@dataclass
class Gas:
    """Gaz parfait à γ constant ; viscosité : inviscid | constant | sutherland."""
    gamma: float = 1.4
    R: float = 287.058                    # air sec, J/(kg·K)
    Pr: float = 0.72
    viscosity: str = "inviscid"
    mu: float = 1.716e-5                  # constante, ou μ_ref de Sutherland (Pa·s)
    T_ref: float = 273.15                 # Sutherland (K)
    S: float = 110.4                      # Sutherland (K)

    def __post_init__(self):
        self.viscosity = str(self.viscosity).lower()
        if self.viscosity in ("euler", "none"):
            self.viscosity = "inviscid"
        if self.viscosity not in ("inviscid", "constant", "sutherland"):
            raise ValueError(f"Viscosité inconnue '{self.viscosity}' : inviscid | constant | "
                             "sutherland.")
        if self.gamma <= 1.0 or self.R <= 0.0 or self.Pr <= 0.0:
            raise ValueError("Gaz : γ > 1, R > 0 et Pr > 0 requis.")

    @property
    def cp(self) -> float:
        return self.gamma * self.R / (self.gamma - 1.0)

    @property
    def viscous(self) -> bool:
        return self.viscosity != "inviscid"

    def mu_of(self, T):
        if self.viscosity == "sutherland":
            T = np.asarray(T, float)
            return self.mu * (T / self.T_ref) ** 1.5 * (self.T_ref + self.S) / (T + self.S)
        if self.viscosity == "constant":
            return self.mu * np.ones_like(np.asarray(T, float))
        return np.zeros_like(np.asarray(T, float))


@dataclass
class State:
    """État uniforme (écoulement amont, état imposé) : ρ, u, v, p (+ T, c, Mach)."""
    rho: float
    u: float
    v: float
    p: float
    gas: Gas

    @property
    def T(self):
        return self.p / (self.rho * self.gas.R)

    @property
    def c(self):
        return float(np.sqrt(self.gas.gamma * self.p / self.rho))

    @property
    def speed(self):
        return float(np.hypot(self.u, self.v))

    @property
    def mach(self):
        return self.speed / self.c

    @property
    def p0(self):
        return self.p * (1.0 + 0.5 * (self.gas.gamma - 1.0) * self.mach ** 2) ** (
            self.gas.gamma / (self.gas.gamma - 1.0))

    @property
    def T0(self):
        return self.T + 0.5 * self.speed ** 2 / self.gas.cp

    def describe(self) -> dict:
        mu = float(self.gas.mu_of(self.T))
        return {"rho": self.rho, "u": self.u, "v": self.v, "p": self.p, "T": self.T,
                "c": self.c, "mach": self.mach, "speed": self.speed, "p0": self.p0,
                "T0": self.T0, "mu": mu}


def make_state(gas: Gas, mach=None, pressure=None, temperature=None, density=None,
               angle_deg: float = 0.0, velocity=None) -> State:
    """État uniforme défini par deux grandeurs parmi (p, T, ρ) et Mach (ou vitesse)."""
    given = [v is not None for v in (pressure, temperature, density)]
    if sum(given) < 2:
        raise ValueError("État : donner deux grandeurs parmi pressure, temperature, density.")
    R = gas.R
    if pressure is None:
        pressure = float(density) * R * float(temperature)
    elif density is None:
        density = float(pressure) / (R * float(temperature))
    p, rho = float(pressure), float(density)
    if p <= 0 or rho <= 0:
        raise ValueError("État : pression et masse volumique doivent être > 0.")
    c = np.sqrt(gas.gamma * p / rho)
    if velocity is not None:
        vel = np.asarray(velocity, float).reshape(-1)
        if vel.size == 2:
            return State(rho, float(vel[0]), float(vel[1]), p, gas)
        speed = float(vel[0])
    else:
        speed = float(mach or 0.0) * c
    a = np.radians(float(angle_deg))
    return State(rho, speed * np.cos(a), speed * np.sin(a), p, gas)


# ============================================================================ flux
def euler_flux(r, u, v, p, nx, ny, g):
    """Flux physique F(Q)·n (par unité de surface)."""
    qn = u * nx + v * ny
    m = r * qn
    H = g / (g - 1.0) * p / r + 0.5 * (u * u + v * v)
    return m, m * u + p * nx, m * v + p * ny, m * H


def roe_flux(WL, WR, nx, ny, g, eps: float = 0.1):
    """Flux de Roe (Blazek éq. 4.89-4.91) + correction d'entropie de Harten
    (ondes acoustiques, δ = eps·c̃). WL, WR : tuples (ρ, u, v, p) de tableaux."""
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
    qL, qR = uL * nx + vL * ny, uR * nx + vR * ny
    dr, dp, dq = rR - rL, pR - pL, qR - qL
    du, dv = uR - uL, vR - vL
    l1, l2, l3 = np.abs(q - c), np.abs(q), np.abs(q + c)
    if eps > 0.0:
        d = eps * c
        # (λ² + δ²)/(2δ) si |λ| < δ, écrit sans branchement : |λ| + max(δ − |λ|, 0)²/(2δ)
        h1, h3 = np.maximum(d - l1, 0.0), np.maximum(d - l3, 0.0)
        l1 = l1 + 0.5 * h1 * h1 / d
        l3 = l3 + 0.5 * h3 * h3 / d
    a1 = l1 * (dp - r * c * dq) / (2.0 * c2)
    a2 = l2 * (dr - dp / c2)
    a3 = l3 * (dp + r * c * dq) / (2.0 * c2)
    a4 = l2 * r
    D0 = a1 + a2 + a3
    D1 = a1 * (u - c * nx) + a2 * u + a3 * (u + c * nx) + a4 * (du - dq * nx)
    D2 = a1 * (v - c * ny) + a2 * v + a3 * (v + c * ny) + a4 * (dv - dq * ny)
    D3 = (a1 * (H - c * q) + a2 * k + a3 * (H + c * q)
          + a4 * (u * du + v * dv - q * dq))
    mL, mR = rL * qL, rR * qR
    return (0.5 * (mL + mR - D0),
            0.5 * (mL * uL + mR * uR + (pL + pR) * nx - D1),
            0.5 * (mL * vL + mR * vR + (pL + pR) * ny - D2),
            0.5 * (mL * HL + mR * HR - D3))


def hllc_flux(WL, WR, nx, ny, g, eps: float = 0.0):
    """Flux HLLC (Toro 2009 § 10.4), vitesses d'onde d'Einfeldt (moyennes de Roe)."""
    rL, uL, vL, pL = WL
    rR, uR, vR, pR = WR
    qL, qR = uL * nx + vL * ny, uR * nx + vR * ny
    cL, cR = np.sqrt(g * pL / rL), np.sqrt(g * pR / rR)
    sL, sR = np.sqrt(rL), np.sqrt(rR)
    iw = 1.0 / (sL + sR)
    HL = g / (g - 1.0) * pL / rL + 0.5 * (uL * uL + vL * vL)
    HR = g / (g - 1.0) * pR / rR + 0.5 * (uR * uR + vR * vR)
    ut = (sL * uL + sR * uR) * iw
    vt = (sL * vL + sR * vR) * iw
    Ht = (sL * HL + sR * HR) * iw
    ct = np.sqrt(np.maximum((g - 1.0) * (Ht - 0.5 * (ut * ut + vt * vt)), 1e-300))
    qt = ut * nx + vt * ny
    SL = np.minimum(qL - cL, qt - ct)
    SR = np.maximum(qR + cR, qt + ct)
    dL, dR = rL * (SL - qL), rR * (SR - qR)
    Sm = (pR - pL + qL * dL - qR * dR) / (dL - dR)
    EL = pL / ((g - 1.0) * rL) + 0.5 * (uL * uL + vL * vL)
    ER = pR / ((g - 1.0) * rR) + 0.5 * (uR * uR + vR * vR)
    FL = euler_flux(rL, uL, vL, pL, nx, ny, g)
    FR = euler_flux(rR, uR, vR, pR, nx, ny, g)
    # états étoile : U*_K = ρ_K (S_K − q_K)/(S_K − S*) [1, u_K + (S* − q_K) n, E_K/ρ_K + …]
    fL = dL / (SL - Sm)
    fR = dR / (SR - Sm)
    UL = (rL, rL * uL, rL * vL, rL * EL)
    UR = (rR, rR * uR, rR * vR, rR * ER)
    sUL = (fL, fL * (uL + (Sm - qL) * nx), fL * (vL + (Sm - qL) * ny),
           fL * (EL + (Sm - qL) * (Sm + pL / dL)))
    sUR = (fR, fR * (uR + (Sm - qR) * nx), fR * (vR + (Sm - qR) * ny),
           fR * (ER + (Sm - qR) * (Sm + pR / dR)))
    out = []
    for k in range(4):
        fsl = FL[k] + SL * (sUL[k] - UL[k])
        fsr = FR[k] + SR * (sUR[k] - UR[k])
        out.append(np.where(SL >= 0.0, FL[k], np.where(Sm >= 0.0, fsl,
                                                       np.where(SR > 0.0, fsr, FR[k]))))
    return tuple(out)


_FLUX_FUNCS = {"roe": roe_flux, "hllc": hllc_flux}


# ============================================================================ réglages
@dataclass
class CompressibleSettings:
    """Paramètres numériques du solveur compressible ([solver] du fichier de cas)."""
    flux: str = "roe"                     # roe | hllc
    order: int = 2                        # 1 | 2 (MUSCL)
    limiter: str = "venkatakrishnan"      # venkatakrishnan | barth_jespersen | none
    venkat_k: float = 0.05                # seuil de Wang : ε = K (max − min global)
    limiter_freeze: int = 0               # > 0 : limiteur gelé après cette itération
    entropy_fix: float = 0.1              # Harten : δ = entropy_fix · c̃ (Roe)
    cfl: float = 0.8
    steady_scheme: str = "rk3"            # rk3 (SSP, pas local) | rk5 | implicit
    cfl_max: float = 1e3                  # implicite : CFL atteint progressivement
    cfl_growth: float = 1.2               # implicite : facteur par itération
    linear_sweeps: int = 2                # implicite : balayages de Gauss-Seidel symétrique
    first_order_iter: int = 0             # itérations d'ordre 1 au départ (démarrage robuste)
    viscous_factor: float = 2.0           # C_v du pas de temps visqueux
    max_iter: int = 5000
    tol: float = 1e-6                     # chute relative des résidus (stationnaire)
    monitor_tol: float | None = None      # arrêt sur efforts stabilisés (relatif)
    monitor_window: int = 200
    log_every: int = 100

    def __post_init__(self):
        self.flux = str(self.flux).lower()
        self.limiter = str(self.limiter).lower().replace("-", "_")
        if self.limiter in ("venkat", "venkatakrishnan_wang"):
            self.limiter = "venkatakrishnan"
        if self.limiter in ("barth", "bj"):
            self.limiter = "barth_jespersen"
        self.steady_scheme = str(self.steady_scheme).lower()
        if self.flux not in FLUXES:
            raise ValueError(f"Flux inconnu '{self.flux}'. Choix : {FLUXES}")
        if self.limiter not in LIMITERS:
            raise ValueError(f"Limiteur inconnu '{self.limiter}'. Choix : {LIMITERS}")
        if self.steady_scheme not in STEADY_SCHEMES:
            raise ValueError(f"Schéma stationnaire inconnu '{self.steady_scheme}'. "
                             f"Choix : {STEADY_SCHEMES}")
        if int(self.order) not in (1, 2):
            raise ValueError("order : 1 ou 2.")
        self.order = int(self.order)


def canonical_bc(kind: str) -> str:
    k = str(kind).lower()
    k = _ALIASES.get(k, k)
    if k not in COMP_BC_TYPES:
        raise ValueError(f"Condition compressible inconnue '{kind}'. Choix : "
                         f"{', '.join(COMP_BC_TYPES)}")
    return k


# ============================================================================ solveur
class CompressibleSolver2D:
    """Solveur compressible (voir l'en-tête du module).

    mesh : Mesh2D ; gas : Gas ; freestream : State (écoulement de référence, sert aussi de
    valeur par défaut aux conditions aux limites) ; boundaries : {patch: {type, …}} ;
    initial : None (écoulement amont partout) ou fonction (x, y) -> (ρ, u, v, p).

    Rangement interne « par variable » : Q (4, nc) = (ρ, ρu, ρv, ρE), W (4, nc) =
    (ρ, u, v, p) ; chaque ligne est contiguë (opérations aux faces rapides). Les accesseurs
    publics (U, p, T, mach, fields…) renvoient des tableaux par cellule (nc, …).
    """

    def __init__(self, mesh: Mesh2D, gas: Gas, freestream: State, boundaries: dict,
                 settings: CompressibleSettings | None = None, initial=None):
        self.mesh = mesh
        self.gas = gas
        self.fs = freestream
        self.settings = settings or CompressibleSettings()
        self.fvm = FVM(mesh, "cpu")
        self.backend, self.xp = self.fvm.backend, np
        self.axisymmetric = False
        # attributs d'interface communs avec Solver2D (interface graphique, échantillonnage)
        self.energy = None
        self.rheology = None
        self.scalars: dict = {}
        self.state: dict = {}
        self.averager = None
        self.U_ref = max(freestream.speed, 1e-300)
        self.model_name = "euler" if not gas.viscous else "laminar"
        # échelles de référence de (ρ, u, v, p) (planchers du limiteur)
        vref = freestream.c + freestream.speed
        self._wref = np.array([freestream.rho, vref, vref, freestream.p])[:, None]
        self._geometry()
        self._setup_bc(boundaries)
        self._grad_ops()
        self._stencil()
        self.history: list[dict] = []
        self.series: list[dict] = []
        self.series_restart: list[dict] = []
        self.monitor: list[dict] = []
        self.iterations = 0
        self.iterations_total = 0
        self.time = 0.0
        self.dt = 0.0
        self.converged = False
        self.wall_time = 0.0
        self._res0 = None
        self._res0_count = 0
        self._psi = None                   # limiteur mémorisé (gel / réutilisation)
        self._psi_it = None
        self._last = {}
        C = mesh.cell_centers
        if initial is None:
            W = np.tile(np.array([freestream.rho, freestream.u, freestream.v,
                                  freestream.p])[:, None], (1, self.nc))
        else:
            vals = initial(C[:, 0], C[:, 1])
            W = np.vstack([np.broadcast_to(np.asarray(a, float), (self.nc,)) for a in vals])
        if np.any(W[0] <= 0) or np.any(W[3] <= 0) or not np.all(np.isfinite(W)):
            raise ValueError("État initial : ρ et p doivent être finis et > 0 partout.")
        self.Q = self.conservative(W)

    # ------------------------------------------------------------------ géométrie
    def _geometry(self):
        m, f = self.mesh, self.fvm
        self.nc, self.ni, self.nb = f.nc, f.ni, f.nb
        ni = self.ni
        self.P, self.N, self.Pb = m.owner[:ni], m.neighbour, m.owner[ni:]
        self.V = m.cell_volumes
        self.magS = m.magSf[:ni]
        self.nx, self.ny = m.nf[:ni, 0].copy(), m.nf[:ni, 1].copy()
        self.magSb = m.magSf[ni:]
        self.nbx, self.nby = m.nf[ni:, 0].copy(), m.nf[ni:, 1].copy()
        self.nb_hat = m.nf[ni:]
        rP, rN = np.asarray(f.rP), np.asarray(f.rN)
        rb = m.face_centers[ni:] - m.cell_centers[self.Pb]
        self.rPx, self.rPy = rP[:, 0].copy(), rP[:, 1].copy()
        self.rNx, self.rNy = rN[:, 0].copy(), rN[:, 1].copy()
        self.rbx, self.rby = rb[:, 0].copy(), rb[:, 1].copy()
        d = m.d_PN
        self.dist = np.linalg.norm(d, axis=1)
        self.ex, self.ey = d[:, 0] / self.dist, d[:, 1] / self.dist
        self.distb = np.linalg.norm(rb, axis=1)
        self.ebx, self.eby = rb[:, 0] / self.distb, rb[:, 1] / self.distb
        self.w = m.weights
        S2 = np.bincount(self.P, self.magS ** 2, self.nc) + np.bincount(
            self.N, self.magS ** 2, self.nc)
        if self.nb:
            S2 = S2 + np.bincount(self.Pb, self.magSb ** 2, self.nc)
        self.sumS2 = S2
        # somme des flux sortants : R = F_i @ D_i + F_b @ D_b (matrices creuses)
        ones = np.ones(ni)
        self.Di = sp.csr_matrix((np.concatenate([ones, -ones]),
                                 (np.concatenate([np.arange(ni), np.arange(ni)]),
                                  np.concatenate([self.P, self.N]))), shape=(ni, self.nc))
        self.Db = sp.csr_matrix((np.ones(self.nb), (np.arange(self.nb), self.Pb)),
                                shape=(self.nb, self.nc))

    def _grad_ops(self):
        """Opérateur de Green-Gauss (identique à FVM.grad) en matrices creuses
        (nc + nb, nc), appliquées à droite de [valeurs de cellules, valeurs frontières]
        rangées par variable (k, nc + nb)."""
        nc, ni, nb = self.nc, self.ni, self.nb
        P, N, Pb, w = self.P, self.N, self.Pb, self.w
        S = self.mesh.Sf
        V = self.V
        ops = []
        for c in range(2):
            s = S[:ni, c]
            rows = np.concatenate([P, P, N, N, Pb])
            cols = np.concatenate([P, N, P, N, nc + np.arange(nb)])
            vals = np.concatenate([w * s / V[P], (1 - w) * s / V[P],
                                   -w * s / V[N], -(1 - w) * s / V[N], S[ni:, c] / V[Pb]])
            ops.append(sp.csr_matrix((vals, (cols, rows)), shape=(nc + nb, nc)))
        self.GxT, self.GyT = ops

    def gradient(self, X, Xb):
        """Gradients de Green-Gauss de k variables : X (k, nc), Xb (k, nb) ->
        (gx, gy) chacun (k, nc)."""
        ext = np.hstack([X, Xb])
        return np.asarray(ext @ self.GxT), np.asarray(ext @ self.GyT)

    def _stencil(self):
        """Pochoir du limiteur : pour chaque cellule, ses voisines par face (cellule ou
        face frontière) ; « place » de chaque côté de face dans ce pochoir."""
        nc, ni, nb = self.nc, self.ni, self.nb
        cells = np.concatenate([self.P, self.N, self.Pb])
        other = np.concatenate([self.N, self.P, nc + np.arange(nb)])
        order = np.argsort(cells, kind="stable")
        counts = np.bincount(cells, minlength=nc)
        start = np.concatenate([[0], np.cumsum(counts)[:-1]])
        slot = np.empty(len(cells), dtype=np.int64)
        slot[order] = np.arange(len(cells)) - start[cells[order]]
        m = int(counts.max())
        idx = np.tile(np.arange(nc), (m, 1))
        idx[slot, cells] = other
        self._st_idx = idx
        self._st_m = m
        self._slots = (slot[:ni], slot[ni:2 * ni], slot[2 * ni:])

    # ------------------------------------------------------------------ conversions
    def conservative(self, W):
        g = self.gas.gamma
        r, u, v, p = W
        return np.vstack([r, r * u, r * v, p / (g - 1.0) + 0.5 * r * (u * u + v * v)])

    def primitive(self, Q=None):
        Q = self.Q if Q is None else Q
        g = self.gas.gamma
        r = Q[0]
        u, v = Q[1] / r, Q[2] / r
        p = (g - 1.0) * (Q[3] - 0.5 * r * (u * u + v * v))
        return np.vstack([r, u, v, p])

    # ------------------------------------------------------------------ conditions
    def _setup_bc(self, boundaries: dict):
        m, gas, fs = self.mesh, self.gas, self.fs
        nb = self.nb
        kind = np.full(nb, -1, dtype=int)
        self.bc_state = np.tile(np.array([fs.rho, fs.u, fs.v, fs.p])[:, None], (1, nb))
        self.bc_p = np.full(nb, fs.p)
        self.bc_p0 = np.full(nb, fs.p0)
        self.bc_T0 = np.full(nb, fs.T0)
        spd = max(fs.speed, 1e-300)
        d0 = [fs.u / spd, fs.v / spd] if fs.speed > 0 else [1.0, 0.0]
        self.bc_dir = np.tile(np.array(d0)[:, None], (1, nb))
        self.bc_Tw = np.full(nb, np.nan)                     # nan : paroi adiabatique
        self.bc_Uw = np.zeros((2, nb))
        self.patch_slices = {}
        self.bc_types = {}
        periodic = {q for pair in m.periodic_pairs for q in pair[:2]}
        names = {p.name for p in m.patches}
        missing = [p.name for p in m.patches if p.name not in boundaries and p.type != "empty"]
        if missing:
            raise ValueError(f"Conditions aux limites manquantes pour : {missing} "
                             f"(types : {', '.join(COMP_BC_TYPES)})")
        for name in boundaries:
            if name not in names and name not in periodic:
                raise ValueError(f"[boundary.{name}] : patch absent du maillage "
                                 f"({sorted(names)}).")
        C = m.face_centers[self.ni:]
        for patch in m.patches:
            sl = slice(patch.start - self.ni, patch.start - self.ni + patch.size)
            self.patch_slices[patch.name] = sl
            spec = dict(boundaries.get(patch.name, {"type": "slip_wall"}))
            if patch.type == "empty":
                spec = {"type": "slip_wall"}
            bkind = canonical_bc(spec.get("type", "wall"))
            if bkind == "wall" and not gas.viscous:
                bkind = "slip_wall"                  # Euler : paroi = glissement
            self.bc_types[patch.name] = bkind
            kind[sl] = {"slip_wall": _SLIP, "symmetry": _SLIP, "wall": _NOSLIP,
                        "farfield": _FAR, "supersonic_inlet": _SUPIN, "inlet": _INLET,
                        "outlet": _OUTLET, "supersonic_outlet": _EXTRAP}[bkind]
            x, y = C[sl, 0], C[sl, 1]
            if bkind in ("farfield", "supersonic_inlet"):
                st = _state_from_spec(spec, fs)
                self.bc_state[:, sl] = np.array([st.rho, st.u, st.v, st.p])[:, None]
            if bkind == "inlet":
                self.bc_p0[sl] = _value(spec.get("p0", spec.get("total_pressure", fs.p0)),
                                        x, y)
                self.bc_T0[sl] = _value(spec.get("T0", spec.get("total_temperature",
                                                                fs.T0)), x, y)
                if "direction" in spec or "angle" in spec:
                    if "direction" in spec:
                        dd = np.asarray(spec["direction"], float)
                    else:
                        a = np.radians(float(spec["angle"]))
                        dd = np.array([np.cos(a), np.sin(a)])
                    self.bc_dir[:, sl] = (dd / np.linalg.norm(dd))[:, None]
            if bkind == "outlet":
                self.bc_p[sl] = _value(spec.get("p", spec.get("pressure", fs.p)), x, y)
            if bkind == "wall":
                if "T" in spec or "temperature" in spec:
                    self.bc_Tw[sl] = _value(spec.get("T", spec.get("temperature")), x, y)
                if "U" in spec:
                    self.bc_Uw[0, sl] = _value(spec["U"][0], x, y)
                    self.bc_Uw[1, sl] = _value(spec["U"][1], x, y)
        if np.any(kind < 0):
            raise ValueError("Faces frontières sans condition aux limites.")
        self.kind = kind
        self._bc = {c: np.nonzero(kind == c)[0] for c in range(7)}
        self._bc = {c: i for c, i in self._bc.items() if len(i)}
        self.is_wall = np.isin(kind, (_SLIP, _NOSLIP))
        self.is_noslip = kind == _NOSLIP
        self.is_iso = self.is_noslip & np.isfinite(self.bc_Tw)
        self._noslip = np.nonzero(self.is_noslip)[0]
        self._adiab = np.nonzero(self.is_noslip & ~self.is_iso)[0]
        self._iso = np.nonzero(self.is_iso)[0]
        self._slip = np.nonzero(kind == _SLIP)[0]

    def ghost(self, W):
        """État fantôme (4, nb) des faces frontières à partir de l'état intérieur W (4, nb)
        (reconstruit à la face ou valeur de cellule)."""
        g = self.gas.gamma
        gm1 = g - 1.0
        nx, ny = self.nbx, self.nby
        G = W.copy()
        for code, s in self._bc.items():
            if code == _EXTRAP:
                continue
            r, u, v, p = W[:, s]
            if code == _SLIP:                                  # miroir
                qn = u * nx[s] + v * ny[s]
                G[1, s] = u - 2.0 * qn * nx[s]
                G[2, s] = v - 2.0 * qn * ny[s]
            elif code == _NOSLIP:
                G[1, s] = 2.0 * self.bc_Uw[0, s] - u
                G[2, s] = 2.0 * self.bc_Uw[1, s] - v
            elif code == _SUPIN:
                G[:, s] = self.bc_state[:, s]
            elif code == _FAR:
                G[:, s] = self._farfield(W[:, s], self.bc_state[:, s], nx[s], ny[s])
            elif code == _INLET:
                qn = u * nx[s] + v * ny[s]
                c = np.sqrt(g * p / r)
                Rp = qn + 2.0 * c / gm1                    # invariant sortant
                d = self.bc_dir[:, s]
                cs = np.maximum(-(d[0] * nx[s] + d[1] * ny[s]), 0.05)
                T0 = self.bc_T0[s]
                H0 = self.gas.cp * T0
                a = 0.25 * gm1 * cs * cs + 0.5
                b = 0.5 * gm1 * Rp * cs
                cc = 0.25 * gm1 * Rp * Rp - H0
                Vm = np.maximum((-b + np.sqrt(np.maximum(b * b - 4 * a * cc, 0.0)))
                                / (2 * a), 0.0)
                Tb = np.maximum(T0 - 0.5 * Vm * Vm / self.gas.cp, 1e-3 * T0)
                pb = self.bc_p0[s] * (Tb / T0) ** (g / gm1)
                G[0, s] = pb / (self.gas.R * Tb)
                G[1, s] = Vm * d[0]
                G[2, s] = Vm * d[1]
                G[3, s] = pb
            elif code == _OUTLET:
                qn = u * nx[s] + v * ny[s]
                c = np.sqrt(g * p / r)
                pb = np.where(qn < c, self.bc_p[s], p)     # supersonique : extrapolation
                rb = r * (pb / p) ** (1.0 / g)
                cb = np.sqrt(g * pb / rb)
                dq = 2.0 * (c - cb) / gm1
                G[0, s] = rb
                G[1, s] = u + dq * nx[s]
                G[2, s] = v + dq * ny[s]
                G[3, s] = pb
        return G

    def _farfield(self, W, Winf, nx, ny):
        """Invariants de Riemann normaux (Blazek § 8.4 ; Jameson & Baker 1983)."""
        g = self.gas.gamma
        gm1 = g - 1.0
        r, u, v, p = W
        ri, ui, vi, pi = Winf
        c, ci = np.sqrt(g * p / r), np.sqrt(g * pi / ri)
        qn, qi = u * nx + v * ny, ui * nx + vi * ny
        Rp = qn + 2.0 * c / gm1                       # sortant (intérieur)
        Rm = qi - 2.0 * ci / gm1                      # entrant (amont)
        qb = 0.5 * (Rp + Rm)
        cb = 0.25 * gm1 * (Rp - Rm)
        out = qb > 0.0
        s = np.where(out, p / r ** g, pi / ri ** g)
        ut = np.where(out, u - qn * nx, ui - qi * nx)
        vt = np.where(out, v - qn * ny, vi - qi * ny)
        rb = (cb * cb / (g * s)) ** (1.0 / gm1)
        pb = rb * cb * cb / g
        G = np.vstack([rb, ut + qb * nx, vt + qb * ny, pb])
        # supersonique : tout de l'amont (entrée) ou de l'intérieur (sortie)
        sup_in = qi <= -ci
        sup_out = qn >= c
        G[:, sup_in] = Winf[:, sup_in]
        G[:, sup_out] = W[:, sup_out]
        return G

    def boundary_values(self, W):
        """Valeurs aux faces frontières (5, nb) de (ρ, u, v, p, T) pour les gradients :
        moyenne cellule / état fantôme ; parois adhérentes : vitesse de paroi, T_w."""
        Wc = W[:, self.Pb]
        Wb = 0.5 * (Wc + self.ghost(Wc))
        s = self._noslip
        if len(s):
            Wb[1:3, s] = self.bc_Uw[:, s]
            Wb[0, s] = Wc[0, s]
            Wb[3, s] = Wc[3, s]
        T = Wb[3] / (Wb[0] * self.gas.R)
        if len(self._iso):
            T[self._iso] = self.bc_Tw[self._iso]
        return np.vstack([Wb, T])

    # ------------------------------------------------------------------ limiteurs
    def _limiter(self, W, Wb, dP, dN, db):
        """Coefficient ψ (4, nc) du limiteur ; dP, dN, db : variations cellule → face
        non limitées (côtés owner, neighbour, frontière)."""
        lim = self.settings.limiter
        if lim == "none":
            return np.ones_like(W)
        ext = np.hstack([W, Wb])[:, self._st_idx]              # (4, m, nc)
        dmax = np.maximum(ext.max(axis=1), W) - W
        dmin = np.minimum(ext.min(axis=1), W) - W
        if lim == "venkatakrishnan":
            # seuil de Wang : ε² = (K (max − min global))² ; plancher pour une variable
            # uniforme : 1e-6 × échelle de référence
            rng = np.maximum(W.max(axis=1, keepdims=True) - W.min(axis=1, keepdims=True),
                             1e-6 * self._wref)
            eps2 = (self.settings.venkat_k * rng) ** 2

            def psi(d, cells):
                lo = np.take(dmin, cells, axis=1)
                D1 = lo + (d > 0.0) * (np.take(dmax, cells, axis=1) - lo)
                D1sq = D1 * D1 + eps2
                dD1 = d * D1
                return (D1sq + 2.0 * dD1) / (D1sq + 2.0 * d * d + dD1)
        else:                                                   # Barth-Jespersen
            def psi(d, cells):
                pos, neg = d > 0.0, d < 0.0
                return np.where(pos, dmax[:, cells] / np.where(pos, d, 1.0),
                                np.where(neg, dmin[:, cells] / np.where(neg, d, -1.0), 1.0))
        out = np.ones((4, self._st_m, self.nc))
        sP, sN, sB = self._slots
        out[:, sP, self.P] = psi(dP, self.P)
        out[:, sN, self.N] = psi(dN, self.N)
        if self.nb:
            out[:, sB, self.Pb] = psi(db, self.Pb)
        return np.minimum(out.min(axis=1), 1.0)

    # ------------------------------------------------------------------ résidu
    def residual(self, Q, order=None, it=None, store=False, reuse_limiter=False):
        """Résidu R(Q) (4, nc) = Σ_f (F_c − F_v)·S_f (sortant) ; dQ/dt = −R/V.

        it : itération courante (gel du limiteur) ; reuse_limiter : réutilise le limiteur
        de la 1re étape du pas (stationnaire : étapes suivantes du Runge-Kutta)."""
        s = self.settings
        g = self.gas.gamma
        order = s.order if order is None else order
        P, N, Pb = self.P, self.N, self.Pb
        W = self.primitive(Q)
        if not (np.all(W[0] > 0.0) and np.all(W[3] > 0.0) and np.all(np.isfinite(W))):
            bad = int(np.sum((W[0] <= 0) | (W[3] <= 0) | ~np.isfinite(W).all(axis=0)))
            raise FloatingPointError(f"État non physique (ρ ≤ 0 ou p ≤ 0) dans {bad} "
                                     "cellule(s) : réduire le CFL, démarrer à l'ordre 1 "
                                     "(first_order_iter) ou changer de flux (hllc).")
        viscous = self.gas.viscous
        Xb = self.boundary_values(W)
        if order == 2 or viscous:
            X = np.vstack([W, W[3] / (W[0] * self.gas.R)])
            gx, gy = self.gradient(X, Xb)
        if order == 2:
            gx4, gy4 = gx[:4], gy[:4]
            dP = _take(gx4, P) * self.rPx + _take(gy4, P) * self.rPy
            dN = _take(gx4, N) * self.rNx + _take(gy4, N) * self.rNy
            db = _take(gx4, Pb) * self.rbx + _take(gy4, Pb) * self.rby
            frozen = s.limiter_freeze and it is not None and it > s.limiter_freeze
            if (frozen or reuse_limiter) and self._psi is not None:
                psi = self._psi
            else:
                psi = self._limiter(W, Xb[:4], dP, dN, db)
                self._psi = psi
            WL = _take(W, P) + _take(psi, P) * dP
            WR = _take(W, N) + _take(psi, N) * dN
            Wf = _take(W, Pb) + _take(psi, Pb) * db
            bad = (WL[0] <= 0) | (WL[3] <= 0) | (WR[0] <= 0) | (WR[3] <= 0)
            if bad.any():
                WL[:, bad], WR[:, bad] = W[:, P[bad]], W[:, N[bad]]
            badb = (Wf[0] <= 0) | (Wf[3] <= 0)
            if badb.any():
                Wf[:, badb] = W[:, Pb[badb]]
        else:
            WL, WR, Wf = _take(W, P), _take(W, N), _take(W, Pb)
        flux = _FLUX_FUNCS[s.flux]
        Fi = np.vstack(flux(WL, WR, self.nx, self.ny, g, s.entropy_fix))
        Fb = np.vstack(flux(Wf, self.ghost(Wf), self.nbx, self.nby, g, s.entropy_fix))
        if viscous:
            Fvi, Fvb = self._viscous_flux(X, Xb, gx, gy)
            Fi -= Fvi
            Fb -= Fvb
        Fi *= self.magS
        Fb *= self.magSb
        R = np.asarray(Fi @ self.Di + Fb @ self.Db)
        if store:
            self._last = {"Fb": Fb, "W": W}
            if viscous:
                self._last["Fvb"] = Fvb
        return R

    def _viscous_flux(self, X, Xb, gx, gy):
        """Flux visqueux (4, ·) par unité de surface aux faces internes et frontières.
        X : (ρ, u, v, p, T) aux cellules ; seules u, v, T (lignes 1, 2, 4) servent."""
        P, N, Pb = self.P, self.N, self.Pb
        w = self.w
        sel = [1, 2, 4]
        phi, gxs, gys = X[sel], gx[sel], gy[sel]
        gax = w * _take(gxs, P) + (1 - w) * _take(gxs, N)
        gay = w * _take(gys, P) + (1 - w) * _take(gys, N)
        ex, ey = self.ex, self.ey
        phP, phN = _take(phi, P), _take(phi, N)
        corr = (phN - phP) / self.dist - (gax * ex + gay * ey)
        gfx, gfy = gax + corr * ex, gay + corr * ey
        uf = w * phP + (1 - w) * phN
        Fi = self._stress_flux(uf, gfx, gfy, self.nx, self.ny)
        # frontières : parois adhérentes (correction selon d_b) ; autres : gradient de cellule
        phib = Xb[sel]
        gbx, gby = _take(gxs, Pb), _take(gys, Pb)
        s = self._noslip
        if len(s):
            ebx, eby = self.ebx[s], self.eby[s]
            cb = (phib[:, s] - phi[:, Pb[s]]) / self.distb[s] - (gbx[:, s] * ebx
                                                                  + gby[:, s] * eby)
            gbx[:, s] += cb * ebx
            gby[:, s] += cb * eby
        a = self._adiab                                 # adiabatique : ∂T/∂n = 0
        if len(a):
            qn = gbx[2, a] * self.nbx[a] + gby[2, a] * self.nby[a]
            gbx[2, a] -= qn * self.nbx[a]
            gby[2, a] -= qn * self.nby[a]
        Fb = self._stress_flux(phib, gbx, gby, self.nbx, self.nby)
        if len(self._slip):
            Fb[:, self._slip] = 0.0
        return Fi, Fb

    def _stress_flux(self, uvT, gx, gy, nx, ny):
        gas = self.gas
        u, v, T = uvT
        ux, vx, Tx = gx
        uy, vy, Ty = gy
        mu = gas.mu_of(np.maximum(T, 1e-3))
        kcond = mu * gas.cp / gas.Pr
        div = ux + vy
        txx = mu * (2.0 * ux - 2.0 / 3.0 * div)
        tyy = mu * (2.0 * vy - 2.0 / 3.0 * div)
        txy = mu * (uy + vx)
        fx = txx * nx + txy * ny
        fy = txy * nx + tyy * ny
        fe = u * fx + v * fy + kcond * (Tx * nx + Ty * ny)
        return np.vstack([np.zeros_like(fx), fx, fy, fe])

    # ------------------------------------------------------------------ pas de temps
    def spectral_radius(self, Q):
        """Λ_i = Λc_i + C_v Λv_i (Blazek § 6.1.4) : Δt_i = CFL V_i / Λ_i."""
        g = self.gas.gamma
        W = self.primitive(Q)
        c = np.sqrt(g * W[3] / W[0])
        P, N, Pb = self.P, self.N, self.Pb
        un = 0.5 * np.abs((W[1, P] + W[1, N]) * self.nx + (W[2, P] + W[2, N]) * self.ny)
        lam = (un + 0.5 * (c[P] + c[N])) * self.magS
        L = np.bincount(P, lam, self.nc) + np.bincount(N, lam, self.nc)
        if self.nb:
            lb = (np.abs(W[1, Pb] * self.nbx + W[2, Pb] * self.nby) + c[Pb]) * self.magSb
            L = L + np.bincount(Pb, lb, self.nc)
        if self.gas.viscous:
            mu = self.gas.mu_of(W[3] / (W[0] * self.gas.R))
            Lv = max(4.0 / 3.0, g / self.gas.Pr) * mu / W[0] * self.sumS2 / self.V
            L = L + self.settings.viscous_factor * Lv
        return L

    def local_dt(self, Q, cfl):
        return cfl * self.V / self.spectral_radius(Q)

    # ------------------------------------------------------------------ intégration
    def _ssp_rk3(self, dt, order, it=None, R0=None, reuse=False):
        """Runge-Kutta SSP d'ordre 3 (Shu & Osher 1988) ; dt global ou local (nc,) ;
        R0 : résidu déjà calculé de l'état courant ; reuse : limiteur de la 1re étape."""
        Q0 = self.Q
        f = dt / self.V
        if R0 is None:
            R0 = self.residual(Q0, order, it)
        Q1 = Q0 - f * R0
        Q2 = 0.75 * Q0 + 0.25 * (Q1 - f * self.residual(Q1, order, it, reuse_limiter=reuse))
        self.Q = Q0 / 3.0 + 2.0 / 3.0 * (Q2 - f * self.residual(Q2, order, it,
                                                                  reuse_limiter=reuse))

    def _rk5(self, dt, order, it=None, R0=None):
        Q0 = self.Q
        f = dt / self.V
        Q = Q0
        for k, a in enumerate(_RK5):
            R = R0 if (k == 0 and R0 is not None) else self.residual(
                Q, order, it, reuse_limiter=k > 0)
            Q = Q0 - a * f * R
        self.Q = Q

    def residual_norms(self, R) -> np.ndarray:
        """Normes RMS des résidus par unité de volume (4 équations)."""
        return np.sqrt(np.mean((R / self.V) ** 2, axis=1))

    def run_steady(self, max_iter=None, tol=None, verbose=False, log_every=None,
                   callback=None, monitor=None):
        """Marche en pseudo-temps (pas local) jusqu'à la chute des résidus de `tol`
        (normalisés par leur maximum sur les 10 premières itérations)."""
        s = self.settings
        max_iter = int(max_iter or s.max_iter)
        tol = float(tol or s.tol)
        log_every = log_every or s.log_every
        t0 = time.perf_counter()
        converged = False
        cfl = s.cfl
        stepper = None
        if s.steady_scheme == "implicit":
            from .compressible_implicit import ImplicitStepper
            stepper = ImplicitStepper(self)
        forces_hist = []
        it = 0
        for it in range(1, max_iter + 1):
            gi = self.iterations_total + it
            order = 1 if gi <= s.first_order_iter else s.order
            R = self.residual(self.Q, order, gi)
            if stepper is not None:
                ok = stepper.step(R, cfl, order)
                cfl = min(cfl * s.cfl_growth, s.cfl_max) if ok else max(0.5 * cfl, 0.1)
            else:
                dt = self.local_dt(self.Q, cfl)
                if s.steady_scheme == "rk5":
                    self._rk5(dt, order, gi, R0=R)
                else:
                    self._ssp_rk3(dt, order, gi, R0=R, reuse=True)
            rel = self._relative(self.residual_norms(R))
            rec = {"iteration": gi, "rho": rel[0], "rhoU": rel[1], "rhoV": rel[2],
                   "rhoE": rel[3]}
            if not all(np.isfinite(v) for v in rel) or not np.all(np.isfinite(self.Q)):
                raise FloatingPointError(f"Divergence à l'itération {gi}.")
            self.history.append(rec)
            if monitor is not None and (it % 10 == 0 or it == 1):
                self.monitor.append({"iteration": gi, **monitor(self)})
            if verbose and (it % log_every == 0 or it == 1):
                print("  it %6d  " % gi + "  ".join(f"{k}={v:.2e}" for k, v in rec.items()
                                                     if k != "iteration")
                      + (f"  CFL={cfl:.3g}" if stepper is not None else ""))
            if s.monitor_tol and it % 10 == 0:
                self._last = {}
                fr = self.forces()
                forces_hist.append(np.concatenate([f["total"] for f in fr.values()])
                                   if fr else np.zeros(1))
                w = max(s.monitor_window // 10, 1)
                if len(forces_hist) > w and gi > s.first_order_iter + s.monitor_window:
                    ref = max(float(np.max(np.abs(forces_hist[-1]))), 1e-300)
                    if np.max(np.abs(forces_hist[-1] - forces_hist[-1 - w])) / ref \
                            < s.monitor_tol:
                        converged = "forces"
            if callback and callback(self, it):
                break
            if max(rel) < tol:
                converged = True
            if converged:
                break
        self.residual(self.Q, s.order, self.iterations_total + it, store=True)
        self.converged = converged
        self.iterations = it
        self.iterations_total += it
        self.wall_time = time.perf_counter() - t0
        if verbose:
            how = " (efforts stabilisés)" if converged == "forces" else ""
            print(f"  {'convergé' if converged else 'NON convergé'}{how} en {it} itérations "
                  f"({self.wall_time:.1f} s)")
        return bool(converged)

    def _relative(self, norms):
        """Résidus normalisés par le maximum observé sur les 10 premières itérations ;
        quantité de mouvement : échelle commune aux deux composantes."""
        if self._res0 is None:
            self._res0 = np.zeros(4)
        if self._res0_count < 10:
            self._res0 = np.maximum(self._res0, norms)
            self._res0[1:3] = max(self._res0[1], self._res0[2])
            self._res0_count += 1
        return [float(v) for v in norms / np.maximum(self._res0, 1e-300)]

    def run_transient(self, t_end, dt=None, cfl=None, verbose=False, log_every=None,
                      callback=None, probes=None, max_steps=10 ** 8):
        """Instationnaire : SSP-RK3, pas global Δt = CFL · min(V/Λ) (ou Δt fixé)."""
        s = self.settings
        cfl = cfl or s.cfl
        log_every = log_every or s.log_every
        t0 = time.perf_counter()
        n = 0
        series = self.series = []
        while self.time < t_end * (1.0 - 1e-12) and n < max_steps:
            n += 1
            step = float(dt) if dt else float(np.min(self.local_dt(self.Q, cfl)))
            step = min(step, t_end - self.time)
            order = 1 if self.iterations_total + n <= s.first_order_iter else s.order
            self._ssp_rk3(step, order)
            self.time += step
            self.dt = step
            if not np.all(np.isfinite(self.Q)):
                raise FloatingPointError(f"Divergence au pas {n} (t = {self.time:.4g}).")
            rec = {"time": self.time, "dt": step}
            if probes:
                rec.update(probes(self))
            series.append(rec)
            if verbose and (n % log_every == 0 or n == 1):
                print(f"  pas {n:6d}  t = {self.time:.5g}  Δt = {step:.3g}")
            if callback and callback(self, n):
                break
        self.residual(self.Q, s.order, None, store=True)
        self.iterations = n
        self.iterations_total += n
        self.wall_time = time.perf_counter() - t0
        if verbose:
            print(f"  t = {self.time:.5g} atteint en {n} pas ({self.wall_time:.1f} s)")
        return series

    # ------------------------------------------------------------------ post-traitement
    @property
    def W(self):
        """Variables primitives par cellule (nc, 4) : ρ, u, v, p."""
        return self.primitive().T.copy()

    @property
    def U(self):
        return self.primitive()[1:3].T.copy()

    @property
    def p(self):
        return self.primitive()[3].copy()

    @property
    def rho(self):
        return self.Q[0].copy()

    @property
    def T(self):
        W = self.primitive()
        return W[3] / (W[0] * self.gas.R)

    @property
    def mach(self):
        W = self.primitive()
        return np.hypot(W[1], W[2]) / np.sqrt(self.gas.gamma * W[3] / W[0])

    def _ensure_last(self):
        if not self._last:
            self.residual(self.Q, self.settings.order, None, store=True)

    def boundary_p(self, p=None):
        """Pression aux faces frontières (nb,) : sur les parois, pression du flux de
        Riemann (quantité de mouvement normale transmise) ; ailleurs moyenne cellule /
        état fantôme."""
        self._ensure_last()
        pb = self.boundary_values(self.primitive())[3].copy()
        w = self.is_wall
        if w.any():
            Fc = self._last["Fb"].copy()
            if "Fvb" in self._last:
                Fc += self._last["Fvb"] * self.magSb
            pw = (Fc[1] * self.nbx + Fc[2] * self.nby) / self.magSb
            pb[w] = pw[w]
        return pb

    def grad_U(self, U=None):
        W = self.primitive()
        Xb = self.boundary_values(W)
        gx, gy = self.gradient(W[1:3], Xb[1:3])
        out = np.empty((self.nc, 2, 2))
        out[:, :, 0], out[:, :, 1] = gx.T, gy.T
        return out

    def _viscous_traction(self, sl):
        """−τ·S sur les faces `sl` (effort visqueux exercé sur la paroi), (n, 2)."""
        if "Fvb" not in self._last:
            return np.zeros((len(self.magSb[sl]), 2))
        return -(self._last["Fvb"][1:3, sl] * self.magSb[sl]).T

    def forces(self, patches=None):
        """Efforts sur des patches (N/m d'envergure) : pression (p − p∞) et frottement."""
        names = patches or [p.name for p in self.mesh.patches
                            if self.bc_types.get(p.name) in ("wall", "slip_wall")]
        self._ensure_last()
        pb = self.boundary_p()
        Sb = self.mesh.Sf[self.ni:]
        out = {}
        for name in names:
            sl = self.patch_slices[name]
            Fp = np.sum((pb[sl] - self.fs.p)[:, None] * Sb[sl], axis=0)
            Fv = np.sum(self._viscous_traction(sl), axis=0)
            out[name] = {"pressure": Fp, "viscous": Fv, "total": Fp + Fv}
        return out

    def moment(self, patch, center=(0.0, 0.0)):
        sl = self.patch_slices[patch]
        pb = self.boundary_p()[sl]
        dF = (pb - self.fs.p)[:, None] * self.mesh.Sf[self.ni:][sl]
        dF = dF + self._viscous_traction(sl)
        r = self.mesh.face_centers[self.ni:][sl] - np.asarray(center, float)
        return float(np.sum(r[:, 0] * dF[:, 1] - r[:, 1] * dF[:, 0]))

    def wall_shear(self, patch):
        """(centres de faces, τ_w signé selon la tangente t = (−n_y, n_x), y⁺)."""
        self._ensure_last()
        sl = self.patch_slices[patch]
        xf = self.mesh.face_centers[self.ni:][sl]
        n = self.nb_hat[sl]
        t = np.column_stack([-n[:, 1], n[:, 0]])
        tau = np.sum(self._viscous_traction(sl) * t, axis=1) / self.magSb[sl]
        W = self.primitive()
        c = self.Pb[sl]
        rho = W[0, c]
        nu = self.gas.mu_of(W[3, c] / (rho * self.gas.R)) / rho
        if not self.gas.viscous:
            return xf, tau, np.zeros(len(tau))
        yplus = self.distb[sl] * np.sqrt(np.abs(tau) / rho) / nu
        return xf, tau, yplus

    def wall_heat_flux(self, patch):
        """(T_paroi, flux de chaleur q entrant dans le fluide, W/m²)."""
        self._ensure_last()
        sl = self.patch_slices[patch]
        Tb = self.boundary_values(self.primitive())[4, sl]
        if "Fvb" not in self._last:
            return Tb, np.zeros(len(Tb))
        return Tb, self._last["Fvb"][3, sl].copy()

    def fields(self) -> dict:
        W = self.primitive()
        g = self.gas.gamma
        fs = self.fs
        r, u, v, p = W
        c = np.sqrt(g * p / r)
        gU = self.grad_U()
        q = 0.5 * fs.rho * max(fs.speed, 1e-300) ** 2
        return {"rho": r.copy(), "U": np.column_stack([u, v]), "U_mag": np.hypot(u, v),
                "p": p.copy(), "T": p / (r * self.gas.R), "Mach": np.hypot(u, v) / c,
                "Cp": (p - fs.p) / q, "entropy": (p / fs.p) / (r / fs.rho) ** g - 1.0,
                "vorticity": gU[:, 1, 0] - gU[:, 0, 1]}

    def mean_fields(self) -> dict:
        return {}

    def to_cpu(self):
        return self

    def totals(self) -> dict:
        """Masse, quantité de mouvement et énergie totales (intégrales sur le domaine)."""
        tot = self.Q @ self.V
        return {"mass": float(tot[0]), "momentum_x": float(tot[1]),
                "momentum_y": float(tot[2]), "energy": float(tot[3])}


def _take(a, idx):
    """a[:, idx] (np.take : 3 à 4 fois plus rapide que l'indexation avancée)."""
    return np.take(a, idx, axis=1)


def _angle(st: State) -> float:
    return float(np.degrees(np.arctan2(st.v, st.u))) if st.speed > 0 else 0.0


def _state_from_spec(spec: dict, default: State) -> State:
    """État imposé d'une condition (farfield, supersonic_inlet) : mach, pressure,
    temperature, density, angle (°) ou velocity ; valeurs manquantes : `default`."""
    keys = ("mach", "pressure", "temperature", "density", "velocity", "angle")
    if not any(k in spec for k in keys):
        return default
    gas = default.gas
    p, T, rho = spec.get("pressure"), spec.get("temperature"), spec.get("density")
    if sum(v is not None for v in (p, T, rho)) < 2:
        if p is None:
            p = default.p
        if sum(v is not None for v in (p, T, rho)) < 2:
            T = default.T if T is None else T
    return make_state(gas, spec.get("mach", default.mach), p, T,
                      rho if (p is None or T is None) else None,
                      spec.get("angle", _angle(default)), spec.get("velocity"))


def _value(v, x, y):
    """Constante ou formule en x, y (safe_expr)."""
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return np.full(len(x), float(v))
    from ..safe_expr import evaluate
    val = evaluate(v, {"x": x, "y": y})
    return np.broadcast_to(np.asarray(val, dtype=float), x.shape).copy()
