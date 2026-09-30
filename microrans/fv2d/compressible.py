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
        l1 = np.where(l1 < d, 0.5 * (l1 * l1 + d * d) / d, l1)
        l3 = np.where(l3 < d, 0.5 * (l3 * l3 + d * d) / d, l3)
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
        self._wref = np.array([freestream.rho, vref, vref, freestream.p])
        self._geometry()
        self._setup_bc(boundaries)
        self._grad_ops()
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
        self._frozen_lim = None
        self._last = {}
        # état initial
        C = mesh.cell_centers
        if initial is None:
            W = np.tile([freestream.rho, freestream.u, freestream.v, freestream.p],
                        (self.nc, 1))
        else:
            r, u, v, p = initial(C[:, 0], C[:, 1])
            W = np.column_stack([np.broadcast_to(np.asarray(a, float), (self.nc,))
                                 for a in (r, u, v, p)])
        if np.any(W[:, 0] <= 0) or np.any(W[:, 3] <= 0):
            raise ValueError("État initial : ρ et p doivent être > 0 partout.")
        self.Q = self.conservative(W)

    # ------------------------------------------------------------------ géométrie
    def _geometry(self):
        m, f = self.mesh, self.fvm
        self.nc, self.ni, self.nb = f.nc, f.ni, f.nb
        ni = self.ni
        self.P, self.N, self.Pb = m.owner[:ni], m.neighbour, m.owner[ni:]
        self.V = m.cell_volumes
        self.magS = m.magSf[:ni]
        self.n = m.nf[:ni]
        self.magSb = m.magSf[ni:]
        self.nb_hat = m.nf[ni:]
        self.rP = np.asarray(f.rP)
        self.rN = np.asarray(f.rN)
        self.rb = m.face_centers[ni:] - m.cell_centers[self.Pb]
        d = m.d_PN
        self.dist = np.linalg.norm(d, axis=1)
        self.e = d / self.dist[:, None]
        self.distb = np.linalg.norm(self.rb, axis=1)
        self.eb = self.rb / self.distb[:, None]
        self.w = m.weights
        S2 = np.bincount(self.P, self.magS ** 2, self.nc) + np.bincount(
            self.N, self.magS ** 2, self.nc)
        if self.nb:
            S2 = S2 + np.bincount(self.Pb, self.magSb ** 2, self.nc)
        self.sumS2 = S2
        self.h = np.sqrt(self.V)

    def _grad_ops(self):
        """Opérateur de Green-Gauss (identique à FVM.grad) en matrices creuses
        Gx, Gy : (nc, nc + nb), appliquées à [valeurs de cellules ; valeurs frontières]."""
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
            ops.append(sp.csr_matrix((vals, (rows, cols)), shape=(nc, nc + nb)))
        self.Gx, self.Gy = ops

    def gradient(self, W, Wb):
        """Gradients de Green-Gauss de plusieurs variables : W (nc, k), Wb (nb, k) ->
        (gx, gy) chacun (nc, k)."""
        ext = np.vstack([W, Wb])
        return self.Gx @ ext, self.Gy @ ext

    # ------------------------------------------------------------------ conversions
    def conservative(self, W):
        g = self.gas.gamma
        r, u, v, p = W[:, 0], W[:, 1], W[:, 2], W[:, 3]
        return np.column_stack([r, r * u, r * v, p / (g - 1.0) + 0.5 * r * (u * u + v * v)])

    def primitive(self, Q=None):
        Q = self.Q if Q is None else Q
        g = self.gas.gamma
        r = Q[:, 0]
        u, v = Q[:, 1] / r, Q[:, 2] / r
        p = (g - 1.0) * (Q[:, 3] - 0.5 * r * (u * u + v * v))
        return np.column_stack([r, u, v, p])

    # ------------------------------------------------------------------ conditions
    def _setup_bc(self, boundaries: dict):
        m, gas, fs = self.mesh, self.gas, self.fs
        nb = self.nb
        self.kind = np.full(nb, -1, dtype=int)
        self.bc_state = np.tile([fs.rho, fs.u, fs.v, fs.p], (nb, 1))     # état imposé
        self.bc_p = np.full(nb, fs.p)                                    # p (sortie)
        self.bc_p0 = np.full(nb, fs.p0)
        self.bc_T0 = np.full(nb, fs.T0)
        sp_ = max(fs.speed, 1e-300)
        self.bc_dir = np.tile([fs.u / sp_, fs.v / sp_] if fs.speed > 0 else [1.0, 0.0],
                              (nb, 1))
        self.bc_Tw = np.full(nb, np.nan)                                 # nan : adiabatique
        self.bc_Uw = np.zeros((nb, 2))
        self.patch_slices = {}
        self.bc_types = {}
        periodic = {q for pair in m.periodic_pairs for q in pair[:2]}
        missing = [p.name for p in m.patches if p.name not in boundaries
                   and p.type != "empty"]
        if missing:
            raise ValueError(f"Conditions aux limites manquantes pour : {missing} "
                             f"(types : {', '.join(COMP_BC_TYPES)})")
        for name in boundaries:
            if name in periodic:
                continue
            if name not in {p.name for p in m.patches}:
                raise ValueError(f"[boundary.{name}] : patch absent du maillage "
                                 f"({[p.name for p in m.patches]}).")
        C = m.face_centers[self.ni:]
        for patch in m.patches:
            sl = slice(patch.start - self.ni, patch.start - self.ni + patch.size)
            self.patch_slices[patch.name] = sl
            spec = dict(boundaries.get(patch.name, {"type": "slip_wall"}))
            if patch.type == "empty":
                spec = {"type": "slip_wall"}
            kind = canonical_bc(spec.get("type", "wall"))
            if kind == "wall" and not gas.viscous:
                kind = "slip_wall"                  # Euler : paroi = glissement
            self.bc_types[patch.name] = kind
            code = {"slip_wall": _SLIP, "symmetry": _SLIP, "wall": _NOSLIP,
                    "farfield": _FAR, "supersonic_inlet": _SUPIN, "inlet": _INLET,
                    "outlet": _OUTLET, "supersonic_outlet": _EXTRAP}[kind]
            self.kind[sl] = code
            x, y = C[sl, 0], C[sl, 1]
            if kind in ("farfield", "supersonic_inlet") and any(
                    k in spec for k in ("mach", "pressure", "temperature", "density",
                                        "velocity", "angle")):
                st = make_state(gas, spec.get("mach", fs.mach),
                                spec.get("pressure", None if "density" in spec and
                                         "temperature" in spec else fs.p),
                                spec.get("temperature", None if "density" in spec
                                         else fs.T),
                                spec.get("density"), spec.get("angle", _angle(fs)),
                                spec.get("velocity"))
                self.bc_state[sl] = [st.rho, st.u, st.v, st.p]
            if kind == "inlet":
                self.bc_p0[sl] = _value(spec.get("p0", spec.get("total_pressure", fs.p0)),
                                        x, y)
                self.bc_T0[sl] = _value(spec.get("T0", spec.get(
                    "total_temperature", fs.T0)), x, y)
                if "direction" in spec or "angle" in spec:
                    if "direction" in spec:
                        d = np.asarray(spec["direction"], float)
                    else:
                        a = np.radians(float(spec["angle"]))
                        d = np.array([np.cos(a), np.sin(a)])
                    self.bc_dir[sl] = d / np.linalg.norm(d)
            if kind == "outlet":
                self.bc_p[sl] = _value(spec.get("p", spec.get("pressure", fs.p)), x, y)
            if kind == "wall":
                if "T" in spec or "temperature" in spec:
                    self.bc_Tw[sl] = _value(spec.get("T", spec.get("temperature")), x, y)
                if "U" in spec:
                    U = spec["U"]
                    self.bc_Uw[sl, 0] = _value(U[0], x, y)
                    self.bc_Uw[sl, 1] = _value(U[1], x, y)
        if np.any(self.kind < 0):
            raise ValueError("Faces frontières sans condition aux limites.")
        self.is_wall = np.isin(self.kind, (_SLIP, _NOSLIP))
        self.is_noslip = self.kind == _NOSLIP
        self.is_iso = self.is_noslip & np.isfinite(self.bc_Tw)

    def ghost(self, W):
        """État fantôme (nb, 4) des faces frontières à partir de l'état intérieur W (nb, 4)
        (reconstruit à la face ou valeur de cellule)."""
        g = self.gas.gamma
        gm1 = g - 1.0
        k = self.kind
        n = self.nb_hat
        nx, ny = n[:, 0], n[:, 1]
        r, u, v, p = W[:, 0], W[:, 1], W[:, 2], W[:, 3]
        qn = u * nx + v * ny
        c = np.sqrt(g * p / r)
        G = W.copy()
        # parois glissantes / symétrie : miroir
        s = k == _SLIP
        if s.any():
            G[s, 1] = u[s] - 2.0 * qn[s] * nx[s]
            G[s, 2] = v[s] - 2.0 * qn[s] * ny[s]
        s = k == _NOSLIP
        if s.any():
            G[s, 1] = 2.0 * self.bc_Uw[s, 0] - u[s]
            G[s, 2] = 2.0 * self.bc_Uw[s, 1] - v[s]
        s = k == _SUPIN
        if s.any():
            G[s] = self.bc_state[s]
        s = k == _FAR
        if s.any():
            G[s] = self._farfield(W[s], self.bc_state[s], nx[s], ny[s])
        s = k == _INLET
        if s.any():
            Rp = qn[s] + 2.0 * c[s] / gm1
            d = self.bc_dir[s]
            cs = np.maximum(-(d[:, 0] * nx[s] + d[:, 1] * ny[s]), 0.05)
            H0 = self.gas.cp * self.bc_T0[s]
            a = 0.25 * gm1 * cs * cs + 0.5
            b = 0.5 * gm1 * Rp * cs
            cc = 0.25 * gm1 * Rp * Rp - H0
            Vm = np.maximum((-b + np.sqrt(np.maximum(b * b - 4 * a * cc, 0.0))) / (2 * a), 0.0)
            Tb = np.maximum(self.bc_T0[s] - 0.5 * Vm * Vm / self.gas.cp, 1e-3 * self.bc_T0[s])
            pb = self.bc_p0[s] * (Tb / self.bc_T0[s]) ** (g / gm1)
            G[s, 0] = pb / (self.gas.R * Tb)
            G[s, 1] = Vm * d[:, 0]
            G[s, 2] = Vm * d[:, 1]
            G[s, 3] = pb
        s = k == _OUTLET
        if s.any():
            sub = qn[s] < c[s]                     # sortie supersonique : extrapolation
            pb = np.where(sub, self.bc_p[s], p[s])
            rb = r[s] * (pb / p[s]) ** (1.0 / g)
            cb = np.sqrt(g * pb / rb)
            qb = qn[s] + 2.0 * (c[s] - cb) / gm1
            G[s, 0] = rb
            G[s, 1] = u[s] + (qb - qn[s]) * nx[s]
            G[s, 2] = v[s] + (qb - qn[s]) * ny[s]
            G[s, 3] = pb
        return G

    def _farfield(self, W, Winf, nx, ny):
        """Invariants de Riemann normaux (Blazek § 8.4 ; Jameson & Baker 1983)."""
        g = self.gas.gamma
        gm1 = g - 1.0
        r, u, v, p = W.T
        ri, ui, vi, pi = Winf.T
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
        G = np.column_stack([rb, ut + qb * nx, vt + qb * ny, pb])
        # supersonique : tout de l'amont (entrée) ou de l'intérieur (sortie)
        sup_in = qi <= -ci
        sup_out = qn >= c
        G[sup_in] = Winf[sup_in]
        G[sup_out] = W[sup_out]
        return G

    def boundary_values(self, W):
        """Valeurs aux faces frontières (nb, 5) de (ρ, u, v, p, T) pour les gradients :
        moyenne cellule / état fantôme ; parois adhérentes : vitesse de paroi, T_w."""
        G = self.ghost(W[self.Pb])
        Wb = 0.5 * (W[self.Pb] + G)
        s = self.is_noslip
        Wb[s, 1:3] = self.bc_Uw[s]
        Wb[s, 0] = W[self.Pb[s], 0]
        Wb[s, 3] = W[self.Pb[s], 3]
        T = Wb[:, 3] / (Wb[:, 0] * self.gas.R)
        iso = self.is_iso
        T[iso] = self.bc_Tw[iso]
        return np.column_stack([Wb, T])

    # ------------------------------------------------------------------ limiteurs
    def _limiter(self, W, Wb, gx, gy):
        """Coefficient ψ (nc, k) du limiteur choisi."""
        lim = self.settings.limiter
        if lim == "none":
            return np.ones_like(W)
        if getattr(self, "_stencil", None) is None:
            # voisins de chaque cellule (même pochoir que fvm.limit_grad), rangés
            # (m, nc) : les réductions sur le pochoir portent sur des tranches contiguës
            # (la dernière colonne du pochoir de fvm est toujours la cellule elle-même : ôtée)
            idx, R = self.fvm._limiter_stencil()
            idx, R = np.asarray(idx)[:, :-1], np.asarray(R)[:, :-1]
            self._stencil = (np.ascontiguousarray(idx.T),
                             np.ascontiguousarray(R[:, :, 0].T)[:, :, None],
                             np.ascontiguousarray(R[:, :, 1].T)[:, :, None])
        idx, Rx, Ry = self._stencil
        ext = np.vstack([W, Wb])[idx]                            # (m, nc, k)
        dmax = np.maximum(ext.max(axis=0), W) - W
        dmin = np.minimum(ext.min(axis=0), W) - W
        d = gx[None] * Rx + gy[None] * Ry                        # variation cellule → face
        if lim == "barth_jespersen":
            pos, neg = d > 0.0, d < 0.0
            r = np.where(pos, dmax[None] / np.where(pos, d, 1.0),
                         np.where(neg, dmin[None] / np.where(neg, d, -1.0), 1.0))
            return np.minimum(r.min(axis=0), 1.0)
        # Venkatakrishnan, seuil de Wang : ε² = (K (max − min global))²
        # (plancher : 1e-6 × échelle de référence, pour une variable uniforme)
        rng = np.maximum(W.max(axis=0) - W.min(axis=0), 1e-6 * self._wref)
        eps2 = (self.settings.venkat_k * rng) ** 2
        D1 = np.where(d > 0.0, dmax[None], dmin[None])
        D1sq = D1 * D1 + eps2
        dD1 = d * D1
        psi = (D1sq + 2.0 * dD1) / (D1sq + 2.0 * d * d + dD1)
        return np.minimum(psi.min(axis=0), 1.0)

    # ------------------------------------------------------------------ résidu
    def residual(self, Q, order=None, it=None, store=False):
        """Résidu R(Q) (nc, 4) = Σ_f (F_c − F_v)·S_f (sortant) ; dQ/dt = −R/V."""
        s = self.settings
        g = self.gas.gamma
        order = s.order if order is None else order
        P, N, Pb = self.P, self.N, self.Pb
        W = self.primitive(Q)
        if np.any(W[:, 0] <= 0.0) or np.any(W[:, 3] <= 0.0) or not np.all(np.isfinite(W)):
            bad = int(np.sum((W[:, 0] <= 0) | (W[:, 3] <= 0) | ~np.isfinite(W).all(axis=1)))
            raise FloatingPointError(f"État non physique (ρ ≤ 0 ou p ≤ 0) dans {bad} "
                                     "cellule(s) : réduire le CFL, démarrer à l'ordre 1 "
                                     "(first_order_iter) ou changer de flux (hllc).")
        viscous = self.gas.viscous
        Wb5 = self.boundary_values(W)
        need_grad = order == 2 or viscous
        if need_grad:
            T = W[:, 3] / (W[:, 0] * self.gas.R)
            W5 = np.column_stack([W, T])
            gx, gy = self.gradient(W5, Wb5)
        if order == 2:
            gx4, gy4 = gx[:, :4], gy[:, :4]
            if s.limiter_freeze and it is not None and it > s.limiter_freeze \
                    and self._frozen_lim is not None:
                psi = self._frozen_lim
            else:
                psi = self._limiter(W, Wb5[:, :4], gx4, gy4)
                if s.limiter_freeze and it is not None:
                    self._frozen_lim = psi
            lx, ly = gx4 * psi, gy4 * psi
            WL = W[P] + lx[P] * self.rP[:, :1] + ly[P] * self.rP[:, 1:]
            WR = W[N] + lx[N] * self.rN[:, :1] + ly[N] * self.rN[:, 1:]
            Wf = W[Pb] + lx[Pb] * self.rb[:, :1] + ly[Pb] * self.rb[:, 1:]
            bad = (WL[:, 0] <= 0) | (WL[:, 3] <= 0) | (WR[:, 0] <= 0) | (WR[:, 3] <= 0)
            if bad.any():
                WL[bad], WR[bad] = W[P[bad]], W[N[bad]]
            badb = (Wf[:, 0] <= 0) | (Wf[:, 3] <= 0)
            if badb.any():
                Wf[badb] = W[Pb[badb]]
        else:
            WL, WR, Wf = W[P], W[N], W[Pb]
        flux = _FLUX_FUNCS[s.flux]
        n = self.n
        Fi = np.column_stack(flux(tuple(WL.T), tuple(WR.T), n[:, 0], n[:, 1], g,
                                  s.entropy_fix))
        G = self.ghost(Wf)
        nb = self.nb_hat
        Fb = np.column_stack(flux(tuple(Wf.T), tuple(G.T), nb[:, 0], nb[:, 1], g,
                                  s.entropy_fix))
        if viscous:
            Fvi, Fvb = self._viscous_flux(W5, Wb5, gx, gy)
            Fi = Fi - Fvi
            Fb = Fb - Fvb
        Fi = Fi * self.magS[:, None]
        Fb = Fb * self.magSb[:, None]
        R = np.empty((self.nc, 4))
        for k in range(4):
            R[:, k] = (np.bincount(P, Fi[:, k], self.nc) - np.bincount(N, Fi[:, k], self.nc)
                       + np.bincount(Pb, Fb[:, k], self.nc))
        if store:
            self._last = {"Fb": Fb, "W": W}
            if viscous:
                self._last["Fvb"] = Fvb
        return R

    def _viscous_flux(self, W5, Wb5, gx, gy):
        """Flux visqueux (par unité de surface) aux faces internes et frontières."""
        gas = self.gas
        P, N, Pb = self.P, self.N, self.Pb
        w = self.w[:, None]
        # variables u, v, T (colonnes 1, 2, 4)
        cols = [1, 2, 4]
        phi = W5[:, cols]
        gax = w * gx[P][:, cols] + (1 - w) * gx[N][:, cols]
        gay = w * gy[P][:, cols] + (1 - w) * gy[N][:, cols]
        ex, ey = self.e[:, :1], self.e[:, 1:]
        corr = (phi[N] - phi[P]) / self.dist[:, None] - (gax * ex + gay * ey)
        gfx, gfy = gax + corr * ex, gay + corr * ey
        uf = w * phi[P] + (1 - w) * phi[N]
        Fi = self._stress_flux(uf, gfx, gfy, self.n)
        # frontières : parois adhérentes (correction selon d_b) ; autres : gradient de cellule
        phib = Wb5[:, cols]
        gbx, gby = gx[Pb][:, cols], gy[Pb][:, cols]
        s = self.is_noslip
        if s.any():
            ebx, eby = self.eb[s, :1], self.eb[s, 1:]
            cb = (phib[s] - phi[Pb[s]]) / self.distb[s, None] - (gbx[s] * ebx + gby[s] * eby)
            gbx[s] = gbx[s] + cb * ebx
            gby[s] = gby[s] + cb * eby
            # adiabatique : flux de chaleur nul
            ad = ~np.isfinite(self.bc_Tw[s])
            if ad.any():
                gxT, gyT = gbx[s, 2], gby[s, 2]
                nn = self.nb_hat[s]
                qn = gxT * nn[:, 0] + gyT * nn[:, 1]
                gxT = np.where(ad, gxT - qn * nn[:, 0], gxT)
                gyT = np.where(ad, gyT - qn * nn[:, 1], gyT)
                gbx[s, 2], gby[s, 2] = gxT, gyT
        Fb = self._stress_flux(phib, gbx, gby, self.nb_hat)
        slip = self.kind == _SLIP
        Fb[slip] = 0.0
        return Fi, Fb

    def _stress_flux(self, uvT, gx, gy, n):
        gas = self.gas
        u, v, T = uvT[:, 0], uvT[:, 1], uvT[:, 2]
        ux, vx, Tx = gx[:, 0], gx[:, 1], gx[:, 2]
        uy, vy, Ty = gy[:, 0], gy[:, 1], gy[:, 2]
        mu = gas.mu_of(np.maximum(T, 1e-3))
        kcond = mu * gas.cp / gas.Pr
        div = ux + vy
        txx = mu * (2.0 * ux - 2.0 / 3.0 * div)
        tyy = mu * (2.0 * vy - 2.0 / 3.0 * div)
        txy = mu * (uy + vx)
        nx, ny = n[:, 0], n[:, 1]
        fx = txx * nx + txy * ny
        fy = txy * nx + tyy * ny
        fe = u * fx + v * fy + kcond * (Tx * nx + Ty * ny)
        return np.column_stack([np.zeros_like(fx), fx, fy, fe])

    # ------------------------------------------------------------------ pas de temps
    def local_dt(self, Q, cfl):
        """Δt_i = CFL V_i / (Λc_i + C_v Λv_i) (rayons spectraux, Blazek § 6.1.4)."""
        g = self.gas.gamma
        W = self.primitive(Q)
        c = np.sqrt(g * W[:, 3] / W[:, 0])
        P, N, Pb = self.P, self.N, self.Pb
        n = self.n
        un = 0.5 * np.abs((W[P, 1] + W[N, 1]) * n[:, 0] + (W[P, 2] + W[N, 2]) * n[:, 1])
        lam = (un + 0.5 * (c[P] + c[N])) * self.magS
        L = np.bincount(P, lam, self.nc) + np.bincount(N, lam, self.nc)
        if self.nb:
            nb = self.nb_hat
            lb = (np.abs(W[Pb, 1] * nb[:, 0] + W[Pb, 2] * nb[:, 1]) + c[Pb]) * self.magSb
            L = L + np.bincount(Pb, lb, self.nc)
        if self.gas.viscous:
            T = W[:, 3] / (W[:, 0] * self.gas.R)
            mu = self.gas.mu_of(T)
            Lv = max(4.0 / 3.0, g / self.gas.Pr) * mu / W[:, 0] * self.sumS2 / self.V
            L = L + self.settings.viscous_factor * Lv
        return cfl * self.V / L

    # ------------------------------------------------------------------ intégration
    def _ssp_rk3(self, dt, order, it=None):
        Q0 = self.Q
        f = (dt / self.V)[:, None] if np.ndim(dt) else dt / self.V[:, None]
        Q1 = Q0 - f * self.residual(Q0, order, it)
        Q2 = 0.75 * Q0 + 0.25 * (Q1 - f * self.residual(Q1, order, it))
        R2 = self.residual(Q2, order, it)
        self.Q = Q0 / 3.0 + 2.0 / 3.0 * (Q2 - f * R2)
        return None

    def _rk5(self, dt, order, it=None):
        Q0 = self.Q
        f = (dt / self.V)[:, None]
        Q = Q0
        R0 = None
        for a in _RK5:
            R = self.residual(Q, order, it)
            if R0 is None:
                R0 = R
            Q = Q0 - a * f * R
        self.Q = Q
        return R0

    def residual_norms(self, R) -> np.ndarray:
        """Normes RMS des résidus par unité de volume (4 équations)."""
        return np.sqrt(np.mean((R / self.V[:, None]) ** 2, axis=0))

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
        if s.steady_scheme == "implicit":
            from .compressible_implicit import ImplicitStepper
            stepper = ImplicitStepper(self)
            cfl = min(s.cfl, s.cfl_max)
        forces_hist = []
        it = 0
        for it in range(1, max_iter + 1):
            gi = self.iterations_total + it
            order = 1 if gi <= s.first_order_iter else s.order
            if s.steady_scheme == "implicit":
                R = stepper.step(cfl, order, gi)
                cfl = min(cfl * s.cfl_growth, s.cfl_max)
            else:
                dt = self.local_dt(self.Q, cfl)
                if s.steady_scheme == "rk5":
                    R = self._rk5(dt, order, gi)
                else:
                    R = self.residual(self.Q, order, gi)
                    self._ssp_rk3_from(R, dt, order, gi)
            norms = self.residual_norms(R)
            rel = self._relative(norms, gi)
            rec = {"iteration": gi, "rho": rel[0], "rhoU": rel[1], "rhoV": rel[2],
                   "rhoE": rel[3]}
            if not all(np.isfinite(v) for v in rec.values()) or \
                    not np.all(np.isfinite(self.Q)):
                raise FloatingPointError(f"Divergence à l'itération {gi}.")
            self.history.append(rec)
            if monitor is not None and (it % 10 == 0 or it == 1):
                mrec = {"iteration": gi, **monitor(self)}
                self.monitor.append(mrec)
            if verbose and (it % log_every == 0 or it == 1):
                print("  it %6d  " % gi + "  ".join(f"{k}={v:.2e}" for k, v in rec.items()
                                                     if k != "iteration")
                      + (f"  CFL={cfl:.3g}" if s.steady_scheme == "implicit" else ""))
            if s.monitor_tol and it % 10 == 0:
                forces_hist.append(np.concatenate([f["total"] for f in
                                                   self.forces().values()] or [np.zeros(1)]))
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

    def _ssp_rk3_from(self, R0, dt, order, it):
        """SSP-RK3 dont le 1er résidu (déjà calculé, sert au suivi) est R0."""
        Q0 = self.Q
        f = (dt / self.V)[:, None]
        Q1 = Q0 - f * R0
        Q2 = 0.75 * Q0 + 0.25 * (Q1 - f * self.residual(Q1, order, it))
        self.Q = Q0 / 3.0 + 2.0 / 3.0 * (Q2 - f * self.residual(Q2, order, it))

    def _relative(self, norms, gi):
        """Résidus normalisés par le maximum observé sur les 10 premières itérations ;
        quantité de mouvement : échelle commune aux deux composantes."""
        if self._res0 is None:
            self._res0 = np.zeros(4)
            self._res0_count = 0
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
            step = dt if dt else float(np.min(self.local_dt(self.Q, cfl)))
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
        return self.primitive()

    @property
    def U(self):
        W = self.primitive()
        return W[:, 1:3].copy()

    @property
    def p(self):
        return self.primitive()[:, 3]

    @property
    def rho(self):
        return self.Q[:, 0]

    @property
    def T(self):
        W = self.primitive()
        return W[:, 3] / (W[:, 0] * self.gas.R)

    @property
    def mach(self):
        W = self.primitive()
        return np.hypot(W[:, 1], W[:, 2]) / np.sqrt(self.gas.gamma * W[:, 3] / W[:, 0])

    def boundary_p(self, p=None):
        """Pression aux faces frontières (nb,) : valeur de paroi issue du flux de Riemann
        sur les parois, moyenne cellule / fantôme ailleurs."""
        W = self.primitive()
        Wb = self.boundary_values(W)
        pb = Wb[:, 3].copy()
        last = self._last.get("Fb")
        if last is not None and self.is_wall.any():
            nb = self.nb_hat
            Fc = last.copy()
            if "Fvb" in self._last:
                Fc = Fc + self._last["Fvb"] * self.magSb[:, None]
            pw = (Fc[:, 1] * nb[:, 0] + Fc[:, 2] * nb[:, 1]) / self.magSb
            pb[self.is_wall] = pw[self.is_wall]
        return pb

    def grad_U(self, U=None):
        W = self.primitive()
        Wb = self.boundary_values(W)
        gx, gy = self.gradient(W[:, 1:3], Wb[:, 1:3])
        out = np.empty((self.nc, 2, 2))
        out[:, :, 0], out[:, :, 1] = gx, gy
        return out

    def forces(self, patches=None):
        """Efforts sur des patches (N/m d'envergure) : pression (p − p∞) et frottement."""
        names = patches or [p.name for p in self.mesh.patches
                            if self.bc_types.get(p.name) in ("wall", "slip_wall")]
        if not self._last:
            self.residual(self.Q, self.settings.order, None, store=True)
        pb = self.boundary_p()
        Sb = self.mesh.Sf[self.ni:]
        out = {}
        for name in names:
            sl = self.patch_slices[name]
            Fp = np.sum((pb[sl] - self.fs.p)[:, None] * Sb[sl], axis=0)
            if "Fvb" in self._last:
                Fv = -np.sum(self._last["Fvb"][sl, 1:3] * self.magSb[sl, None], axis=0)
            else:
                Fv = np.zeros(2)
            out[name] = {"pressure": Fp, "viscous": Fv, "total": Fp + Fv}
        return out

    def moment(self, patch, center=(0.0, 0.0)):
        sl = self.patch_slices[patch]
        pb = self.boundary_p()[sl]
        Sb = self.mesh.Sf[self.ni:][sl]
        dF = (pb - self.fs.p)[:, None] * Sb
        if "Fvb" in self._last:
            dF = dF - self._last["Fvb"][sl, 1:3] * self.magSb[sl, None]
        r = self.mesh.face_centers[self.ni:][sl] - np.asarray(center, float)
        return float(np.sum(r[:, 0] * dF[:, 1] - r[:, 1] * dF[:, 0]))

    def wall_shear(self, patch):
        """(centres de faces, τ_w signé selon la tangente t = (−n_y, n_x), y⁺)."""
        sl = self.patch_slices[patch]
        n = self.nb_hat[sl]
        xf = self.mesh.face_centers[self.ni:][sl]
        if "Fvb" not in self._last:
            return xf, np.zeros(len(n)), np.zeros(len(n))
        t = np.column_stack([-n[:, 1], n[:, 0]])
        tau = -np.sum(self._last["Fvb"][sl, 1:3] * t, axis=1)
        W = self.primitive()
        c = self.Pb[sl]
        rho = W[c, 0]
        T = W[c, 3] / (rho * self.gas.R)
        nu = self.gas.mu_of(T) / rho
        yplus = self.distb[sl] * np.sqrt(np.abs(tau) / rho) / nu
        return xf, tau, yplus

    def wall_heat_flux(self, patch):
        """(T_paroi, flux de chaleur q entrant dans le fluide, W/m²)."""
        sl = self.patch_slices[patch]
        Tb = self.boundary_values(self.primitive())[sl, 4]
        if "Fvb" not in self._last:
            return Tb, np.zeros(len(Tb))
        return Tb, self._last["Fvb"][sl, 3].copy()

    def fields(self) -> dict:
        W = self.primitive()
        g = self.gas.gamma
        fs = self.fs
        U = W[:, 1:3]
        c = np.sqrt(g * W[:, 3] / W[:, 0])
        gU = self.grad_U()
        q = 0.5 * fs.rho * max(fs.speed, 1e-300) ** 2
        return {"rho": W[:, 0], "U": U, "U_mag": np.hypot(U[:, 0], U[:, 1]), "p": W[:, 3],
                "T": W[:, 3] / (W[:, 0] * self.gas.R), "Mach": np.hypot(U[:, 0], U[:, 1]) / c,
                "Cp": (W[:, 3] - fs.p) / q,
                "entropy": (W[:, 3] / fs.p) / (W[:, 0] / fs.rho) ** g - 1.0,
                "vorticity": gU[:, 1, 0] - gU[:, 0, 1]}

    def mean_fields(self) -> dict:
        return {}

    def to_cpu(self):
        return self

    def totals(self) -> dict:
        """Masse, quantité de mouvement et énergie totales (intégrales sur le domaine)."""
        tot = np.sum(self.Q * self.V[:, None], axis=0)
        return {"mass": float(tot[0]), "momentum_x": float(tot[1]),
                "momentum_y": float(tot[2]), "energy": float(tot[3])}


def _angle(st: State) -> float:
    return float(np.degrees(np.arctan2(st.v, st.u))) if st.speed > 0 else 0.0


def _value(v, x, y):
    """Constante ou formule en x, y (safe_expr)."""
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return np.full(len(x), float(v))
    from ..safe_expr import evaluate
    val = evaluate(v, {"x": x, "y": y})
    return np.broadcast_to(np.asarray(val, dtype=float), x.shape).copy()
