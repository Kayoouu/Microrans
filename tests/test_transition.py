"""Transition laminaire-turbulent : modèle γ à une équation de Menter et al. (2015) couplé
au k-ω SST (`sst_gamma`). Plaque plane de type T3A (ERCOFTAC) sur maillage grossier,
registre des modèles, conditions aux limites de γ, chemin GPU (faux GPU)."""
import warnings

import numpy as np
import pytest

from fake_device import register_fake_gpu
from microrans.cases import run_rans_channel
from microrans.fv2d import Settings, Solver2D
from microrans.mesh2d import Circle, Rectangle, flat_plate_mesh
from microrans.mesh2d.unstructured import triangulate
from microrans.models import MODELS, TURBULENT_MODELS, canonical_name, get_model
from microrans.models.transition_gamma import freestream_decay

U_T3A, NU_T3A = 5.4, 1.5e-5
# entrée à 0.04 m en amont du bord d'attaque : Tu = 3.7 %, ν_t/ν = 12.4 → Tu ≈ 3.35 % au
# bord d'attaque, décroissance ajustée sur les mesures T3A (voir docs/notes_transition.md)
T3A_INFLOW = {"intensity": 0.037, "viscosity_ratio": 12.4}
PLATE_BCS = {"inlet": {"type": "inlet", "U": [U_T3A, 0.0]},
             "outlet": {"type": "outlet", "p": 0.0}, "top": {"type": "outlet", "p": 0.0},
             "symmetry": {"type": "symmetry"}, "plate": {"type": "wall"}}


def blasius(rex):
    return 0.664 / np.sqrt(rex)


def turbulent(rex):
    return 0.0576 * rex ** -0.2


@pytest.fixture(scope="module")
def t3a_coarse():
    """T3A grossier : 4 100 cellules, 1re maille 5e-5 m (y⁺ ≤ 1.4), convection des
    équations de turbulence au 2e ordre limité (comme l'exemple), ~7 s."""
    m = flat_plate_mesh(length=1.6, upstream=0.04, height=0.3, nx_up=6, nx_plate=70, ny=50,
                        first_height=5e-5, le_fraction=3e-3 / 1.6)
    s = Solver2D(m, NU_T3A, PLATE_BCS, model="sst_gamma", initial_U=(U_T3A, 0.0),
                 turbulence_inflow=T3A_INFLOW, reference_velocity=U_T3A,
                 settings=Settings(convection_turb="linearUpwindLimited"))
    assert s.run_steady(max_iter=1500, tol=1e-6)
    x, tau, yp = s.wall_shear("plate")
    o = np.argsort(x[:, 0])
    x = x[o, 0]
    return s, U_T3A * x / NU_T3A, 2.0 * tau[o] / U_T3A ** 2, yp[o]


def test_registry_and_freestream():
    assert canonical_name("sst_gamma") == "sst_gamma"
    assert canonical_name("transition") == "sst_gamma"
    assert "sst_gamma" in MODELS and "sst_gamma" not in TURBULENT_MODELS
    m = flat_plate_mesh(1.0, 0.1, 0.5, 4, 8, 8)
    model = get_model("sst_gamma", m, 1e-5)
    assert model.variables == ("k", "omega", "gamma")
    fs = model.freestream_values(10.0, 0.01, 10.0)
    assert fs["gamma"] == 1.0 and fs["k"] == pytest.approx(1.5 * 0.1 ** 2)
    # normale pariétale : +y au-dessus de la plaque, radiale devant le bord d'attaque
    n = m.wall_normal()
    C = m.cell_centers
    assert np.allclose(n[C[:, 0] > 0.0], [0.0, 1.0])
    up = C[:, 0] < 0.0
    assert np.allclose(n[up], C[up] / np.linalg.norm(C[up], axis=1)[:, None])


def test_gamma_boundary_conditions():
    """γ : valeur amont 1 en entrée, gradient nul à la paroi (pas de valeur imposée)."""
    m = flat_plate_mesh(1.0, 0.1, 0.5, 4, 8, 8)
    s = Solver2D(m, 1e-5, {**PLATE_BCS, "inlet": {"type": "inlet", "U": [1.0, 0.0]}},
                 model="sst_gamma")
    a, b, g, d = s.scalar_bc("gamma")
    wall = s.is_wall
    assert np.all(a[wall] == 1.0) and np.all(g[wall] == 0.0) and np.all(d[wall] == 0.0)
    inlet = s.patch_slices["inlet"]
    assert np.allclose(b[inlet], 1.0) and np.all(g[inlet] < 0.0)
    # k reste nul à la paroi (Dirichlet)
    a, b, g, d = s.scalar_bc("k")
    assert np.all(a[wall] == 0.0) and np.all(g[wall] < 0.0)
    with pytest.raises(ValueError, match="transition"):
        Solver2D(m, 1e-5, PLATE_BCS, model="sst_gamma",
                 settings=Settings(wall_treatment="wall_function"))


def test_t3a_laminar_region_near_blasius(t3a_coarse):
    """Avant la transition (Re_x < 1.2e5), C_f reste proche de Blasius mais au-dessus :
    +9 à +21 % mesurés, dus à la viscosité turbulente de l'écoulement libre (ν_t/ν ≈ 12)
    qui pénètre la couche limite laminaire (même maillage en laminaire : +1.5 à +2.5 % ;
    expérience T3A : 0 à +11 % sur cette plage). Moins de 20 % de l'écart
    laminaire → turbulent."""
    s, rex, cf, yp = t3a_coarse
    sel = (rex > 1.5e4) & (rex < 1.2e5)
    lam, turb = blasius(rex[sel]), turbulent(rex[sel])
    ratio = cf[sel] / lam
    assert np.all(ratio > 1.0) and np.all(ratio < 1.25)
    assert np.all((cf[sel] - lam) / (turb - lam) < 0.2)
    g = s.state["gamma"]
    assert g.min() >= 0.0 and g.max() <= 1.0


def test_t3a_turbulent_downstream(t3a_coarse):
    """Loin en aval (Re_x > 4.5e5), C_f à ±5 % de la corrélation turbulente 0.0576 Re_x^-0.2
    (mesuré : −2.4 à −3 % ; expérience T3A à x = 1.495 m : 0.00408, corrélation 0.00411)."""
    s, rex, cf, yp = t3a_coarse
    sel = rex > 4.5e5
    assert np.allclose(cf[sel] / turbulent(rex[sel]), 1.0, atol=0.05)
    assert yp.max() < 1.5


def test_t3a_transition_location(t3a_coarse):
    """Début de transition (minimum de C_f) et point à mi-chemin laminaire → turbulent :
    expérience T3A (Savill 1993) : minimum vers Re_x ≈ 1.4e5, mi-transition Re_x ≈ 2.3e5.
    Maillage grossier mesuré : 1.35e5 et 1.86e5 (90 % à 2.67e5) ; bandes larges car le
    maillage et le schéma déplacent la transition de ±10 %."""
    s, rex, cf, yp = t3a_coarse
    i0 = np.argmin(np.where(rex > 2e4, cf, 1.0))
    assert 1.0e5 < rex[i0] < 2.2e5
    g = (cf - blasius(rex)) / (turbulent(rex) - blasius(rex))
    j = i0 + np.argmax(g[i0:] > 0.5)
    assert 1.5e5 < rex[j] < 3.0e5
    k = i0 + np.argmax(g[i0:] > 0.9)
    assert rex[k] < 3.5e5


def test_freestream_decay_matches_analytic(t3a_coarse):
    """Tu hors couche limite = décroissance analytique du SST (β₂, β*) utilisée pour régler
    la turbulence amont : écart < 0.02 point de Tu (mesuré : ~0.003 sur le maillage moyen)."""
    s = t3a_coarse[0]
    C = s.mesh.cell_centers
    ycol = np.unique(C[:, 1])
    sel = np.isclose(C[:, 1], ycol[np.argmin(np.abs(ycol - 0.15))]) & (C[:, 0] > -0.03)
    tu = 100 * np.sqrt(2 * s.state["k"][sel] / 3) / U_T3A
    tu_ref, ratio = freestream_decay(C[sel, 0] + 0.04, U_T3A, NU_T3A, **T3A_INFLOW)
    assert np.max(np.abs(tu - 100 * tu_ref)) < 0.02
    assert 100 * freestream_decay(0.04, U_T3A, NU_T3A, **T3A_INFLOW)[0] == pytest.approx(
        3.35, abs=0.01)                                  # Tu au bord d'attaque


def test_channel_1d_close_to_sst():
    """Canal établi (1D, Re_τ = 395) : écoulement turbulent partout, le modèle γ reste
    proche du SST (γ ≈ 0.03 seulement dans la sous-couche visqueuse)."""
    ref = run_rans_channel("sst", 395.0).summary["Ub_plus"]
    r = run_rans_channel("sst_gamma", 395.0)
    assert r.summary["converged"]
    assert r.summary["Ub_plus"] == pytest.approx(ref, rel=0.02)
    g = r.solution.state["gamma"]
    assert g.max() == pytest.approx(1.0, abs=1e-6) and g[0] == pytest.approx(g[1], rel=0.01)


@pytest.mark.parametrize("conv", ["upwind", "linearUpwindLimited"])
def test_fake_gpu_matches_cpu(conv):
    """Triangles, entrée / sortie / symétrie / paroi : mêmes résultats CPU et faux GPU."""
    register_fake_gpu()
    dom = Rectangle(-2, -1.5, 4, 1.5, names={"left": "inlet", "right": "outlet",
                                             "bottom": "side", "top": "side"})
    m = triangulate(dom - Circle((0, 0), 0.3, "body").as_wall(), 0.35,
                    [{"shape": Circle((0, 0), 0.3), "h": 0.1, "growth": 0.3}], max_iter=60)
    out = {}
    for be in ("cpu", "fakegpu"):
        s = Solver2D(m, 1e-3, {"inlet": {"type": "inlet", "U": [1, 0]},
                               "outlet": {"type": "outlet"}, "side": {"type": "symmetry"},
                               "body": {"type": "wall"}}, model="sst_gamma", initial_U=(1, 0),
                     turbulence_inflow={"intensity": 0.03, "viscosity_ratio": 10.0},
                     settings=Settings(backend=be, solver_p="direct", convection_turb=conv))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s.run_steady(max_iter=25, tol=1e-12)
            ex = s.backend.to_host(s.model.extra_fields(s.state, s.flow()))
        s.to_cpu()
        out[be] = s, ex
    (c, ec), (g, eg) = out["cpu"], out["fakegpu"]
    assert np.max(np.abs(c.U - g.U)) < 1e-6
    for k in c.state:
        assert np.max(np.abs(c.state[k] - g.state[k])) <= 1e-6 * np.max(np.abs(c.state[k]))
    assert np.allclose(ec["Re_theta_c"], eg["Re_theta_c"])
    assert np.any(c.state["gamma"] < 0.5)          # la couche limite reste laminaire


@pytest.mark.parametrize("scheme", ["linearUpwind", "linearUpwindLimited"])
def test_second_order_turbulence_convection_on_orthogonal_mesh(scheme):
    """convection_turb = linearUpwind(Limited) agit aussi sur maillage orthogonal (elle y
    retombait silencieusement sur upwind, faute de gradient) ; indispensable à la
    transition : l'upwind avance la transition du cas T3A- de 36 %."""
    m = flat_plate_mesh(1.0, 0.1, 0.5, 4, 16, 16, first_height=1e-3)
    assert m.quality()["non_orthogonality_max_deg"] < 1e-6
    out = {}
    for conv in ("upwind", scheme):
        s = Solver2D(m, 1e-4, {**PLATE_BCS, "inlet": {"type": "inlet", "U": [1.0, 0.0]}},
                     model="sst_gamma", initial_U=(1.0, 0.0),
                     turbulence_inflow={"intensity": 0.03, "viscosity_ratio": 10.0},
                     settings=Settings(convection_turb=conv))
        s.run_steady(max_iter=20, tol=1e-12)
        out[conv] = s.state["k"]
    assert np.max(np.abs(out[scheme] - out["upwind"])) > 1e-3 * np.max(out["upwind"])
