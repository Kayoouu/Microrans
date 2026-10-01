"""Résumé lisible d'un calcul (audit utilisateur U7 : JSON brut à 16 chiffres)."""
from microrans.fv2d.report import summary_text


def _base(**kw):
    s = {"mode": "steady", "model": "laminar", "n_cells": 4096, "axisymmetric": False,
         "nu": 0.01, "reference_velocity": 1.0, "reference_length": 1.0,
         "angle_of_attack": 0.0, "converged": True, "iterations": 351, "wall_time_s": 10.31,
         "U_mean": [7.896858213370121e-05, 9.41440091064883e-06], "backend": "cpu",
         "lid": {"Cd": -0.3835019112103236, "Cl": 0.0628097540249018,
                 "Cm": -0.5388562077760898, "Cd_pressure": 0.0,
                 "Cd_viscous": -0.3835019112103236, "yplus_max": 0.764593863242399,
                 "yplus_mean": 0.31}}
    s.update(kw)
    return s


def test_steady_summary():
    t = summary_text(_base())
    assert t.startswith("Calcul stationnaire, laminaire — 4 096 cellules.")
    assert "Convergé en 351 itérations (10.3 s)." in t
    assert "Re = U L / ν = 100." in t
    assert "  lid         -0.3835   0.06281   -0.5389   0.7646" in t
    assert "lid : Cd = 0 (pression) − 0.3835 (frottement)." in t
    assert "0.3835019112103236" not in t                     # 4 chiffres, pas 16


def test_not_converged_and_unknown_keys():
    t = summary_text(_base(converged=False, iterations=3000, nouvelle_cle=1.23456789))
    assert "NON CONVERGÉ après 3 000 itérations" in t
    assert "Autres valeurs : nouvelle_cle = 1.235" in t


def test_transient_with_strouhal():
    s = _base(mode="transient", time=50.0, steps=5000)
    s.pop("converged")
    s["lid"].update(Cd_mean=1.3286, Cl_rms=0.2399, Cl_amplitude=0.3395, strouhal=0.16011,
                    periods_used=15)
    t = summary_text(s)
    assert "Temps simulé t = 50 en 5 000 pas (10.3 s)." in t
    assert "Strouhal 0.1601 (15 périodes)" in t


def test_compressible_summary():
    s = {"solver": "compressible", "mode": "steady", "model": "euler", "n_cells": 12288,
         "flux": "roe", "order": 2, "limiter": "venkatakrishnan", "converged": True,
         "iterations": 510, "wall_time_s": 254.48, "angle_of_attack": 1.25,
         "freestream": {"mach": 0.8, "speed": 272.24, "p": 101325.0, "T": 288.15,
                        "rho": 1.225},
         "final_residuals": {"rho": 8.9e-07},
         "airfoil": {"Cd": 0.022076590579364997, "Cl": 0.3352655151046101,
                     "Cm": -0.03405076987914659, "Cd_pressure": 0.0220765, "Cd_viscous": 0.0}}
    t = summary_text(s)
    assert t.startswith("Calcul compressible stationnaire, Euler (non visqueux)")
    assert "Amont : Mach 0.8, vitesse 272.2 m/s" in t
    assert "force / (½ ρ U² L)" in t and "airfoil     0.02208    0.3353  -0.03405" in t
