"""Indice de convergence de maillage (procédure de Celik et al. 2008)."""
import pytest

from microrans.gci import gci


def test_published_example():
    """Exemple de calcul publié par Celik et al. (2008), tableau 1, première colonne
    (maillages 2D de 18 000, 8 000 et 4 500 cellules)."""
    r = gci((6.063, 5.972, 5.863), (18000, 8000, 4500), 2)
    assert r["r21"] == pytest.approx(1.5) and r["r32"] == pytest.approx(4 / 3)
    assert round(r["p"], 2) == 1.53
    assert round(r["phi_ext"], 4) == 6.1685
    assert round(100 * r["e_a"], 1) == 1.5
    assert round(100 * r["e_ext"], 1) == 1.7
    assert round(100 * r["gci_fine"], 1) == 2.2
    assert not r["oscillating"]


def test_exact_second_order_data():
    """φ(h) = φ₀ + C h² exactement : ordre 2 et valeur extrapolée φ₀ (2D et 3D)."""
    for dim, cells in ((2, (1600, 400, 100)), (3, (8000, 1000, 125))):
        phi = [1.25 + 0.3 * h ** 2 for h in (0.25, 0.5, 1.0)]
        r = gci(phi, cells, dim)
        assert r["r21"] == pytest.approx(2.0)
        assert r["p"] == pytest.approx(2.0, abs=1e-9)
        assert r["phi_ext"] == pytest.approx(1.25, abs=1e-12)


def test_oscillating_and_invalid_inputs():
    assert gci((1.0, 1.1, 0.95), (400, 100, 25), 2)["oscillating"]
    with pytest.raises(ValueError, match="même valeur"):
        gci((1.0, 1.0, 1.2), (400, 100, 25), 2)
    with pytest.raises(ValueError, match="strictement décroissants"):
        gci((1.0, 1.1, 1.2), (25, 100, 400), 2)
