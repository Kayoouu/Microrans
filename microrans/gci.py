"""Incertitude de discrétisation en maillage : indice de convergence de maillage (GCI) de
Roache, procédure de Celik et al., « Procedure for estimation and reporting of uncertainty
due to discretization in CFD applications », J. Fluids Eng. 130 (2008) 078001.

Trois maillages (1 : fin, 2, 3 : grossier) raffinés de façon systématique, une grandeur φ
calculée sur chacun : ordre apparent p, valeur extrapolée (Richardson) et GCI du maillage
fin (bande d'incertitude relative, facteur de sécurité 1.25).
"""
from __future__ import annotations

import math


def gci(phi, cells, dim: int, safety: float = 1.25, tol: float = 1e-12) -> dict:
    """phi = (φ1, φ2, φ3) du plus fin au plus grossier ; cells = (N1, N2, N3) nombres de
    cellules ; dim = 2 ou 3 (r = (N_fin / N_grossier)^(1/dim)).

    Renvoie r21, r32, p (ordre apparent), phi_ext (valeur extrapolée de 1 et 2), e_a
    (écart relatif 1–2), e_ext (écart relatif de φ1 à la valeur extrapolée), gci_fine
    (GCI relatif du maillage fin) et `oscillating` (convergence oscillante : écarts de
    signes opposés, p calculé avec leurs valeurs absolues, résultat à prendre avec
    prudence comme le recommande Celik)."""
    f1, f2, f3 = (float(v) for v in phi)
    n1, n2, n3 = (float(n) for n in cells)
    if not n1 > n2 > n3 > 0:
        raise ValueError("cells : nombres de cellules du plus fin au plus grossier, "
                         "strictement décroissants.")
    r21, r32 = (n1 / n2) ** (1.0 / dim), (n2 / n3) ** (1.0 / dim)
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0.0 or e32 == 0.0:
        raise ValueError("deux maillages donnent la même valeur : ordre apparent indéfini "
                         "(grandeur déjà convergée, ou insensible au maillage).")
    s = math.copysign(1.0, e32 / e21)
    ratio = abs(e32 / e21)
    p = abs(math.log(ratio)) / math.log(r21)       # point de départ : q = 0
    for _ in range(200):                           # point fixe de Celik (éq. 5)
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(ratio) + q) / math.log(r21)
        if abs(p_new - p) < tol:
            p = p_new
            break
        p = p_new
    k = r21 ** p
    phi_ext = (k * f1 - f2) / (k - 1.0)
    e_a = abs((f1 - f2) / f1)
    return {"r21": r21, "r32": r32, "p": p, "phi_ext": phi_ext, "e_a": e_a,
            "e_ext": abs((phi_ext - f1) / phi_ext), "gci_fine": safety * e_a / (k - 1.0),
            "oscillating": s < 0}
