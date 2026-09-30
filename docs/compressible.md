# Écoulements compressibles : solveur en densité (Euler, Navier-Stokes laminaire)

Document technique (le README en reprend l'essentiel). Tous les chiffres ci-dessous ont été
mesurés avec ce code ; les écarts aux références sont donnés tels quels. Tableaux
reproductibles : `python -m microrans.fv2d.compressible_validation --out docs`
(4.5 min sans le NACA, mesuré ; le NACA ajoute ~15 min, estimé d'après les temps
par itération du § 4).

## 1. Ce que c'est, ce que ce n'est pas

- **2D plan**, gaz parfait (γ, R, Pr constants ; air sec par défaut), grandeurs **SI
  dimensionnelles** (Pa, K, kg/m³, m/s). Mêmes maillages que le solveur incompressible
  (structurés, triangles, hybrides, importés, périodicité comprise).
- **Euler** (non visqueux) ou **Navier-Stokes laminaire** (μ constante ou loi de
  Sutherland, conduction de Fourier, dissipation visqueuse).
- Stationnaire (pas de temps local : Runge-Kutta ou **implicite**) et instationnaire
  (Runge-Kutta SSP d'ordre 3, pas global).
- **Pas de turbulence** (ni RANS ni LES) : un écoulement compressible turbulent n'est pas
  calculable. **Pas d'axisymétrique**, **pas de carte graphique** (CPU seulement), pas de
  gaz réel, de combustion, ni de multi-espèces.

Activation : `[physics] compressible = true` ; `run_case` délègue alors tout le calcul à
`fv2d/compressible_case.py`. Mêmes sorties que l'incompressible (`summary.json`,
`history.csv`, `wall_<patch>.csv` avec Cp, C_f, y⁺, T, q ; `fields.vtk` avec ρ, U, p, T,
Mach, Cp, entropie ; figures ; `checkpoint.npz` et reprise exacte ou interpolée).

## 2. Méthodes (et pourquoi)

| Élément | Choix | Référence |
|---|---|---|
| Flux convectifs | Roe (défaut) avec correction d'entropie de Harten (δ = 0.1 c̃) ; HLLC en option | Roe 1981 ; Blazek (2015) éq. 4.89-4.91 ; Toro (2009) § 10.4 |
| Ordre 2 | reconstruction MUSCL des variables primitives, gradient de Green-Gauss | Blazek § 5.3 |
| Limiteur | Venkatakrishnan (seuil ε de Wang : indépendant des unités), Barth-Jespersen, aucun ; gel optionnel (`limiter_freeze`, **à éviter** : voir § 3.5) | Venkatakrishnan 1995 ; Wang 2000 ; SU2 |
| Flux visqueux | gradient aux faces moyen + correction selon PN | Blazek § 5.4 ; Weiss et al. 1999 |
| Temps (instationnaire) | SSP-RK3, pas global au CFL | Shu & Osher 1988 |
| Stationnaire | pas local ; RK3, RK5 (van Leer, Tai & Powell 1989) ou **implicite** : Euler implicite linéarisé, jacobienne d'ordre 1 (dissipation de Roe exacte), GMRES préconditionné par Gauss-Seidel symétrique par blocs, CFL adaptatif (plafond divisé par 2 seulement si la solution oscille sans converger : § 3.5) | Blazek § 6.2 ; SU2 |
| Conditions | champ lointain (invariants de Riemann) + **correction de tourbillon ponctuel** optionnelle, entrée subsonique (p0, T0, direction), sortie subsonique (p), entrée / sortie supersoniques, paroi glissante, paroi adhérente adiabatique ou isotherme (mobile), symétrie, périodicité | Blazek § 8.4-8.5 ; Thomas & Salas 1986 |

**Correction de tourbillon ponctuel** (`[boundary.<champ lointain>] vortex = [x, y]`, point
d'application, ex. quart de corde) : la vitesse induite par la circulation
Γ = L′/(ρ∞U∞) (portance actuelle des parois, Kutta-Joukowski) est ajoutée à l'état amont
des faces de champ lointain, avec le facteur de Prandtl-Glauert ; pression et masse
volumique suivent à enthalpie totale et entropie amont. Mesuré (NACA 0012, M = 0.5,
α = 1.25°, 64 cellules autour du profil) :

| Distance du champ lointain | 10 cordes | 30 cordes | 100 cordes | écart |
|---|---:|---:|---:|---:|
| C_l sans correction | 0.1590 | 0.1634 | 0.1651 | 3.8 % |
| C_l avec correction | 0.1655 | 0.1657 | 0.1658 | 0.2 % |

## 3. Validation

### 3.1 Tube à choc de Sod (instationnaire, Euler)

Bande de 1 maille d'épaisseur, t = 0.2 L/√(p_L/ρ_L), solution exacte (problème de Riemann,
Toro 2009 chap. 4).

| Flux | Mailles | L1(ρ) | L1(u)/a | L1(p)/p_L | choc (exact 0.8504) | contact (0.6855) | milieu de détente (0.3747) | masse |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Roe | 100 | 6.03e-3 | 1.07e-2 | 4.67e-3 | 0.8514 | 0.6848 | 0.3771 | −9e-10 |
| Roe | 200 | 3.33e-3 | 5.49e-3 | 2.41e-3 | 0.8510 | 0.6858 | 0.3759 | −1e-14 |
| Roe | 400 | 1.92e-3 | 2.90e-3 | 1.26e-3 | 0.8508 | 0.6861 | 0.3753 | −2e-14 |
| HLLC | 100 | 6.32e-3 | 1.12e-2 | 4.90e-3 | 0.8515 | 0.6845 | 0.3769 | −8e-10 |
| HLLC | 200 | 3.48e-3 | 5.74e-3 | 2.54e-3 | 0.8510 | 0.6856 | 0.3758 | −1e-14 |
| HLLC | 400 | 1.99e-3 | 3.02e-3 | 1.32e-3 | 0.8507 | 0.6860 | 0.3752 | −2e-14 |

Ordre observé en L1(ρ) : 0.86 puis 0.80 — attendu inférieur à 1 avec des discontinuités
(l'ordre 2 ne vaut que dans les zones régulières). Roe est un peu plus précis que HLLC.
Petit dépassement (~1 %) en tête et en queue de détente (limiteur de Venkatakrishnan,
qui n'est pas strictement monotone).

![Sod](compressible_sod.png)

### 3.2 Rampe de compression supersonique (stationnaire, Euler)

M = 2, θ = 10°, théorie du choc oblique (NACA Report 1135) : β = 39.31°, p2/p1 = 1.7066,
M2 = 1.6405.

| Cellules | β | p2/p1 | M2 | C_d de la rampe (exact 0.0445) | Itérations (implicite) |
|---:|---:|---:|---:|---:|---:|
| 1 350 | 39.28° | 1.708 | 1.639 | 0.04449 | 55 |
| 5 400 | 39.33° | 1.707 | 1.640 | 0.04449 | 70 |
| 21 600 | 39.30° | 1.707 | 1.640 | 0.04450 | 87 |

![Rampe](compressible_rampe.png)

### 3.3 Plaque plane laminaire (Navier-Stokes)

M = 0.2, Re_L = 1e5, paroi adiabatique ; Blasius C_f = 0.664/√Re_x (effet de
compressibilité < 1 % à ce Mach).

| Cellules | C_d (Blasius 0.00420) | écart C_f x = 0.1 | 0.4 | 0.8 | 0.95 | écart moyen (0.1 < x < 0.95) |
|---:|---:|---:|---:|---:|---:|---:|
| 3 840 | 0.004236 | −2.0 % | −0.3 % | +1.3 % | +1.8 % | 1.0 % |
| 15 360 | 0.004231 | −1.4 % | +0.1 % | +1.7 % | +2.7 % | 0.8 % |

L'écart près du bord de fuite (+2.7 %) **ne diminue pas** avec le maillage : c'est la
sortie à pression imposée placée au bord de fuite (x = 1), qui perturbe la couche limite
sur les derniers pourcents ; il faudrait prolonger le domaine en aval.

![Plaque](compressible_plaque.png)

### 3.4 Écoulement de Couette avec dissipation visqueuse

μ constante, parois isothermes, T = T_w + Pr U²/(2c_p) η(1 − η) (White, « Viscous Fluid
Flow », § 3-2) : vitesse exacte à 1e-8, température à 0.2 % de ΔT_max, frottement et flux
de chaleur pariétaux à 3e-8 (201 itérations implicites, 0.8 s).

### 3.5 Profil NACA 0012 (Euler)

Profil à bord de fuite fermé (coefficient −0.1036 de l'équation NACA), maillage en O,
champ lointain à 30 cordes avec correction de tourbillon, Roe + MUSCL + Venkatakrishnan
(limiteur **non gelé**), schéma implicite ; hauteur de 1re maille 2e-3 × 192 / n_autour.
Tous les calculs ci-dessous sont convergés (résidus relatifs < 1e-8).

**M = 0.8, α = 1.25°** (transsonique : choc à l'extrados vers x = 0.63, choc faible à
l'intrados vers x = 0.35)

| Maillage | Cellules | C_l | C_d | C_m (¼ de corde, cabreur > 0) | Itérations |
|---|---:|---:|---:|---:|---:|
| 96 × 32 | 3 072 | 0.3139 | 0.02217 | −0.0284 | 377 |
| 192 × 64 | 12 288 | 0.3353 | 0.02208 | −0.0341 | 1 011 |
| 384 × 128 | 49 152 | 0.3345 | 0.02186 | −0.0337 | 1 726 |
| valeurs publiées (voir ci-dessous) | | ≈ 0.35 | ≈ 0.022-0.023 | | |

![NACA 0012 M = 0.8](compressible_naca0012_M0.8.png)

**C_l plafonne vers 0.335, environ 4.5 % sous la valeur publiée ; C_d est ~3 % en
dessous. L'écart n'est pas expliqué.** Ce n'est ni la convergence (résidus < 1e-8) ni la
finesse du maillage (−0.2 % entre 192 × 64 et 384 × 128), ni les choix numériques testés
sur 192 × 64 :

| Variante (192 × 64) | C_l | C_d | écart de C_l |
|---|---:|---:|---:|
| réglages de l'exemple (Roe, correction d'entropie 0.1, K = 0.05, 30 cordes) | 0.3353 | 0.02208 | — |
| champ lointain à 100 cordes (72 mailles radiales) | 0.3355 | 0.02216 | +0.06 % |
| flux HLLC | 0.3352 | 0.02208 | −0.03 % |
| sans correction d'entropie | 0.3351 | 0.02207 | −0.06 % |
| limiteur plus souple (K = 0.3) | 0.3337 | 0.02190 | −0.5 % |

Non testé : la forme du bord de fuite (les études de référence utilisent soit le bord de
fuite ouvert d'origine — AGARD-AR-211, 1985 —, soit un bord de fuite pointu obtenu en
prolongeant la corde — Vassberg & Jameson 2010, J. Aircraft 47(4)), la topologie du
maillage au bord de fuite, un maillage plus fin que 384 × 128. **Les valeurs de référence
elles-mêmes n'ont pas pu être revérifiées** (articles inaccessibles depuis la machine de
développement) : « ≈ 0.35 » et « 0.022-0.023 » sont des ordres de grandeur connus, pas
des chiffres relus.

**M = 0.5, α = 1.25°** (subsonique : traînée nulle en Euler, toute traînée calculée est
une erreur numérique)

| Maillage | Cellules | C_l | C_d (exact : 0) | Itérations |
|---|---:|---:|---:|---:|
| 96 × 32 | 3 072 | 0.1715 | 0.00277 | 157 |
| 192 × 64 | 12 288 | 0.1771 | 0.00055 | 286 |
| 384 × 128 | 49 152 | 0.1791 | 0.00017 | 533 |

La traînée parasite est divisée par 5.0 puis 3.2 à chaque raffinement (ordre 2.3 puis
1.7) ; C_l converge (écarts 0.0056 puis 0.0020 : ordre 1.5, valeur extrapolée ≈ 0.180).

**Convergence : deux pièges, corrigés** (mesurés sur 192 × 64, M = 0.8) :

1. **Limiteur gelé trop tôt = solution fausse.** Avec `limiter_freeze = 200` (ancien
   réglage de l'exemple), le limiteur est figé alors que le choc n'est pas encore en place
   (C_l = 0.296 à l'itération 200) : les résidus tombent ensuite à 1e-10, mais vers
   **C_l = 0.3238 (−3.4 %)**. Le limiteur de Venkatakrishnan converge ici sans être gelé ;
   ne geler (`limiter_freeze`) que si les résidus cyclent, et seulement une fois les
   efforts stabilisés.
2. **Plafond de CFL réduit pendant le démarrage.** Le schéma implicite divisait le
   plafond de CFL par 2 dès que les résidus stagnaient 25 itérations — ce qui arrive
   normalement au démarrage, pendant que le choc se déplace. Le CFL tombait à 6.25 et
   n'en remontait pas : pas de convergence en 3 000 itérations. Le plafond n'est plus
   réduit que si la solution **oscille** pendant la stagnation :
   ‖Q_n − Q_début‖ / Σ‖ΔQ_k‖ < 0.2 (mesuré : 0.38 à 0.92 aux stagnations du démarrage,
   0.03 à 0.11 dans le vrai cycle limite à CFL 100 ; seuil calibré sur ce cas et sur la
   plaque plane seulement). `cfl_cuts = 0` supprime toute réduction.

| Pilotage du CFL (limiteur non gelé) | C_l à ±0.001 dès l'itération | résidus < 1e-6 dès |
|---|---:|---:|
| ancienne règle, plafond 50 (réduit à 6.25) | non atteint en 1 500 | non atteint en 1 500 |
| CFL fixe 25 | 690 | non atteint en 1 000 (3e-6) |
| CFL fixe 50 | 390 | 630 |
| CFL fixe 100 | 250 | jamais (cycle limite, ~3e-3) |
| **nouvelle règle, plafond 100 (exemple)** | **250** | **510** |

Test de non-régression : `test_implicit_cfl_cap_kept_during_transonic_startup` (échoue
avec l'ancienne règle : 2 réductions en 120 itérations).

**Deux solutions stationnaires selon le schéma d'avancement** (96 × 32, M = 0.8) :
implicite et RK3 convergent tous deux (résidus 1e-8 et 1e-10) mais vers des solutions
différentes, C_l = 0.3139 et 0.3162 (0.7 %), C_d = 0.02217 et 0.02229. Ils redonnent
exactement la même solution sans correction de tourbillon (0.3069 / 0.3069) ou à l'ordre 1
(0.2627 / 0.2627) : c'est la combinaison reconstruction limitée d'ordre 2 + champ lointain
dépendant de la portance qui admet deux états stationnaires. Mécanisme non établi ; non
mesuré sur les maillages plus fins (RK3 trop lent). N'explique pas l'écart de 4.5 % ci-dessus.

## 4. Performances (un cœur)

Machine de développement (4 cœurs logiques Intel Xeon 2.1 GHz virtualisés), un seul
calcul à la fois ; mémoire = maximum résident du processus (Python et bibliothèques
compris, ~120 Mo).

| Cas | Cellules | Schéma | ms / itération | µs / itération / cellule | Mémoire |
|---|---:|---|---:|---:|---:|
| NACA 0012 (Euler) | 3 072 | implicite | 23.8 | 7.7 | 144 Mo |
| NACA 0012 (Euler) | 12 288 | implicite | 98.7 | 8.0 | 200 Mo |
| NACA 0012 (Euler) | 49 152 | implicite | 449 | 9.1 | 422 Mo |
| NACA 0012 (Euler) | 12 288 | RK3 (pas local) | 44.2 | 3.6 | 151 Mo |
| Plaque plane (Navier-Stokes) | 3 840 | implicite | 34.9 | 9.1 | 149 Mo |
| Plaque plane (Navier-Stokes) | 3 840 | RK3 (pas local) | 19.3 | 5.0 | 131 Mo |

Une itération implicite coûte ~2 itérations RK3 (assemblage de la jacobienne, GMRES),
mais il en faut beaucoup moins : NACA 96 × 32 convergé (résidus < 1e-6) en 213
itérations implicites (~5 s) contre 16 136 itérations RK3 (177 s, mesuré avec un autre
calcul en parallèle). Durées totales (implicite, résidus
< 1e-8) : ~1 min 40 s pour 12 288 cellules, ~13 min pour 49 152 cellules.

## 5. Limites

1. **Pas de turbulence** (ni RANS ni LES), **pas d'axisymétrique**, **CPU seulement**,
   gaz parfait. Un écoulement compressible turbulent (profil à Reynolds de vol, tuyère
   réelle) n'est pas calculable.
2. **Transsonique : C_l du NACA 0012 à M = 0.8 environ 4.5 % sous les valeurs publiées**,
   écart non expliqué (§ 3.5) ; la position du choc et la forme de C_p sont correctes
   qualitativement, pas quantitativement vérifiées.
3. **`limiter_freeze`** : un gel avant que la solution soit établie donne une solution
   convergée mais fausse, sans que les résidus le signalent (§ 3.5).
4. **Pilotage du CFL implicite** calibré sur deux cas (NACA transsonique, plaque plane).
   Si les résidus cyclent : baisser `cfl_max` ; `cfl_cuts` règle le nombre de réductions
   automatiques.
5. **Plaque plane** : C_f +2.7 % près du bord de fuite, qui ne diminue pas avec le maillage
   (sortie placée au bord de fuite, § 3.3).
6. **Bas Mach** : flux de Roe / HLLC sans préconditionnement. Mesuré correct à M = 0.2
   (plaque plane) ; plus bas, non mesuré — la littérature (Guillard & Viozat 1999)
   montre une perte de précision quand M → 0. En dessous de M ≈ 0.3, le solveur
   incompressible est le bon outil.
7. **Coût** : ~8-9 µs par itération et par cellule en implicite ; au-delà de ~50 000
   cellules, compter plus de 10 minutes par calcul stationnaire.
8. **Solution stationnaire pas toujours unique** : avec la correction de tourbillon, le
   schéma implicite et RK3 convergent vers deux solutions distantes de 0.7 % en C_l
   (96 × 32, § 3.5).
9. **Couverture de la validation** : pas de cas visqueux avec choc (interaction choc /
   couche limite), un seul cas instationnaire (Sod).
