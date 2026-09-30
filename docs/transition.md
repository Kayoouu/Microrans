# Transition laminaire-turbulent : modèle γ de Menter (2015) couplé au k-ω SST

Document technique détaillé (le README en reprend l'essentiel). Tous les chiffres
ci-dessous ont été mesurés avec ce code ; les écarts aux références sont donnés tels quels.

## 1. Ce qui a été ajouté

- Modèle `sst_gamma` (alias `sst-gamma`, `gamma`, `transition`), libellé
  « k-ω SST + transition γ (Menter 2015) » : `microrans/models/transition_gamma.py`.
  Utilisable en 2D (CPU et GPU, vérifié avec le faux GPU) et en 1D (canal).
- Option `[physics] model_options = { kato_launder = true }` : production de k en ν_t S Ω.
- Fonction `freestream_decay(distance, U, ν, Tu, ν_t/ν)` (même fichier) : décroissance de la
  turbulence amont prédite par le SST, pour régler Tu et ν_t/ν d'entrée (voir § 4).
- Exemple `microrans/examples/plaque_plane_transition_t3a.toml` (T3A, ~30-50 s).
- Tests `tests/test_transition.py` (11 tests, ~11 s au total).
- Interface graphique : entrée ajoutée à `MODEL_LABELS` (visible aussi dans l'onglet 1D) ;
  page « Numérique » : liste « Convection turbulence » (upwind par défaut, inchangé ;
  linearUpwindLimited conseillé en transition). Test d'aller-retour dans
  `tests/test_gui_forms.py`.
- Figure `docs/transition_plaques_t3.png`.
- Petits ajouts hors du modèle :
  - `models/base.py` : attribut `wall_zero_gradient` (variables à gradient nul aux parois) ;
  - `fv2d/solver.py`, `scalar_bc` : respecte `wall_zero_gradient` (γ) ;
  - `fv2d/solver.py`, `Solver2D.__init__` : lois de paroi refusées avec `sst_gamma` ;
  - `fv2d/solver.py`, `Step2D.solve` : **correction d'un défaut existant** — sur maillage
    orthogonal, `convection_turb = "linearUpwind"` retombait silencieusement sur `upwind`
    (le gradient n'était calculé que pour la correction non orthogonale) ; en outre
    `linearUpwindLimited` (limiteur de Barth-Jespersen, déjà utilisé pour les scalaires) est
    maintenant accepté pour les équations de turbulence. Défaut inchangé (`upwind`), donc
    aucun résultat existant ne change ;
  - `mesh2d/mesh.py` : `Mesh2D.wall_normal()` (direction du point de paroi le plus proche,
    = ∇d), sans changer `distance_to_patches`.

## 2. Équations implémentées (forme de l'article, corrélation fondée sur Tu)

Réf. : F.R. Menter, P.E. Smirnov, T. Liu, R. Avancha, *A One-Equation Local
Correlation-Based Transition Model*, Flow Turbulence Combust. 95 (2015) 583-619. L'article
n'était pas accessible depuis la machine de travail : les équations ont été écrites de
mémoire puis **vérifiées ligne à ligne contre l'implémentation libre de SU2** (option
« SLM / MENTER_SLM », branche de la PR su2code/SU2#1901 : `trans_sources.hpp`,
`trans_correlations.hpp`, couplage dans `turb_sources.hpp`) — mêmes formes et constantes.

```
Dγ/Dt = P_γ − E_γ + ∇·[(ν + ν_t/σ_f) ∇γ]
P_γ = F_length S γ (1 − γ) F_onset          E_γ = c_a2 Ω F_turb γ (c_e2 γ − 1)
F_onset = max(min(Re_V/(2.2 Re_θc), 2) − max(1 − (R_T/3.5)³, 0), 0)
F_turb = exp(−(R_T/2)⁴)    Re_V = d² S/ν    R_T = k/(ν ω)
Re_θc = 100 + 1000 exp(−Tu_L F_PG(λ_θL))    Tu_L = min(100 √(2k/3)/(ω d), 100)
λ_θL = clip(−7.57e-3 (dV/dy) d²/ν + 0.0128, −1, 1),  dV/dy = n·∇(n·U) = nᵢ nⱼ ∂uᵢ/∂xⱼ
F_PG = min(1 + 14.68 λ, 1.5) si λ ≥ 0 ; min(1 − 7.34 λ, 3) sinon ; F_PG ≥ 0
F_length = 100, c_a2 = 0.06, c_e2 = 50, σ_f = 1
k :  P̃_k = γ P_k + P_k^lim,   D̃_k = max(γ, 0.1) β* k ω
     P_k^lim = 5 max(γ − 0.2, 0)(1 − γ) F_on^lim max(3ν − ν_t, 0) S Ω
     F_on^lim = min(max(Re_V/(2.2·1100) − 1, 0), 3)
F1 = max(F1_SST, exp(−(R_y/120)⁸)), R_y = d √k/ν ;  équation de ω inchangée (SST 2003)
```

Choix de mise en œuvre : n = vecteur unitaire du point de paroi le plus proche vers le
centre de la cellule (calculé une fois, géométriquement ; n·∇n = 0 le long de la normale,
d'où dV/dy = nᵢnⱼ∂uᵢ/∂xⱼ). Sources de γ linéarisées par Newton (`linearize_source`), γ
borné à [0, 1]. Conditions aux limites : γ = 1 en entrée / champ lointain, gradient nul à
la paroi (en 1D : valeur du premier nœud intérieur, décalée d'une itération). En 1D,
dV/dy = 0 (canal). Non implémentés : rugosité, écoulement transverse (crossflow), variante
Spalart-Allmaras, corrections de courbure.

## 3. Texte proposé pour le README

### § 4, tableau « méthodes » (nouvelle ligne)

| Transition | γ à une équation de Menter et al. (2015) couplé au SST : corrélation locale Re_θc(Tu_L, λ_θL), P_k × γ, D_k × max(γ, 0.1), P_k^lim ; γ à gradient nul en paroi | Menter, Smirnov, Liu & Avancha (2015) ; implémentation SU2 « SLM » |

Et, dans la ligne « Turbulence » : ajouter « SST + transition γ (Menter 2015) ».

### § 6, tableau de validation 2D (nouvelles lignes)

Maillage 13 760 cellules (1re maille 5e-5 m ; y⁺ = 1.6 au bord d'attaque, ≤ 0.7 pour
x > 1 cm, ≤ 0.5 pour x > 5 cm), convection de
la turbulence `linearUpwindLimited`, entrée à 0.04 m en amont du bord d'attaque.

| Cas | Grandeur | microrans | Référence |
|-----|----------|-----------|-----------|
| Plaque T3A (ERCOFTAC), U = 5.4 m/s, Tu_BA = 3.35 % (décroissance ajustée sur les mesures), SST + γ | Re_x du minimum de C_f / de mi-transition / fin (max C_f) | 1.47e5 / 1.91e5 / 2.9e5 | ≈ 1.42e5 / 2.28e5 / 2.9-3.2e5 (Savill 1993, un point tous les 0.1 m) |
| idem | C_f max ; C_f(x = 1.495 m) | 0.00449 ; 0.00403 | 0.00486 ; 0.00408 (corrélation 0.0576 Re_x^-0.2 : 0.00411) |
| idem | C_f / Blasius avant transition (x ≤ 0.4 m) | 1.10 à 1.26 | mesures : 1.00 à 1.19 |
| idem, maillage fin (27 720 cellules, y⁺ ≤ 0.9) | minimum / mi-transition / C_f max | 1.50e5 / 1.95e5 / 0.00455 | écart de maillage ≤ 2 % |
| Plaque T3A-, U = 19.8 m/s, Tu_BA = 0.85 % (entrée 0.874 %, ν_t/ν = 8.72) | Re_x du minimum de C_f ; 90 % de la transition | 1.36e6 (maillage fin : 1.42e6) ; 1.45e6 (fin : 1.50e6) | ≈ 1.45e6 ; > 2.0e6 (C_f ≈ 0.0015 à 2.0e6, soit ~38 %) — lecture graphique |
| Plaque T3B, U = 9.4 m/s, Tu_BA = 6.1 % (entrée 6.5 %, ν_t/ν = 100) | Re_x du minimum / du maximum de C_f ; C_f min | 7.5e4 / 1.5e5 ; 0.0045 | ≈ 6e4 / 1.25e5 ; ≈ 0.0034 — lecture graphique |

Sources des données : T3A — fichier `T3A.dat` (C_f et Tu tous les 0.1 m, données ERCOFTAC
de Savill 1993) distribué avec le tutoriel OpenFOAM `simpleFoam/T3A`. T3A- et T3B — **pas de
fichier de données accessible** : valeurs lues sur les figures de vérification du modèle
Langtry-Menter de SU2 (`vandv_files/LM_model/T3Am/All_Cf.png`, `.../T3B/All_Cf.png`, points
« Exp » ERCOFTAC), précision de lecture ≈ ±5 %.

Canal 1D Re_τ = 395 (écoulement entièrement turbulent) : U_b⁺ = 17.18 avec SST + γ contre
17.38 avec le SST seul (−1.1 %) et 17.20 (corrélation de Dean) ; γ ≈ 0.03 dans la
sous-couche visqueuse, 1 ailleurs.

### § 7, performances (un cœur, machine partagée chargée : temps majorés)

| Cas | Cellules | Itérations | Temps |
|---|---:|---:|---:|
| T3A (exemple) | 13 760 | 449 | 35 s (28 s mesurés avec `upwind`, machine moins chargée) |
| T3A fin | 27 720 | 720 | 139 s |
| T3A- | 13 760 / 27 720 | 1 343 / 1 424 | 110 s / 238 s |
| T3B | 13 760 | 550 | 62 s |

### § 8, limites (texte proposé, remplace « pas de transition » du point 4)

**Transition (modèle γ, `sst_gamma`)** : validée seulement sur plaques planes sans gradient
de pression (T3A, T3A-, T3B) ; pas de cas avec gradient de pression, décollement laminaire
ou profil (le terme P_k^lim et λ_θL ne sont donc pas validés). Mesuré :
- le début de transition est bien placé sur T3A (+3 % en Re_x) et T3A- (−6 %, −2 % sur
  maillage fin), mais la
  transition prédite est **trop raide** : sur T3A la mi-transition arrive 16 % trop tôt et
  le maximum de C_f est 8 % trop bas ; sur T3A- la transition est terminée à Re_x ≈ 1.5e6
  alors que les mesures montrent une montée lente encore inachevée à 2.0e6 ;
- avant la transition, C_f est 10 à 26 % au-dessus de Blasius sur T3A (mesures : 0 à 19 %) :
  la viscosité turbulente de l'écoulement libre (ν_t/ν ≈ 12) pénètre la couche limite
  laminaire (le γ ne réduit que la production de k). Sur T3B, le « creux » laminaire n'est
  pas reproduit (C_f min 0.0045 contre ≈ 0.0034) ; à faible Tu (T3A-), C_f laminaire = Blasius
  à +2 à +4 % ;
- **très sensible à la turbulence amont** : 0.3 point de Tu au bord d'attaque (3.04 → 3.35 %)
  déplace le début de transition de 17 % (T3A). Il faut régler Tu et ν_t/ν d'entrée pour
  retrouver la décroissance mesurée de Tu (§ 4 ci-dessous) ; sans mesure de décroissance
  (T3A-, T3B ici), le résultat dépend du ν_t/ν choisi ;
- **sensible au schéma de convection de la turbulence** : `upwind` (défaut du code) avance
  la transition de T3A- de 36 % (Re_x 8.8e5 au lieu de 1.36e6) et la retarde de ~10 % sur
  T3A ; utiliser `convection_turb = "linearUpwindLimited"` (converge aussi bien ;
  `linearUpwind` non limité ne converge pas sur T3A- : résidus bloqués à 1e-3) ;
- y⁺ ≈ 1 obligatoire (lois de paroi refusées) ; en turbulent établi, γ ≈ 0.03 dans la
  sous-couche visqueuse (propriété du modèle) : C_f(Re_θ) ≈ 2 à 3 % au-dessus du SST seul
  au même Re_θ (T3A aval, maillage 13 760 cellules, `upwind` : +0.9 à +1.9 % sur la
  corrélation de Kármán-Schoenherr, SST seul −0.2 à −2.7 %).

### § 9, feuille de route

Remplacer « transition γ-Re_θ » par : « transition : cas avec gradient de pression et
décollement laminaire (T3C, profils à bas Reynolds), variantes rugosité / crossflow ».

### § 10, « Ajouter un modèle de turbulence »

Ajouter : « une variable à gradient nul aux parois (au lieu d'une valeur imposée) se déclare
dans `wall_zero_gradient` (ex. γ du modèle de transition). »

## 4. Régler la turbulence amont (Tu, ν_t/ν) pour un cas de transition

Dans l'écoulement libre uniforme, le SST (F1 = 0, donc β = β₂ = 0.0828) donne, à une
distance x de l'entrée :

```
Tu(x) = Tu₀ (1 + β ω₀ x/U)^(−β*/(2β)),   ω₀ = k₀ / (ν · (ν_t/ν)₀),   k₀ = 1.5 (Tu₀ U)²
```

(`microrans.models.transition_gamma.freestream_decay`, vérifié contre le calcul 2D : écart
< 0.004 point de Tu sur T3A, T3B et T3A-). Méthode : placer l'entrée à une distance connue du bord d'attaque,
puis choisir (Tu₀, ν_t/ν) pour que Tu(x) passe par les mesures (moindres carrés si une
décroissance mesurée existe, sinon au moins Tu au bord d'attaque). Pour T3A, entrée à
x = −0.04 m :

| Réglage | Tu entrée | ν_t/ν entrée | Tu bord d'attaque | Tu à x = 0.395 m (mesure 2.00 %) | Écart max à la décroissance mesurée |
|---|---:|---:|---:|---:|---:|
| Langtry & Menter (2009), tutoriel OpenFOAM | 3.3 % | 12 | 3.04 % | 1.90 % | −0.23 point |
| ajusté (utilisé dans l'exemple) | 3.7 % | 12.4 | 3.35 % | 1.98 % | 0.05 point (x ≤ 1.2 m) |

Valeurs de référence citées dans la littérature (à vérifier sur la source avant usage) :
T3A U∞ = 5.4 m/s, Tu ≈ 3 à 3.5 % au bord d'attaque selon les sources ; T3B U∞ = 9.4 m/s,
Tu ≈ 6 à 6.5 % ; T3A- U∞ = 19.8 m/s, Tu ≈ 0.9 % ; ν = 1.5e-5 m²/s. Les valeurs d'entrée
T3B (6.5 %, 100) et T3A- (0.874 %, 8.72) utilisées ici sont celles de Langtry & Menter,
appliquées 0.04 m en amont du bord d'attaque, **sans ajustement** faute de mesures de
décroissance accessibles.

## 5. Sensibilités mesurées (T3A, Re_x du minimum de C_f / de mi-transition / 90 %)

| Variante | Minimum C_f | Mi-transition | 90 % |
|---|---:|---:|---:|
| **Référence** : 13 760 cellules, `linearUpwindLimited`, Tu ajusté | 1.47e5 | 1.91e5 | 2.59e5 |
| maillage fin 27 720 cellules | 1.50e5 | 1.95e5 | 2.56e5 |
| maillage grossier du test (4 100 cellules) | 1.35e5 | 1.86e5 | 2.67e5 |
| `upwind` (13 760 / 27 720 cellules) | 1.64e5 / 1.62e5 | 2.08e5 / 2.06e5 | 2.71e5 / 2.64e5 |
| `linearUpwind` non limité (13 760) | 1.40e5 | 1.86e5 | 2.55e5 |
| Tu d'entrée de Langtry-Menter (3.3 %, 12) | 1.72e5 | 2.17e5 | 2.79e5 |
| production de Kato-Launder (`upwind`) | 1.68e5 (contre 1.64e5) | — | — |
| Expérience (Savill 1993) | ≈ 1.42e5 | ≈ 2.28e5 | ≈ 2.70e5 |

« Mi-transition » et « 90 % » : Re_x où (C_f − C_f,Blasius)/(C_f,turb − C_f,Blasius) atteint
0.5 et 0.9 après le minimum (C_f,turb = 0.0576 Re_x^-0.2) ; pour l'expérience, interpolation
linéaire entre les points mesurés (un tous les 0.1 m, soit ΔRe_x = 3.6e4).

## 6. Ce qui n'a pas marché / reste ouvert

- T3A- et T3B : pas de données numériques accessibles hors ligne ; comparaison par lecture
  graphique seulement.
- La transition trop raide (T3A, T3A-) et le C_f laminaire trop élevé à Tu élevé ne
  s'expliquent pas par le maillage (écarts ≤ 2 % entre 13 760 et 27 720 cellules) ; ils
  viennent du modèle et/ou du réglage amont. Je n'ai pas pu comparer à des résultats
  publiés du modèle γ sur ces cas (article non consulté).
- `linearUpwind` non limité sur les équations de turbulence : non convergé sur T3A- ;
  plateau de résidus ~1.5e-6 (cycle limite) sur T3A.
