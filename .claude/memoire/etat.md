# État du projet (mis à jour le 2026-10-15)

## Où on en est

- microrans : outil RANS / URANS 1D, 2D, 3D (périmètre réduit) avec mailleur, interface
  PySide6 et exécutables PyInstaller Windows / Linux (CI GitHub Actions).
- Dernier lot : **G1 (C24), partie 2** (2026-10-15) : chiffres du README remesurés après
  la correction (domaine et O / C du NACA, polaire, GCI, de Vahl Davis, Ghia) ; polaire :
  C_d +0.4 % à +3.6 %, mais +24 % d'itérations (×2 à ×3 de 10° à 14°) : plateau de
  résidus dû au limiteur non orthogonal → tâche suivante. Partie 1 (2026-10-14) : flux de
  Rhie-Chow du stationnaire sur la diagonale physique → SIMPLE / SIMPLEC = couplé à 1e-8,
  quels que soient relax_U / pseudo_cfl ; NACA C_d +0.89 %, C_l −0.24 %. Avant : **audit 3, partie 1**
  (2026-10-13) : campagnes sans régression, c17 (21 exemples honnêtes), C24, figures L8 à
  L11. Avant : maillage en C (10 au 12 octobre), GCI (9), F4 à F1, audit 2 : `journal.md`.
- Exécutables de 5b0d45d (G1 partie 2) : run https://github.com/Kayoouu/Microrans/actions/runs/38683284836
  (Windows : artefact 11949073059, Linux : 11949285978 ; expirent le 2027-01-13). CI verte
  (tests https://github.com/Kayoouu/Microrans/actions/runs/38683278975).
- Suite de tests : 514 (512 réussis, 2 ignorés), 5 à 10 min en série sur la machine de session.
- Objectif accepté le 2026-10-05 : jalon A (2D académique crédible), puis B (3D
  académique, petite géométrie) ; C (pré-industriel) non visé. Détail : `a_faire.md`.
- Méthode outillée (2026-10-05) : skill `lot`, `tools/dev/`, hooks avant poussée et fin de
  tour, CI au démarrage, `bilan_lots.md` (estimé / réel, rétrospectives). Commit b81d766, CI
  verte : https://github.com/Kayoouu/Microrans/actions/runs/37388788153
- Aucune release publiée depuis cet environnement (tag refusé) : c'est à l'utilisateur.
- Autonomie : routine quotidienne `trig_01RCTjquBZw82GuKW83rPuSG` (cron `53 3 * * *` UTC,
  dans cette session), 1 lot par jour sauf exception ; première relance 2026-10-03.

## Chiffres mesurés de référence (machine virtuelle 4 cœurs, jusqu'à ×1.7 selon le jour)

| Cas | Mesure |
|---|---|
| 3D, 10⁶ hexaèdres, laminaire, exe Linux | 13.8 s / itération, pic 3.09 Go, 461 s pour maillage + distance paroi + 10 it + VTK |
| Même cas en Python, côte à côte | 13.6 s / it, 3.18 Go, 456 s → **exe = Python** |
| Même cas, mesure antérieure (machine plus rapide ce jour-là) | 8.5 s / it, 2.9 Go, préparation 31 s + 172 s |
| 3D 64³ (262 144) dans l'interface | maillage + vue 11.6 s, 3 it en 65.5 s (solveur 10 s), figures < 1 s, pic 1.19 Go |
| fields.vtk 10⁶ cellules (E3, côte à côte) | texte 14.0–15.1 s, 173 Mo → binaire 0.31–0.35 s, 121 Mo |
| Distance à la paroi 10⁶ cellules (E2, côte à côte) | 206 s → 25.5 s, identique au bit près |
| Calcul complet 10⁶ cellules, 10 it (E2, côte à côte) | 375 s → 195 s (préparation 254 s → 78 s) |
| Distance à la paroi 2D, 490 000 cellules | 91 s → 2.7 s |
| NACA 0012 Euler M = 0.8, exemple (venkat_k 0.3), résidus < 1e-8 | 192 × 64 : 452 it, 50 s ; 384 × 128 : 868 it, 409 s (avant, K 0.05 : 1 011 it / 112 s ; ~13 min) |
| Même cas, C_l ; C_d (384 × 128, bord de fuite fermé) | 0.3333 ; 0.02178 (référence hp-DG lue dans un résumé : 0.333 ; 0.02135) |

Préparation 3D maintenant dominée par la construction du maillage (~35 s à 10⁶) ; un
calcul stationnaire 10⁶ cellules = quelques centaines d'itérations ≈ 1 à 1.5 h.

## Décisions prises (et pourquoi)

- Pas de changement de langage : l'exécutable va aussi vite que Python (mesuré) ; le coût
  est dans les algorithmes (distance à la paroi, SIMPLE), pas dans l'emballage.
- Distance à la paroi : exacte, k-d (8 voisins, `workers=-1`) + arbre de boîtes
  (`mesh2d/bvh.py`) + borne par face ; une descente « vers la boîte la plus proche » sans
  k-d a été essayée et rejetée (64³ : 228 s, vecteurs changés aux égalités).
- Numba désactivé par défaut (`[solver] numba = false`) et exclu de l'exécutable : ~10 % de
  gain mesuré, multi-fil plus lent sur la machine de test.
- Vitesse de référence (C15, 2026-10-04) : `fv2d.solver.choose_reference_velocity`, une
  seule règle pour ν = U L / Re, coefficients, γ̇_ref, C_T : `reference_velocity` donnée,
  sinon vitesse d'entrée (débitante, la plus grande + avertissement), sinon paroi mobile,
  sinon 1 ; jamais la vitesse initiale ; interface : champ vide = auto.
- BLAS : 1 fil par calcul, fixé dans `microrans/__init__.py` avant l'import de NumPy
  (F1 ; 6 cas mesurés A/B/B/A, aucun plus lent ; 2 calculs simultanés 9 s au lieu de
  137 s). La mesure « Numba multi-fil plus lent » datait d'avant : à refaire.
- 3D : figures en coupe plane — interface : x, y ou z = cte (`mesh3d/slice.py` :
  `PlaneSlice` polygones d'intersection, tous types de cellules, 0.4 s à 10⁶ ; `ZSlice`
  inchangé pour z, maillages en couches) ; ligne de commande : `[output] slice_axis` /
  `slice_value` (F3, défaut plan z médian) ; champs
  complets dans `fields.vtk` (ParaView). Compressible, couplé, poreux, swirl, disques,
  animations, axisymétrique : refusés en 3D avec message.
- VTK : binaire par défaut (`[output] vtk_format = "ascii"` = ancien texte, identique
  octet pour octet) ; `read_vtk` dans `mesh2d/io.py` (lecteur minimal, tests) ; pas de
  dépendance à VTK/meshio (test avec le lecteur VTK officiel ignoré s'il est absent).
- Compressible : seuil du limiteur `venkat_k` = 0.3 dans l'exemple NACA (0.05 = défaut du
  solveur, comme SU2 : deux solutions stationnaires et C_l −5 % sur 96 × 32) ; défaut non
  changé (rampe M = 2, Sod non remesurés). NACA : `trailing_edge = closed | open | sharp`.
- Rhie-Chow en stationnaire (C24, 2026-10-14) : coefficient V / a_P physique (sans
  relax_U ni pseudo-temps) : correction de Majumdar pour relax_U (exacte, sans gradient
  explicite), terme explicite pour le pseudo-temps (`_pressure_correction(steady_ref=…)`) ;
  plateau de résidus plus haut sur maillage non orthogonal (limiteur, voir a_faire) ;
  pas d'option pour revenir à l'ancien comportement (résultat dépendant des réglages) ;
  instationnaire inchangé (ddtCorr).
- Extrusion multi-couches ≠ 2D exact (Rhie-Chow : diffusion z dans a_P) ; une couche entre
  deux plans de symétrie = 2D exact. Documenté (README § 8, limite 17).

## Limites ouvertes importantes (détail : README § 8)

- **Défauts connus non encore corrigés (audit 3)** : figures L8 à L11 (lot G2). C24 corrigé
  le 2026-10-14 (chiffres publiés à remesurer, G1 partie 2).

- Un calcul = un cœur (BLAS limité à 1 fil depuis F1) ; SIMPLE lent sur maillages fins étirés.
- NACA 0012 transsonique : ancien « −4.5 % » expliqué (référence d'une autre géométrie +
  seuil du limiteur) ; restent C_d +2 % et une valeur citée de Vassberg & Jameson (≈ 0.347,
  bord de fuite pointu) 3.7 % au-dessus, non vérifiable ici (articles bloqués par le proxy).
- Transition γ validée seulement sur plaques sans gradient de pression.
- 3D : pas de mailleur général ni d'import 3D.

## Fichiers clés

- Solveur incompressible 2D/3D : `microrans/fv2d/solver.py` ; cas : `fv2d/case.py` ;
  vérification des cas : `fv2d/validate.py` ; compressible : `fv2d/compressible*.py`.
- Maillage 2D : `microrans/mesh2d/` ; 3D : `microrans/mesh3d/` (`mesh.py` dont distance à
  la paroi, `generators.py` pavé et extrusion, `slice.py` coupes planes).
- Interface : `microrans/gui/app.py` (auto-test `_selftest`, lancé par la CI sur les exe),
  `gui/widgets.py`. Exemples : `microrans/examples/*.toml`, catalogue `microrans/catalog.py`.
- CI : `.github/workflows/tests.yml`, `build.yml` ; spec PyInstaller dans `packaging/`.
