# État du projet (mis à jour le 2026-10-05)

## Où on en est

- microrans : outil RANS / URANS 1D, 2D, 3D (périmètre réduit) avec mailleur, interface
  PySide6 et exécutables PyInstaller Windows / Linux (CI GitHub Actions).
- Dernier lot terminé : **F2 — messages** (2026-10-06, relance quotidienne) : M13 à M22,
  U16, L6 ; 12 tests (tous « échoue avant : OK ») ; résultats de calcul inchangés au bit
  près (ab.py, 5 exemples). Avant : **F1b — plantages et reprises** (2026-10-05) :
  C16 (perturbation 3D), C17 (reprise compressible implicite exacte au bit près), C19
  (Ctrl-C écrit tout, code 130), C22 (polaire 3D), C23 (multigrille sur extrusion fine) ;
  campagne 7 : 18/18 exemples repris au bit près. Avant : **F1** (2026-10-04) : P1 (1 fil BLAS), C11, C13, C14, C18, C20, C21, U14, U17 corrigés, un test
  par point (chacun échoue sur l'ancien code). Puis **C15** (l'utilisateur m'a laissé
  choisir « la meilleure solution ») : une seule vitesse de référence. Avant : audit 2
  approfondi (576a0bc, e5deb38).
- Exécutables de 898cffa (F1b) : run https://github.com/Kayoouu/Microrans/actions/runs/37263377558
  (Windows : artefact 11325098571, Linux : 11324344878 ; expirent le 2027-01-03). CI verte
  (tests https://github.com/Kayoouu/Microrans/actions/runs/37263373865). Exécutable F1
  construit en local : 1 fil BLAS effectif (plaque : CPU 9.4 s pour 9.6 s).
- Suite de tests : 491 (489 réussis, 2 ignorés), 5 à 8 min en série sur la machine de session.
- Objectif accepté le 2026-10-05 : jalon A (2D académique crédible), puis B (3D
  académique, petite géométrie) ; C (pré-industriel) non visé. Détail : `a_faire.md`.
- Méthode outillée (2026-10-05) : skill `lot`, `tools/dev/`, hooks avant poussée et fin de
  tour, CI au démarrage, `bilan_lots.md` (estimé / réel, rétrospectives). Commit b81d766, CI
  verte : https://github.com/Kayoouu/Microrans/actions/runs/37388788153
- Aucune release publiée depuis cet environnement (tag refusé) : c'est à l'utilisateur.
- Autonomie : routine quotidienne `trig_01RCTjquBZw82GuKW83rPuSG` (cron `53 3 * * *` UTC,
  dans cette session), 1 lot par jour sauf exception ; première relance 2026-10-03.

## Chiffres mesurés de référence (machine virtuelle 4 cœurs, ±50 % selon le jour)

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
  inchangé pour z, maillages en couches) ; ligne de commande : plan z médian ; champs
  complets dans `fields.vtk` (ParaView). Compressible, couplé, poreux, swirl, disques,
  animations, axisymétrique : refusés en 3D avec message.
- VTK : binaire par défaut (`[output] vtk_format = "ascii"` = ancien texte, identique
  octet pour octet) ; `read_vtk` dans `mesh2d/io.py` (lecteur minimal, tests) ; pas de
  dépendance à VTK/meshio (test avec le lecteur VTK officiel ignoré s'il est absent).
- Compressible : seuil du limiteur `venkat_k` = 0.3 dans l'exemple NACA (0.05 = défaut du
  solveur, comme SU2 : deux solutions stationnaires et C_l −5 % sur 96 × 32) ; défaut non
  changé (rampe M = 2, Sod non remesurés). NACA : `trailing_edge = closed | open | sharp`.
- Extrusion multi-couches ≠ 2D exact (Rhie-Chow : diffusion z dans a_P) ; une couche entre
  deux plans de symétrie = 2D exact. Documenté (README § 8, limite 17).

## Limites ouvertes importantes (détail : README § 8)

- **Défauts connus non encore corrigés (lot F1b)** : Ctrl-C perd tout (C19) ; reprise compressible
  implicite non exacte (C17) ; plantages : perturbation 3D (C16), polaire 3D (C22),
  multigrille sur extrusion fine (C23).

- Un calcul = un cœur (BLAS limité à 1 fil depuis F1) ; SIMPLE lent sur maillages fins étirés.
- NACA 0012 transsonique : ancien « −4.5 % » expliqué (référence d'une autre géométrie +
  seuil du limiteur) ; restent C_d +2 % et une valeur citée de Vassberg & Jameson (≈ 0.347,
  bord de fuite pointu) 3.7 % au-dessus, non vérifiable ici (articles bloqués par le proxy).
- Transition γ validée seulement sur plaques sans gradient de pression.
- 3D : pas de mailleur général ni d'import 3D ; plan des figures non réglable en ligne de
  commande (z médian).

## Fichiers clés

- Solveur incompressible 2D/3D : `microrans/fv2d/solver.py` ; cas : `fv2d/case.py` ;
  vérification des cas : `fv2d/validate.py` ; compressible : `fv2d/compressible*.py`.
- Maillage 2D : `microrans/mesh2d/` ; 3D : `microrans/mesh3d/` (`mesh.py` dont distance à
  la paroi, `generators.py` pavé et extrusion, `slice.py` coupes planes).
- Interface : `microrans/gui/app.py` (auto-test `_selftest`, lancé par la CI sur les exe),
  `gui/widgets.py`. Exemples : `microrans/examples/*.toml`, catalogue `microrans/catalog.py`.
- CI : `.github/workflows/tests.yml`, `build.yml` ; spec PyInstaller dans `packaging/`.
