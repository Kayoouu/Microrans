# État du projet (mis à jour le 2026-10-04)

## Où on en est

- microrans : outil RANS / URANS 1D, 2D, 3D (périmètre réduit) avec mailleur, interface
  PySide6 et exécutables PyInstaller Windows / Linux (CI GitHub Actions).
- Dernier lot terminé : **audit 2 approfondi** (2026-10-04, demande : « n'hésite pas à
  approfondir l'audit ») : contrôles de cohérence (interface = ligne de commande, reprise
  exacte, 2D = 3D une couche, maillages exportés, CSV pariétaux), fichiers Windows, Ctrl-C,
  cas limites, commandes de la doc, arrêts de l'interface, exécutable construit en local.
  Constats C15 à C23, M20 à M22, U20, P1, P2, L7, D8 dans `docs/audit_utilisateur.md`
  (partie « Audit 2 approfondi ») ; lots F1, F1b, F2 à F4 dans a_faire. Avant : audit 2
  (6608f4d), #30 transsonique (688df2a).
- Exécutables de 688df2a (lot #30) : run https://github.com/Kayoouu/Microrans/actions/runs/37145810874
  (Windows : artefact 11281804767, Linux : 11282345490 ; expirent le 2027-01-01). CI verte
  (tests https://github.com/Kayoouu/Microrans/actions/runs/37145806547).
- Suite de tests : 455 (454 réussis, 1 ignoré), 5 à 8 min en série sur la machine de session.
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

- **Défauts connus non encore corrigés (audits 2, lots F1 et F1b)** : polaire / balayage /
  multigrille en échec hors du dossier du cas (C11) ; interface : canal_turbulent_3d
  calculé faux (C13), maillage périmé utilisé (U14), vecteurs de bruit (U17) ; frontière
  inexistante dans `forces` → calcul perdu (C14) ; vitesse de référence par défaut
  incohérente (C15 : C_d × 9 entre interface et ligne de commande, Re faux si U ≠ 1) ;
  contour CSV « x;y » à virgule décimale → 14 Go (C18) ; balayage d'une clé mal écrite
  accepté (C21) ; Ctrl-C perd tout (C19) ; reprise compressible implicite non exacte
  (C17) ; plantages : perturbation 3D (C16), balayage mesh.nx (C20), polaire 3D (C22),
  multigrille sur extrusion fine (C23) ; **BLAS multi-fil par défaut : deux calculs
  simultanés jusqu'à 13 fois plus lents (P1)**.

- Un calcul = un cœur de solveur, mais BLAS occupe les 4 cœurs sans gain (P1) ; SIMPLE lent sur maillages fins étirés.
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
