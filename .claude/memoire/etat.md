# État du projet (mis à jour le 2026-10-02)

## Où on en est

- microrans : outil RANS / URANS 1D, 2D, 3D (périmètre réduit) avec mailleur, interface
  PySide6 et exécutables PyInstaller Windows / Linux (CI GitHub Actions).
- Dernier lot terminé : **Lot E — 3D dans l'exécutable** (interface et ligne de commande),
  commits f08fb62 et b289a91. CI verte (tests + exécutables).
- Exécutables de b289a91 : run https://github.com/Kayoouu/Microrans/actions/runs/36941858951
  (Windows : artefact 11201355109, Linux : 11201055273 ; expirent le 2026-12-30).
- Suite de tests : 436 réussis, 1 ignoré (437), ~8 min en série sur la machine de session.
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
| fields.vtk ASCII 10⁶ cellules | 169 Mo |

La préparation (distance à la paroi exacte) domine le temps des petits nombres
d'itérations en 3D ; un calcul stationnaire 10⁶ cellules = quelques centaines d'itérations
≈ 1 à 1.5 h.

## Décisions prises (et pourquoi)

- Pas de changement de langage : l'exécutable va aussi vite que Python (mesuré) ; le coût
  est dans les algorithmes (distance à la paroi, SIMPLE), pas dans l'emballage.
- Numba désactivé par défaut (`[solver] numba = false`) et exclu de l'exécutable : ~10 % de
  gain mesuré, multi-fil plus lent sur la machine de test.
- 3D : figures en coupe z = constante (`mesh3d/slice.py`, maillages en couches selon z) ;
  champs complets dans `fields.vtk` (ParaView). Compressible, couplé, poreux, swirl, disques,
  animations, axisymétrique : refusés en 3D avec message.
- Extrusion multi-couches ≠ 2D exact (Rhie-Chow : diffusion z dans a_P) ; une couche entre
  deux plans de symétrie = 2D exact. Documenté (README § 8, limite 17).

## Limites ouvertes importantes (détail : README § 8)

- Un calcul = un cœur ; SIMPLE lent sur maillages fins étirés.
- C_l NACA 0012 transsonique ~4.5 % bas, non expliqué (tâche #30).
- Transition γ validée seulement sur plaques sans gradient de pression.
- 3D : pas de mailleur général ni d'import 3D ; coupes x / y absentes de l'interface.

## Fichiers clés

- Solveur incompressible 2D/3D : `microrans/fv2d/solver.py` ; cas : `fv2d/case.py` ;
  vérification des cas : `fv2d/validate.py` ; compressible : `fv2d/compressible*.py`.
- Maillage 2D : `microrans/mesh2d/` ; 3D : `microrans/mesh3d/` (`mesh.py` dont distance à
  la paroi, `generators.py` pavé et extrusion, `slice.py` coupe z).
- Interface : `microrans/gui/app.py` (auto-test `_selftest`, lancé par la CI sur les exe),
  `gui/widgets.py`. Exemples : `microrans/examples/*.toml`, catalogue `microrans/catalog.py`.
- CI : `.github/workflows/tests.yml`, `build.yml` ; spec PyInstaller dans `packaging/`.
