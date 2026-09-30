# Calcul sur plusieurs cœurs

Toutes les valeurs sont MESURÉES (machine de développement : 4 cœurs logiques Intel Xeon
2.1 GHz virtualisés, Linux, Python 3.11, NumPy 2.4, SciPy 1.17), machine sans autre charge.

## 1. Balayages et polaires en parallèle (livré)

### Utilisation

```bash
microrans polar naca0012_polaire --alpha -4 14 2 --jobs 4     # 4 points à la fois
microrans sweep cavite_re100 --param physics.nu --values 0.005 0.01 0.02 0.04 -j 4
microrans polar naca0012_polaire --alpha -4 14 2 -j 0         # 0 = tous les cœurs
```

Dans le fichier de cas : `[sweep] jobs = 4`. Interface : page « 6. Calcul », cadre
« Polaire / balayage », champ **« Calculs en parallèle »**. Défaut : 1 (calcul séquentiel,
comportement inchangé).

### Fonctionnement

- Chaque point du balayage est un calcul complet, exécuté dans un **processus séparé**
  (`concurrent.futures.ProcessPoolExecutor`, démarrage « spawn » : identique sous Windows,
  Linux, macOS et dans l'exécutable PyInstaller). Le maillage est construit une fois et
  transmis aux processus.
- Bibliothèques multi-fils (OpenBLAS, OpenMP, MKL, Numba) limitées à **1 fil par
  processus** pendant le balayage (variables d'environnement héritées) : pas de
  sursouscription des cœurs.
- **Sans continuation** (`--no-continuation`) : les points sont indépendants et distribués
  dynamiquement (un processus libre prend le point suivant). Résultats **identiques au bit
  près** au calcul séquentiel (testé).
- **Avec continuation** (défaut) : les valeurs sont découpées en N blocs **contigus** (N =
  nombre de processus), chacun parcouru dans l'ordre avec continuation. Seul le premier
  point de chaque bloc part de l'état initial. Pourquoi des blocs contigus plutôt qu'un
  point sur N : l'écart entre deux points successifs d'une chaîne reste le pas du balayage,
  donc le point de départ est aussi bon qu'en séquentiel. Conséquence : les points convergés
  sont égaux au calcul séquentiel **à la tolérance de convergence près**, pas au bit près
  (écarts mesurés ci-dessous) ; au-delà du décrochage (solutions multiples, hystérésis), un
  point en tête de bloc peut tomber sur une autre branche que le balayage séquentiel.
- Les points arrivent dans le désordre (affichage de l'interface mis à jour au fil de
  l'eau) ; le tableau final (`balayage.csv`, `balayage.json`) est dans l'ordre des valeurs.
- Journal de chaque point : `journal.txt` dans son sous-dossier (les sorties de plusieurs
  calculs simultanés seraient illisibles dans une seule console). La console affiche une
  ligne par point terminé.
- Bouton « Arrêter » (interface) : plus aucun point n'est lancé, les points en cours
  s'arrêtent à l'itération suivante (marqués non convergés).

### Gains mesurés

| Balayage | Points | 1 processus | 2 processus | 4 processus | Résultats |
|---|---:|---:|---:|---:|---|
| Cavité Re = 25 à 250 (ν), 64 × 64, sans continuation | 8 | 30.1 s | 15.8 s (×1.90) | 8.7 s (×3.46) | identiques au bit près |
| Polaire NACA 0012 SA, −4° à 14°, 8 192 cellules, avec continuation | 10 | 227 s | — | 88 s (×2.58) | C_l à 0.1 %, C_d à 0.45 % près (10° et 12°) |

Polaire : gain inférieur à 4 parce que chaque bloc repart « à froid » (itérations totales
6 465 → 7 176, +11 %) et que les blocs (3, 3, 2, 2 points) ne durent pas le même temps. Les
écarts de C_l / C_d viennent du point de départ différent de la tête de chaque bloc : ils
restent dans la marge laissée par le critère d'arrêt (efforts stabilisés).

**Script Python** : `run_sweep(..., jobs=N)` lance des processus en mode « spawn » ; le
script appelant doit protéger son code par `if __name__ == "__main__":` (règle de
`multiprocessing`), sinon chaque processus relance le script entier et le balayage ne se
termine jamais. La ligne de commande, l'interface et les exécutables sont déjà protégés.

### Limites

- Mémoire : un maillage et un solveur par processus (≈ 140 Mo par processus pour 8 192
  cellules, mesuré) ; 4 processus ≈ 4 fois la mémoire d'un calcul.
- Démarrage : chaque processus réimporte Python, NumPy, SciPy (≈ 0.4 s ici ; plus sous
  Windows et depuis l'exécutable) : inutile pour des points de quelques secondes.
- Gain plafonné par le nombre de points et leur équilibre : 10 points sur 4 processus avec
  continuation = blocs de 3, 3, 2, 2 points ; un bloc qui contient un point lent (près du
  décrochage) retarde la fin. Sans continuation, la répartition est dynamique.
- Un portable a rarement 4 « vrais » cœurs à pleine fréquence : la fréquence turbo baisse
  quand tous les cœurs travaillent et la bande passante mémoire est partagée. Le gain par
  cœur supplémentaire est donc inférieur à 1 (voir le débit de la machine ci-dessous).
- Le rappel `callback(solveur, n)` de `run_sweep` n'est appelé qu'en séquentiel (les
  solveurs sont dans d'autres processus) ; l'arrêt passe par `should_stop()`.
- Exécutables : l'intégration continue (build.yml) lance un balayage `--jobs 2` et
  l'auto-test de l'interface (balayage sur 2 processus) avec les exécutables Windows et
  Linux construits par PyInstaller : réussis. `multiprocessing.freeze_support()` est appelé en tête des deux
  points d'entrée (`packaging/*_entry.py`), ce qu'exige PyInstaller pour les processus
  « spawn ».

## 2. Un seul calcul sur plusieurs cœurs : NON fait

Accélérer UN calcul (opérateurs multi-fils, par ex. Numba) n'a pas été étudié : la tâche a
été interrompue avant. Le profilage (NACA 0012, 8 192 cellules) montre un temps réparti
sur de nombreuses opérations NumPy (solveurs linéaires ≈ 30 % seulement) : le gain attendu
d'une parallélisation partielle est modeste (Amdahl) ; c'est un chantier à part.
