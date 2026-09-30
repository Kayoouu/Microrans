# Calcul sur plusieurs cœurs : notes (à intégrer au README)

Texte proposé pour le README (sections 3 « Ligne de commande », 5 « Méthodes numériques » et
7 « Performances mesurées »). Toutes les valeurs sont MESURÉES (machine de développement :
4 cœurs logiques Intel Xeon 2.1 GHz virtualisés, Linux, Python 3.11, NumPy 2.4, SciPy 1.17,
pyamg 5.3, Numba 0.67). Voir « Conditions de mesure » à la fin.

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

MESURES_PARTIE_1

### Limites

- Mémoire : un maillage et un solveur par processus (≈ MEMOIRE Mo par processus pour
  8 192 cellules) ; 4 processus ≈ 4 fois la mémoire d'un calcul.
- Démarrage : chaque processus réimporte Python, NumPy, SciPy (≈ DEMARRAGE s, plus sous
  Windows et depuis l'exécutable) : inutile pour des points de moins de quelques secondes.
- Gain plafonné par le nombre de points et leur équilibre : 10 points sur 4 processus avec
  continuation = blocs de 3, 3, 2, 2 points ; un bloc qui contient un point lent (près du
  décrochage) retarde la fin. Sans continuation, la répartition est dynamique.
- Un portable a rarement 4 « vrais » cœurs à pleine fréquence : la fréquence turbo baisse
  quand tous les cœurs travaillent et la bande passante mémoire est partagée. Le gain par
  cœur supplémentaire est donc inférieur à 1 (voir le débit de la machine ci-dessous).
- Le rappel `callback(solveur, n)` de `run_sweep` n'est appelé qu'en séquentiel (les
  solveurs sont dans d'autres processus) ; l'arrêt passe par `should_stop()`.
- Exécutable : vérifié sous Linux (ligne de commande `--jobs 3` et auto-test de l'interface
  avec 2 processus, exécutable PyInstaller construit avec `packaging/microrans.spec`). Sous
  Windows, non exécuté sur ce poste : l'intégration continue (build.yml) lance maintenant
  l'auto-test de l'interface (balayage sur 2 processus) et un balayage `--jobs 2` avec
  l'exécutable Windows. `multiprocessing.freeze_support()` est appelé en tête des deux
  points d'entrée (`packaging/*_entry.py`), ce qu'exige PyInstaller pour les processus
  « spawn ».

## 2. Opérateurs multi-fils pour UN calcul (étudié, NON livré)

MESURES_PARTIE_2

## Conditions de mesure

CONDITIONS
