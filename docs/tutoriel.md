# Premier calcul

Le premier calcul suit le même chemin dans l'interface et en ligne de commande :
1. calculer la cavité entraînée (Re = 100) ;
2. comparer le résultat à la référence publiée ;
3. changer le nombre de Reynolds ;
4. calculer la traînée d'un cylindre ;
5. vérifier que le résultat ne dépend plus du maillage.

Comptez 15 minutes. Les sorties montrées ici ont été obtenues en suivant ce texte. Les
durées dépendent de la machine.

Termes : [`glossaire.md`](glossaire.md) ; en cas de problème :
[`depannage.md`](depannage.md) ; toutes les clés du fichier de cas :
[`reference_cas.md`](reference_cas.md).

## 1. Le cas : une cavité carrée

Une boîte carrée de côté 1, remplie de fluide. Le couvercle glisse vers la droite à la
vitesse 1 et entraîne le fluide, qui tourne en un grand tourbillon.

Avec ν = 0.01 m²/s, le nombre de Reynolds vaut Re = U L / ν = 1 × 1 / 0.01 = 100 :
l'écoulement est laminaire et stationnaire. Ghia, Ghia & Shin (1982) ont publié la vitesse
le long de l'axe vertical : c'est la référence classique pour vérifier un code.

## 2. Dans l'interface

Lancer `microrans gui`, ou l'exécutable `microrans-gui` (`microrans-gui.exe` sous Windows).

1. **Accueil.** Choisir, dans le groupe « Commencer ici », « Cavité entraînée, Re = 100 »
   (colonne « Durée » : 6 s). Sa description s'affiche sous la liste. Cliquer « Ouvrir
   l'exemple ».
2. **2. Maillage.** L'interface y passe d'elle-même. Cliquer « Générer le maillage » :
   ```
   4096 cellules (4096 quadrilatères)
   non-orthogonalité max 0.0°, moyenne 0.0°
   asymétrie max 0.00, allongement max 1
   ```
   C'est un maillage régulier de 64 × 64 carrés, donc de qualité parfaite.
3. **3. Physique** et **4. Conditions limites** : rien à changer.
   - Page 3 : ν = 0.01 et modèle laminaire.
   - Page 4 : le couvercle « lid » est une paroi qui se déplace à U = (1, 0) ; les trois
     autres côtés (« walls ») sont des parois fixes.
4. **6. Calcul.** Cliquer « Lancer le calcul ». Les résidus s'affichent pendant le calcul
   et doivent baisser.
5. **7. Résultats.** L'interface y passe à la fin. Le résumé commence par :
   ```
   Calcul stationnaire, laminaire — 4 096 cellules.
   Convergé en 351 itérations.
   ν = 0.01 m²/s, U = 1, L = 1 → Re = U L / ν = 100.
   ```
   Choisir un champ dans la liste (|U|, U_x, pression, vorticité…) et cliquer « Tracer
   le champ ». « Ouvrir le dossier de résultats » montre les fichiers écrits.
6. **Changer le Reynolds.**
   - Page 3 : ν = 0.0025 (soit Re = 400).
   - Page 6 : « Lancer le calcul ».
   - Le résumé affiche `Re = U L / ν = 400` et « Convergé en 312 itérations ».
7. **Garder son cas** : menu Fichier, « Enregistrer sous… ». Le fichier `.toml` enregistré
   est le même que celui de la ligne de commande (§ 3).

Les sections 4 et 5 se font aussi dans l'interface :
- **Profil sur une ligne** : page « 7. Résultats », « Tracer le profil ».
- **Cylindre** : Accueil, « Cylindre, Re = 20 ».
- **Raffiner le maillage** : page « 2. Maillage », nombres de mailles.

## 3. En ligne de commande

```
microrans examples                    # liste des exemples, avec leur durée
microrans examples cavite_re100       # copie modifiable dans le dossier courant
```
```
Copié : cavite_re100.toml
Modifier cavite_re100.toml, puis : microrans run2d cavite_re100.toml
```

Le fichier de cas est un texte en sections :
- `[mesh]` : rectangle 64 × 64, et noms des quatre côtés ;
- `[physics]` : `nu = 0.01` ;
- `[boundary.lid]` : paroi mobile, `U = [1.0, 0.0]` ;
- `[boundary.walls]` : parois fixes ;
- `[solver]` : `max_iter = 2000`, `tol = 1e-6` ;
- `[output]` : dossier de résultats.

```
microrans run2d cavite_re100.toml
```
```
Cas 2D : 4096 cellules, modèle Laminaire (ν_t = 0), ν = 0.01, stationnaire
Maillage : non-orthogonalité max 0.0° (moy. 0.0°), asymétrie max 0.00
  it     1  Ux=1.00e+00  Uy=0.00e+00  p=1.00e+00  continuity=1.94e-05
  it   100  Ux=4.76e-04  Uy=7.49e-04  p=5.58e-04  continuity=5.05e-09
  …
  convergé en 351 itérations (3.5 s)

Calcul stationnaire, laminaire — 4 096 cellules.
Convergé en 351 itérations (3.5 s).
ν = 0.01 m²/s, U = 1, L = 1 → Re = U L / ν = 100.

Efforts (coefficients : force / (½ U² L), par unité de profondeur) :
  frontière        Cd        Cl        Cm   y⁺ max
  walls        0.3799  -0.06326    0.4984    0.497
  lid         -0.3835   0.06281   -0.5389   0.7646
…
Résultats dans …/results/cavite_re100
```

Lire cette sortie :
- **Les colonnes Ux, Uy, p, continuity** sont les résidus : de combien chaque équation est
  encore loin d'être satisfaite. Le calcul s'arrête quand ils passent sous `tol`
  (« convergé »).
- **Les efforts** : sur une cavité fermée, Cd et Cl n'ont pas de sens aérodynamique. Ce
  sont les efforts du fluide sur chaque frontière, adimensionnés.

Fichiers écrits dans `results/cavite_re100` :

| Fichier | Contenu |
|---|---|
| `U.png`, `p.png`, `vorticity.png`, `mesh.png` | figures des champs et du maillage |
| `convergence.png`, `history.csv` | résidus au fil des itérations |
| `summary.json` | toutes les valeurs du résumé |
| `wall_lid.csv`, `wall_walls.csv` | le long de chaque paroi : Cp, Cf, y⁺ |
| `fields.vtk` | champs complets, à ouvrir dans ParaView |
| `checkpoint.npz` | sauvegarde, pour poursuivre le calcul |

## 4. Comparer à la référence

Ajouter à la fin de `cavite_re100.toml` un profil le long de l'axe vertical x = 0.5 :

```toml
[[output.lines]]
name = "axe"
start = [0.5, 0.0]
end = [0.5, 1.0]
n = 101
```

Relancer `microrans run2d cavite_re100.toml`. Le fichier `line_axe.csv` (colonnes
`s, x, y, Ux, Uy, U_mag, p`) et la figure `line_axe.png` apparaissent. Comparaison avec le
tableau de Ghia et al. :

| y | u de Ghia et al. (1982) | u de microrans (64 × 64) |
|---|---|---|
| 0.8516 | +0.2315 | +0.2357 |
| 0.5000 | −0.2058 | −0.2086 |
| 0.4531 | −0.2109 | −0.2132 |
| 0.1719 | −0.1015 | −0.1016 |

Sur les 17 points du tableau, l'écart maximal est de 0.004. Le tableau complet est dans
`microrans/fv2d/benchmarks.py`.

## 5. Changer un réglage sans modifier le fichier

```
microrans run2d cavite_re100.toml --set physics.nu=0.0025 -o results/cavite_re400
```
```
Convergé en 312 itérations (2.6 s).
ν = 0.0025 m²/s, U = 1, L = 1 → Re = U L / ν = 400.
```

`--set` remplace une clé pour ce calcul seulement. Plusieurs clés se suivent ; un texte
s'écrit tel quel : `--set solver.convection_U=upwind solver.max_iter=400`.

**Faute de frappe.** Une clé mal écrite est signalée, puis **ignorée** :
```
ATTENTION : [solver] max_iterr : clé inconnue, ignorée — vouliez-vous dire « max_iter » ?
```

## 6. Un résultat physique : la traînée d'un cylindre

```
microrans examples cylindre_re20
microrans run2d cylindre_re20.toml
```
```
Convergé en 182 itérations (2.6 s).
ν = 0.05 m²/s, U = 1, L = 1 → Re = U L / ν = 20.

  frontière        Cd        Cl        Cm   y⁺ max
  cylinder      2.037         0         0  0.06368
  cylinder : Cd = 1.22 (pression) + 0.8172 (frottement).
```

Les grandeurs de référence sont U = 1 et le diamètre L = 1, donc Cd = traînée / (½ U² L).
La référence de Dennis & Chang (1970) donne Cd ≈ 2.05 : l'écart est de 0.6 %.

**Le maillage est-il assez fin ?** Refaire le calcul avec deux fois plus de mailles dans
chaque direction et une première maille deux fois plus fine :
```
microrans run2d cylindre_re20.toml --set mesh.n_around=192 mesh.n_radial=128 mesh.first_height=0.005 -o results/cyl_fin
```
```
Cas 2D : 24576 cellules, …
Convergé en 545 itérations (32.2 s).
  cylinder      2.033 …
```

Cd passe de 2.037 à 2.033, soit 0.2 % d'écart avec 4 fois plus de cellules. Le premier
maillage suffisait. Ce contrôle est à faire pour tout nouveau cas.

## Ensuite

- **Autres exemples** (`microrans examples`) : turbulence (plaque plane, profil NACA 0012,
  polaire), thermique, fluides non newtoniens, compressible. Chaque fichier d'exemple
  commence par sa référence.
- **Instationnaire** : `cylindre_re100_urans` (lâcher de tourbillons, 2,5 min). À Re = 100,
  un calcul stationnaire du cylindre converge vers une solution fausse
  (voir [`depannage.md`](depannage.md), § 4).
- **Turbulence** : vérifier y⁺ (colonne `y⁺ max` du résumé) : ≈ 1 en traitement résolu,
  30 à 300 avec lois de paroi.
- **Un calcul qui diverge ou ne converge pas** : voir [`depannage.md`](depannage.md).
