# Dépannage

Que faire quand un cas est refusé, quand le calcul diverge, ne converge pas, ou donne un
résultat douteux. Les messages cités sont ceux qu'affiche le logiciel. Chaque piste
indique le cas sur lequel elle a été essayée et ce qui a été mesuré. Une piste qui n'a pas
marché sur ce cas est signalée aussi.

Termes : [`glossaire.md`](glossaire.md) ; toutes les clés du fichier de cas :
[`reference_cas.md`](reference_cas.md) ; premier calcul : [`tutoriel.md`](tutoriel.md).

**Où lire les messages.**
- En ligne de commande : `Erreur : …` (calcul refusé ou arrêté) et `ATTENTION : …` (calcul
  lancé, point à vérifier).
- Dans l'interface :
  - fenêtre « Réglages à corriger » : le calcul n'est pas lancé ;
  - fenêtre « À vérifier avant de lancer » : on peut lancer quand même ;
  - fin du calcul : résumé de la page « 7. Résultats ».
- Les noms de réglages des messages (`[solver] max_iter`…) sont ceux du fichier de cas,
  visibles dans l'onglet « Fichier de cas (TOML) ».

## 1. Le cas est refusé avant le calcul

Le cas est vérifié avant le maillage. Toutes les erreurs sont listées ensemble. Rien n'est
calculé tant qu'il en reste une.

| Message (extrait) | Cause | Que faire |
|---|---|---|
| `max_iterr : clé inconnue, ignorée — vouliez-vous dire « max_iter » ?` | faute de frappe : la valeur est **ignorée**, la valeur par défaut s'applique | corriger le nom |
| `nombre entier attendu (écrire 10, sans « .0 »)` | `nx = 10.0` | écrire `10` |
| `[physics] : donner la viscosité nu (m²/s) ou le nombre de Reynolds reynolds` | ni ν ni Re | donner l'un des deux |
| `nu = … et reynolds = … donnés ensemble` | les deux | n'en garder qu'un |
| `Ce fichier décrit seulement un maillage (ni [physics] ni [boundary])` | exemples `mesh_*` | `microrans mesh <fichier>`, ou partir d'un exemple de calcul (`microrans examples`) |
| `[[bodies]] « a » et « b » se recouvrent ou se touchent : impossible en maillage hybride` | corps superposés, maillage `hybrid` | écarter les corps, ou maillage `unstructured` |
| `… se recouvrent ou se touchent : ils formeront un seul obstacle` (ATTENTION) | même chose en `unstructured` | voulu (obstacle composé) ou corps à écarter |
| `écart … < 2 × épaisseur des couches de paroi` | couches de paroi qui se croisent (`hybrid`) | écarter les corps, ou réduire `[mesh.layers]` |
| `le maillage en O entoure exactement un corps` | plusieurs corps en `ogrid` | maillage `unstructured` ou `hybrid` |
| `lois de paroi incompatibles avec le modèle ke` (ou `sst_gamma`) | ces modèles exigent y⁺ ≈ 1 | `wall_treatment = "resolved"`, ou SA, k-ω, SST |
| `non newtonien en laminaire uniquement` | loi de viscosité + modèle de turbulence | `model = "laminar"` |
| `[mesh] 1 000 000 cellules : prévoir ~1.2 Go de mémoire et ~15 s par itération` (ATTENTION) | maillage structuré de plus de 500 000 cellules | régler d'abord le cas sur un maillage plus grossier |
| `3 composantes attendues [x, y, z] (maillage 3D …)` | cas 3D (`type = "box"` ou `[mesh.extrude]`) avec un vecteur à 2 composantes | écrire `U = [ux, uy, uz]`, `body_force = [fx, fy, fz]`… |
| `Maillage 3D (…) : … disponible(s) en 2D seulement` | option sans version 3D (axisymétrique, swirl, poreux, disques, couplé, animation, compressible) | la retirer, ou calculer en 2D |
| `reprise d'un calcul 2D sur un maillage 3D impossible` | `restart` d'un calcul d'une autre dimension | repartir d'un calcul de même dimension |

## 2. Le calcul diverge

Le calcul s'arrête et l'affiche. Les résultats de ce calcul ne sont pas écrits.

### Incompressible, stationnaire

```
Erreur : Le calcul a divergé à l'itération 28 (vitesses infinies, non définies ou
démesurées). Pistes : démarrer en convection_U = "upwind", puis repasser en linearUpwind et
poursuivre (--continue, ou « Continuer le calcul précédent ») ; qualité du maillage (ligne
« Maillage » en tête du calcul) ; conditions aux limites ; sous-relaxations plus faibles.
```

Les pistes ont été essayées sur `mesh_naca_multi` : profil à volet, 17 880 triangles,
laminaire Re = 100. Ce cas est l'exemple de maillage `mesh_naca_multi.toml`, complété par
ν = 0.01, une entrée U = (1, 0), une sortie, des symétries haut et bas et deux parois.

| Essai | Résultat mesuré |
|---|---|
| `linearUpwind` (défaut) | diverge à l'itération 28 |
| `linearUpwindLimited` | diverge à l'itération 40 |
| `relax_U = 0.5`, ou pseudo-transitoire (`pseudo_cfl`) | divergent aussi |
| `upwind` seul | ne diverge pas ; Cl = 0.665 après 400 itérations, non convergé |
| `upwind` 400 itérations, puis `linearUpwind` en continuant | ne diverge pas ; Cl = 0.755 ; résidus en baisse régulière (3.5·10⁻⁵ à 800 itérations) |

Ce qu'il faut retenir :
- Le résultat en `upwind` seul n'est pas un résultat final : Cl diffère de 12 % de celui
  obtenu en `linearUpwind`. `upwind` sert à passer le démarrage.
- Sur ce maillage, la ligne de qualité affichait déjà le problème :
  `ATTENTION : asymétrie max 19.25 > 4`.

**Procédure en ligne de commande** :

```
microrans run2d cas.toml -o res --set solver.convection_U=upwind solver.max_iter=400
microrans run2d cas.toml -o res --continue
```

La 2e commande reprend `res/checkpoint.npz` avec les réglages du fichier
(`linearUpwind`) et fait `max_iter` itérations **de plus**.

**Procédure dans l'interface** :
1. Page « 5. Numérique » : « Convection U » = upwind, « Itérations max » = 400.
2. Page « 6. Calcul » : « Lancer le calcul ».
3. Revenir en linearUpwind (page 5).
4. Page 6 : « Continuer le calcul précédent ».

**Qualité du maillage.** Chaque calcul affiche en tête :

```
Maillage : non-orthogonalité max 63.2° (moy. 8.2°), asymétrie max 19.25
  ATTENTION : asymétrie max 19.25 > 4 (seuil usuel : précision et convergence dégradées)
```

La même alerte apparaît dans la page « 2. Maillage » de l'interface et dans
`microrans mesh`. Les seuils (70° et 4) sont ceux d'OpenFOAM. Ce sont des seuils d'alerte
usuels, pas une limite au-delà de laquelle le calcul échoue à coup sûr. Pour améliorer le
maillage :
- mailles plus fines là où la géométrie est fine (`h_surface`) ;
- `growth` plus faible ;
- plus d'itérations du mailleur (`[mesh] max_iter`).

### Incompressible, instationnaire

Le message propose : pas de temps plus petit, ou `adjust_dt = true`, qui ajuste Δt au
nombre de Courant `max_co`. Avec un schéma en temps explicite, un Δt trop grand est
signalé **avant** le calcul :

```
ATTENTION : Δt = … au-delà de la limite de stabilité de rk3 (Co = …, Dn = …, limites … / …) :
risque de divergence. Réduire Δt ou activer adjust_dt.
```

### Compressible

```
Erreur : État non physique (ρ ≤ 0 ou p ≤ 0) dans 12 cellule(s). Pistes : réduire [solver]
cfl (Runge-Kutta : les exemples utilisent 0.8 à 2) ; en stationnaire, steady_scheme =
"implicit" ; démarrer à l'ordre 1 (first_order_iter) ou flux hllc ; vérifier les conditions
aux limites et la qualité du maillage.
```

Essais sur l'exemple `compressible_rampe_mach2` (5 400 cellules ; l'exemple est réglé en
Runge-Kutta RK3, `cfl = 2`) :

| Réglage | Résultat mesuré |
|---|---|
| `cfl = 2` (exemple) | converge, 722 itérations, 12 s |
| `cfl = 4` ou `8` | s'arrête dès les premières itérations (« État non physique ») |
| `cfl = 4` ou `8` avec `first_order_iter = 300` | s'arrête aussi : ne rattrape pas un CFL trop grand |
| `cfl = 8`, `steady_scheme = "implicit"` | converge, 63 itérations, 1.9 s ; même solution (Cp à 5·10⁻⁹ près) |
| `limiter = "none"`, flux roe ou hllc | converge (488 itérations) |

En stationnaire, le schéma implicite est donc à essayer en premier. Il est déjà utilisé
par les exemples NACA transsonique et plaque plane.

En instationnaire, `cfl` est rapporté à la somme des flux sur les faces d'une cellule. Sur
le tube de Sod (une seule rangée de cellules carrées), `cfl = 3` reste stable et précis :
erreur sur ρ identique à celle de `cfl = 0.8`. La limite dépend donc du maillage.

## 3. « NON CONVERGÉ »

```
NON CONVERGÉ après 50 itérations (0.6 s) : résultats à vérifier (augmenter max_iter, voir les résidus).
```

Le calcul stationnaire a atteint `max_iter` avant que les résidus passent sous `tol`. Les
résultats sont écrits, mais ce n'est pas une solution convergée. En ligne de commande, le
code de sortie vaut 1.

À faire :
1. Regarder `convergence.png` (ou `history.csv`, ou le bouton « Historique » de
   l'interface).
2. Si les résidus baissent encore : continuer avec `--continue` ou « Continuer le calcul
   précédent », qui ajoutent `max_iter` itérations.
3. Si les résidus baissent lentement : essayer l'algorithme couplé
   (`[solver] algorithm = "coupled"`). Mesuré sur les exemples (README, § 5.5) : 2 à 27 fois
   moins de temps, sauf pour la transition SST-γ et la convection naturelle, où il est un
   peu plus lent.
4. Si les résidus stagnent ou oscillent sans baisser, l'écoulement réel est peut-être
   instationnaire : essayer `mode = "transient"`. L'inverse n'est pas vrai : un calcul
   stationnaire peut converger alors que l'écoulement réel est instationnaire (§ 4).

## 4. « Convergé », mais résultat douteux

### y⁺ trop grand pour le traitement de paroi

Avec un modèle de turbulence et `wall_treatment = "resolved"` (défaut), la première maille
doit être à y⁺ ≲ 1. Mesuré sur la plaque plane SST (exemple `plaque_plane_loi_de_paroi`,
première maille à y⁺ ≈ 50) :

| Traitement de paroi | Cd (frottement) |
|---|---|
| `"wall_function"` (exemple) | 0.00547 |
| `"resolved"` | **0.00144**, affiché « Convergé » |
| référence : maillage fin, y⁺ ≈ 0.9 | 0.00552 |

Le résumé de fin de calcul le signale désormais :

```
ATTENTION : y⁺ max = 50 sur « plate » avec le traitement résolu (il faut y⁺ ≲ 1 à la 1re
maille) : frottement et traînée sous-estimés …
```

À faire : affiner la première maille (`first_height`, `[mesh.layers]`, ou `grading` des
blocs), ou passer en lois de paroi (SA, k-ω, SST seulement).

Où lire y⁺ :
- colonne `y⁺ max` du tableau des efforts ;
- colonne `yplus` de `wall_<paroi>.csv`.

Dans les exemples fournis, y⁺ max vaut au plus 1.8 en traitement résolu, et 143 avec lois
de paroi.

### Écoulement réellement instationnaire calculé en stationnaire

Mesuré sur le cylindre à Re = 100 : exemple `cylindre_re100_urans`, relancé avec
`mode = "steady"`.

| Calcul | Cd | Cl |
|---|---|---|
| stationnaire | 1.104, « Convergé » en 127 itérations | ≈ 0 |
| références instationnaires (lâcher de tourbillons) | 1.33 à 1.35 en moyenne | oscille, amplitude ≈ 0.33 |

Le calcul stationnaire a trouvé la solution symétrique. Elle existe mathématiquement, mais
l'écoulement réel ne la suit pas : Cd est 17 % trop bas, sans aucun signal.

Derrière un obstacle non profilé (cylindre, plaque en travers, marche), le lâcher de
tourbillons apparaît en 2D vers Re ≈ 47 pour un cylindre (valeur de la littérature).
Au-delà, calculer en `mode = "transient"`.

### Bonnes pratiques générales (non mesurées ici)

Ce sont des règles usuelles en CFD, valables pour tout logiciel :
- **Indépendance au maillage** : refaire le calcul sur un maillage environ 2 fois plus fin
  dans chaque direction. Si Cd ou Cl change encore beaucoup, le premier maillage était trop
  grossier.
- **Taille du domaine** : les frontières d'entrée, de sortie et de champ lointain doivent
  être loin de l'obstacle, sinon elles le « sentent ».
- **Nombre de Reynolds** : vérifier Re = U L / ν. Le résumé l'affiche :
  `ν = 0.01 m²/s, U = 1, L = 1 → Re = U L / ν = 100`. Une erreur d'unité sur ν change tout.

## 5. Calcul trop long ou trop gros

- **Ordre de grandeur** : ~15 µs par cellule et par itération, ce qui donne ~15 s par
  itération pour 1 million de cellules. Ce chiffre est variable selon la machine.
  L'avertissement apparaît au-delà de 500 000 cellules.
- **3D** : ~8.5 à 14 s par itération et ~3 Go pour 10⁶ hexaèdres en laminaire (SST : ~40 %
  de plus), ~1 à 2 min de préparation (maillage, distance à la paroi). Avertissement au-delà
  de 500 000 cellules (README § 7).
- **Algorithme couplé** : voir § 3. La mémoire mesurée va de 0.3 Go (16 000 cellules) à
  1.65 Go (160 000 cellules).
- **Démarrage multigrille** (`fmg_levels`) : le calcul part d'une solution obtenue sur des
  maillages plus grossiers.
- **Plusieurs cœurs** :
  - les points d'un balayage ou d'une polaire se calculent en parallèle (`-j 4`, ou « jobs »
    dans l'interface) ;
  - un calcul seul utilise un cœur (bibliothèque BLAS limitée à un fil : deux calculs
    lancés en même temps ne se gênent plus). Les noyaux Numba (`numba = true`) gagnent
    ~10 % ; le multi-fil était plus lent sur la machine de test (README).
- **Arrêter puis reprendre** : le calcul est sauvegardé (`checkpoint.npz`) toutes les 5
  minutes et à la fin. On le poursuit avec `--continue` ou « Continuer le calcul
  précédent ».

## 6. Signaler un défaut

En ligne de commande, `microrans --debug …` affiche la trace Python complète. Pour signaler
un défaut, joindre :
- le fichier de cas ;
- la commande lancée ;
- le message obtenu.
