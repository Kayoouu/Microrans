# Liste de travail (ordre = priorité ; mise à jour le 2026-10-04 ; audit 2 fait)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

(rien)

## À faire

Lots de correction de l'audit 2 (détail et reproduction : `docs/audit_utilisateur.md`,
partie « Audit 2 » ; campagnes à relancer après chaque lot : `tools/audit/`). Demande de
l'utilisateur : fouiller les bugs et l'expérience **avant** de continuer la 3D.

1. **F1 — résultats faux, calculs perdus** : C11 (polaire / balayage / multigrille hors du
   dossier du cas : dossier de sortie relatif → le rendre absolu), C13 (interface : faces
   d'extrusion renommées → périodicité effacée ; `_load_extrude` compare à back / front
   littéraux), C14 (frontière inexistante dans `[output] forces` → vérifier juste après le
   maillage, avant les itérations ; interface : retirer les frontières disparues), U14
   (interface : maillage périmé → avertir ou remailler), U17 (vecteurs des coupes : échelle
   rapportée à |U| complet, rien si composantes dans le plan négligeables). Un test par
   point, plus relance des campagnes 1, 2 et 4.
2. **F2 — messages** : M13 à M19, U16 (virgule décimale des sondes), L6 (bruit d'arrondi
   du résumé).
3. **F3 — interface 3D et figures** : U15 (sondes / profils suivent la dimension ; case
   Extruder grisée pour axisymétrique / poreux / animation), U18 (« Enregistrer sous » et
   fichiers relatifs), U19 (petits défauts 3D), L5 (plan des figures 3D en ligne de
   commande, noms d'axes ; remplace le petit travail correspondant).
4. **F4 — textes et documentation** : T1, T2, D5 à D7.
5. **Import de maillage 3D (Gmsh .msh).** `Mesh3D` connaît déjà tétraèdres, pyramides,
   prismes, hexaèdres (`VTK_TYPES`) ; vérifier par des tests que le solveur les traite
   (seuls hexaèdres et prismes sont validés aujourd'hui). Frontières depuis les groupes physiques.
   Les coupes des figures (`PlaneSlice`) gèrent déjà tétraèdres et pyramides (testé) ; la
   coupe z de l'interface utilise `ZSlice` (maillages en couches) : passer à PlaneSlice pour z
   si le maillage n'est pas en couches.
   Attention (mesuré au lot E3) : sur des hexaèdres à faces gauches, les volumes de VTK
   diffèrent des nôtres par cellule (1.4 % pour un gauchissement de 0.5 % de la maille,
   proportionnel ; total égal) : découpage différent des faces, pas une erreur ; à
   documenter si l'import amène de tels maillages.
   Fait quand : lecture testée sur un petit fichier versionné, calcul court, doc.
6. **Un calcul sur plusieurs cœurs.** Numba multi-fil mesuré plus lent sur la machine de
   test : mesurer d'abord où part le temps (AMG pyamg mono-fil, assemblage) avant de choisir.
   Ne rien promettre sans mesure.
7. **Maillage en C pour les profils + comparaison NASA TMR** (vérifier d'abord que les
   données TMR sont accessibles depuis l'environnement ; sinon le noter et passer).
8. **Étude de convergence en maillage (GCI)** sur 2 ou 3 cas de validation du README.

## Plus tard (feuille de route du README § 9, non prioritaire)

Transition avec gradient de pression (T3C), corrections de courbure / rotation, loi de paroi
thermique, k-ε haut-Reynolds, compressible turbulent et axisymétrique, couplé avec énergie et
turbulence, viscoélasticité.

## Petits travaux

- Seuil du limiteur compressible par défaut (`venkat_k` 0.05) : remesurer Sod, rampe M = 2,
  plaque avec 0.1 et 0.3 ; changer le défaut seulement si ces cas restent bons (lot #30 :
  0.05 donne deux solutions stationnaires et biaise C_l sur le NACA).

- Audit 2, observations : `mesh_naca_multi` 14 233 triangles → pic 1.2 Go (cause à
  mesurer) ; pas de garde mémoire en ligne de commande pour un très gros maillage 3D ;
  `microrans mesh` écrit toujours dans results/mesh.


## Décisions qui appartiennent à l'utilisateur (ne pas trancher seul)

- Publier une release (tag refusé depuis l'environnement).
- Changer de langage / d'architecture, supprimer une fonction, nouvelle dépendance lourde.
