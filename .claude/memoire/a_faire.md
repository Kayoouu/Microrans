# Liste de travail (ordre = priorité ; mise à jour le 2026-10-05 ; audit 2 et audit 2 approfondi faits)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

(rien ; F1b fait le 2026-10-05)

## À faire

Lots de correction de l'audit 2 et de l'audit 2 approfondi (détail et reproduction :
`docs/audit_utilisateur.md`, parties « Audit 2 » et « Audit 2 approfondi » ; campagnes à
relancer après chaque lot : `tools/audit/`). Demande de l'utilisateur : fouiller les bugs et
l'expérience **avant** de continuer la 3D. Un test par point corrigé.

1. **F2 — messages** : M13 à M22 (BOM accepté : `utf-8-sig` ; Latin-1, sortie = fichier,
   une seule maille, dossier masquant un exemple), U16 (virgule décimale des sondes), L6.
2. **F3 — interface 3D, figures, sorties** : U15, U18, U19, U20 (Arrêter pendant le
   maillage), L5, L7 (aire et normale dans les CSV pariétaux), P2 (interface 2 fois plus
   lente : mesurer la cause d'abord).
3. **F4 — textes et documentation** : T1, T2, D5 à D8.
4. **Import de maillage 3D (Gmsh .msh).** `Mesh3D` connaît déjà tétraèdres, pyramides,
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
5. **Un calcul sur plusieurs cœurs.** Numba multi-fil mesuré plus lent sur la machine de
   test (peut-être à cause de P1 : refaire la mesure après F1) : mesurer d'abord où part le temps (AMG pyamg mono-fil, assemblage) avant de choisir.
   Ne rien promettre sans mesure.
6. **Maillage en C pour les profils + comparaison NASA TMR** (vérifier d'abord que les
   données TMR sont accessibles depuis l'environnement ; sinon le noter et passer).
7. **Étude de convergence en maillage (GCI)** sur 2 ou 3 cas de validation du README.

## Plus tard (feuille de route du README § 9, non prioritaire)

Transition avec gradient de pression (T3C), corrections de courbure / rotation, loi de paroi
thermique, k-ε haut-Reynolds, compressible turbulent et axisymétrique, couplé avec énergie et
turbulence, viscoélasticité.

## Petits travaux

- Seuil du limiteur compressible par défaut (`venkat_k` 0.05) : remesurer Sod, rampe M = 2,
  plaque avec 0.1 et 0.3 ; changer le défaut seulement si ces cas restent bons (lot #30 :
  0.05 donne deux solutions stationnaires et biaise C_l sur le NACA).

- Balayage d'une clé `bodies.0.radius` : refusé clairement depuis F1 (« clé impossible à
  modifier ») ; à rendre possible si un besoin apparaît.
- Ctrl-C sous Windows (C19) : non vérifié (test ignoré sous win32) ; à essayer avec
  l'exécutable si l'occasion se présente.

- Audit 2, observations : `mesh_naca_multi` 14 233 triangles → pic 1.2 Go (cause à
  mesurer) ; pas de garde mémoire en ligne de commande pour un très gros maillage 3D ;
  `microrans mesh` écrit toujours dans results/mesh.


## Décisions qui appartiennent à l'utilisateur (ne pas trancher seul)

- Publier une release (tag refusé depuis l'environnement).
- Changer de langage / d'architecture, supprimer une fonction, nouvelle dépendance lourde.
