# Liste de travail (ordre = priorité ; mise à jour le 2026-10-04 ; audit 2 et audit 2 approfondi faits)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

(rien)

## À faire

Lots de correction de l'audit 2 et de l'audit 2 approfondi (détail et reproduction :
`docs/audit_utilisateur.md`, parties « Audit 2 » et « Audit 2 approfondi » ; campagnes à
relancer après chaque lot : `tools/audit/`). Demande de l'utilisateur : fouiller les bugs et
l'expérience **avant** de continuer la 3D. Un test par point corrigé.

1. **F1 — résultats faux, calculs perdus** :
   - P1 d'abord (petit, gros effet) : 1 fil BLAS (OPENBLAS / OMP / MKL_NUM_THREADS) fixé aux
     points d'entrée (cli, gui, exe) avant l'import de numpy, sauf si l'utilisateur l'a
     fixé ; mesurer côte à côte un calcul seul et deux simultanés (plaque : 9.7 → 137 s).
   - C15 vitesse de référence : une seule convention (proposition : vitesse imposée
     maximale partout, ν depuis Re compris ; interface sans `reference_velocity = 1`
     forcé ; U_ref affiché) ; vérifier les exemples avec `reynolds` et U ≠ 1 ; doc.
   - C18 contour CSV « x;y » à virgule décimale : détecter le séparateur `;` et la virgule
     décimale, refuser clairement un contour absurde (points hors d'échelle).
   - C21 balayage : refuser une clé inconnue (liste des clés du schéma de cas).
   - C11 (sortie relative des polaires / balayages / multigrille → absolue), C13
     (`_load_extrude` compare à back / front littéraux), C14 (frontière inconnue de
     `[output] forces` vérifiée après le maillage), U14 (maillage périmé dans l'interface),
     U17 (échelle des vecteurs des coupes).
2. **F1b — plantages et reprises** : C16 (perturbation en 3D), C17 (reprise compressible
   implicite / NS : reprendre l'état du pilotage CFL, sinon ne plus écrire « exact »), C19
   (Ctrl-C : écrire checkpoint + résumé + champs comme le bouton Arrêter), C20 (balayage
   de mesh.nx : entier), C22 (polaire 3D : rotation autour de z), C23 (multigrille sur
   extrusion fine : refuser proprement ou interpoler par couche).
3. **F2 — messages** : M13 à M22 (BOM accepté : `utf-8-sig` ; Latin-1, sortie = fichier,
   une seule maille, dossier masquant un exemple), U16 (virgule décimale des sondes), L6.
4. **F3 — interface 3D, figures, sorties** : U15, U18, U19, U20 (Arrêter pendant le
   maillage), L5, L7 (aire et normale dans les CSV pariétaux), P2 (interface 2 fois plus
   lente : mesurer la cause d'abord).
5. **F4 — textes et documentation** : T1, T2, D5 à D8.
6. **Import de maillage 3D (Gmsh .msh).** `Mesh3D` connaît déjà tétraèdres, pyramides,
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
7. **Un calcul sur plusieurs cœurs.** Numba multi-fil mesuré plus lent sur la machine de
   test (peut-être à cause de P1 : refaire la mesure après F1) : mesurer d'abord où part le temps (AMG pyamg mono-fil, assemblage) avant de choisir.
   Ne rien promettre sans mesure.
8. **Maillage en C pour les profils + comparaison NASA TMR** (vérifier d'abord que les
   données TMR sont accessibles depuis l'environnement ; sinon le noter et passer).
9. **Étude de convergence en maillage (GCI)** sur 2 ou 3 cas de validation du README.

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
