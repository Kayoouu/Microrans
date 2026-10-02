# Liste de travail (ordre = priorité ; mise à jour le 2026-10-02 ; E2 fait)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

(rien)

## À faire

1. **E3 — VTK binaire.** fields.vtk ASCII = 169 Mo pour 10⁶ cellules. Écrire le format
   « legacy » binaire (gros-boutiste) en option ou par défaut ; vérifier la relecture (lecteur
   de test dans le dépôt ; meshio n'est pas une dépendance et ne doit pas le devenir pour
   ça). Fait quand : taille et temps d'écriture mesurés, test de relecture exacte, doc.
2. **E4 — Coupes x = cte et y = cte dans l'interface (3D).** Au moins pour les pavés
   (cellules alignées) ; coupe générale d'hexaèdres / prismes par un plan si le coût reste
   raisonnable. Fait quand : choix de l'axe et de la cote dans la page Résultats, tests hors
   écran, étape ajoutée à `_selftest`, captures regardées.
3. **#30 — Écart transsonique NACA 0012 (C_l ~4.5 % bas).** Pistes déjà écartées dans
   `docs/compressible.md` § 3.5. Limiter à une journée de travail ; conclure même si l'écart
   reste inexpliqué (dire ce qui a été vérifié).
4. **Import de maillage 3D (Gmsh .msh).** `Mesh3D` connaît déjà tétraèdres, pyramides,
   prismes, hexaèdres (`VTK_TYPES`) ; vérifier par des tests que le solveur les traite
   (seuls hexaèdres et prismes sont validés aujourd'hui). Frontières depuis les groupes physiques.
   Fait quand : lecture testée sur un petit fichier versionné, calcul court, doc.
5. **Un calcul sur plusieurs cœurs.** Numba multi-fil mesuré plus lent sur la machine de
   test : mesurer d'abord où part le temps (AMG pyamg mono-fil, assemblage) avant de choisir.
   Ne rien promettre sans mesure.
6. **Maillage en C pour les profils + comparaison NASA TMR** (vérifier d'abord que les
   données TMR sont accessibles depuis l'environnement ; sinon le noter et passer).
7. **Étude de convergence en maillage (GCI)** sur 2 ou 3 cas de validation du README.

## Plus tard (feuille de route du README § 9, non prioritaire)

Transition avec gradient de pression (T3C), corrections de courbure / rotation, loi de paroi
thermique, k-ε haut-Reynolds, compressible turbulent et axisymétrique, couplé avec énergie et
turbulence, viscoélasticité.

## Petits travaux

- Erreur ruff E731 préexistante `tests/test_mesh2d.py:74` (lambda assignée).

## Décisions qui appartiennent à l'utilisateur (ne pas trancher seul)

- Publier une release (tag refusé depuis l'environnement).
- Changer de langage / d'architecture, supprimer une fonction, nouvelle dépendance lourde.
