# Liste de travail (ordre = priorité ; mise à jour le 2026-10-05 ; jalons A, B, C)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

(rien ; F4 fait le 2026-10-08)

## À faire

Organisé en jalons, ordre accepté par l'utilisateur le 2026-10-05 (« oui écris ») : A, puis
B ; C non visé. Durées = estimations, pas des mesures. Demande de l'utilisateur : fouiller
les bugs et l'expérience **avant** de continuer la 3D. Un test par point corrigé ; détail
des points d'audit : `docs/audit_utilisateur.md` ; campagnes à relancer après chaque lot :
`tools/audit/`.

### Jalon A — 2D académique crédible (≈ 9 à 16 lots ; visé fin octobre à mi-novembre 2026)

Usages : TP, projets d'étudiants, études paramétriques 2D laminaire / RANS sur cas
classiques. **Atteint quand** : un audit complet ne trouve plus aucun résultat faux
silencieux, et les points 1 et 2 sont faits (ou leur impossibilité écrite). Côté
utilisateur (je ne peux pas le faire) : publier une release, faire essayer l'exécutable
Windows par un humain sur un vrai PC (jamais fait : la CI vérifie seulement qu'il démarre).

1. **Étude de convergence en maillage (GCI)** sur 2 ou 3 cas de validation du README.
2. **Maillage en C pour les profils + comparaison NASA TMR** (plaque plane, NACA 0012, SA /
   SST). Vérifier d'abord que les données TMR sont accessibles depuis l'environnement ;
   sinon l'écrire (validation partielle) et passer.
3. **Audit 3** (même méthode que l'audit 2 approfondi : campagnes `tools/audit/` + nouvelles),
   puis lots de correction. S'il trouve encore un résultat faux silencieux : corriger puis
   refaire un audit ; le jalon n'est pas atteint avant.

### Jalon B — 3D académique, petite géométrie (≈ 1.5 à 3 mois de plus ; janvier–février 2027)

Périmètre : incompressible RANS, jusqu'à ~10⁶ cellules, maillage fait dans Gmsh. Limite
mesurée : 10⁶ cellules = 8.5 à 14 s / itération, ~3 Go → un stationnaire prend des heures.

7. **Import de maillage 3D (Gmsh .msh).** `Mesh3D` connaît déjà tétraèdres, pyramides,
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
8. **Cas 3D turbulent comparé à des données publiées** (vérifier d'abord que les données
   sont accessibles depuis l'environnement).
9. **Un calcul sur plusieurs cœurs.** Numba multi-fil mesuré plus lent avant P1 : refaire
   la mesure ; mesurer où part le temps (AMG pyamg mono-fil, assemblage) avant de choisir.
   Si pas de gain mesuré : l'écrire et s'arrêter. Ne rien promettre sans mesure.
10. **Compressible turbulent** : seulement si l'utilisateur le demande (5 à 10 lots).

### Jalon C — pré-industriel : non visé

Pas atteignable avec l'architecture actuelle : 10⁷ cellules ≈ 30 Go et 1.5 à 2 min /
itération (extrapolé, non mesuré) ; pas de mailleur 3D à couches prismatiques ni d'import
de CAO ; la confiance demande une validation par d'autres. MPI / cœur compilé = changement
d'architecture : décision de l'utilisateur, ne pas commencer seul. Position recommandée :
pré-dimensionnement et enseignement ; OpenFOAM ou SU2 pour les décisions.

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
- Après F3 : la page Résultats de l'interface ne lit pas `[output] slice_axis` /
  `slice_value` (plan des figures en ligne de commande) ; un plan hors du domaine n'est
  signalé qu'après le calcul (figures omises, ATTENTION), pas avant.


## Décisions qui appartiennent à l'utilisateur (ne pas trancher seul)

- Publier une release (tag refusé depuis l'environnement).
- Changer de langage / d'architecture (dont MPI ou cœur compilé, jalon C), supprimer une
  fonction, nouvelle dépendance lourde.
- Compressible turbulent (jalon B, point 10) : seulement sur demande.
- Jalon A, côté utilisateur : release, essai de l'exécutable Windows sur un vrai PC.
