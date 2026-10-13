# Liste de travail (ordre = priorité ; mise à jour le 2026-10-13 ; jalons A, B, C)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

Rien (audit 3, partie 1, terminé le 2026-10-13).

## À faire

Organisé en jalons, ordre accepté par l'utilisateur le 2026-10-05 (« oui écris ») : A, puis
B ; C non visé. Durées = estimations, pas des mesures. Demande de l'utilisateur : fouiller
les bugs et l'expérience **avant** de continuer la 3D. Un test par point corrigé ; détail
des points d'audit : `docs/audit_utilisateur.md` ; campagnes à relancer après chaque lot :
`tools/audit/`.

### Jalon A — 2D académique crédible (≈ 9 à 16 lots ; visé fin octobre à mi-novembre 2026)

Usages : TP, projets d'étudiants, études paramétriques 2D laminaire / RANS sur cas
classiques. **Atteint quand** : un audit complet ne trouve plus aucun résultat faux
silencieux, l'étude GCI (faite le 2026-10-09) et le maillage en C (fait du 10 au 12
octobre ; comparaison NASA TMR impossible depuis l'environnement, écrite au README) sont
faits. Côté
utilisateur (je ne peux pas le faire) : publier une release, faire essayer l'exécutable
Windows par un humain sur un vrai PC (jamais fait : la CI vérifie seulement qu'il démarre).

1. **Lot G1 — C24 : solution stationnaire indépendante de `relax_U` / `pseudo_cfl`**
   (audit 3 ; détail et reproduction : docs/audit_utilisateur.md, Audit 3, C24). Mesuré :
   NACA SA, `pseudo_cfl = 20` déplace la solution convergée de +0.75 % sur C_d (point fixe
   différent, pas un arrêt prématuré) ; README : `relax_U` 0.7 → 0.5, C_d −0.6 %. Étapes :
   (a) mesurer l'effet de `relax_U` et `pseudo_cfl` sur 3 à 4 exemples (cylindre Re 20,
   cavité, NACA SA, plaque SA) ; (b) corriger l'interpolation de Rhie-Chow (correction de
   relaxation type Majumdar 1988 ; en pseudo-temps, terme comme `ddtCorr`) ; test : même
   solution convergée à 1e-6 près pour deux `relax_U` et avec / sans `pseudo_cfl` ;
   (c) cas à solution exacte et GCI inchangés en ordre ; (d) remesurer les chiffres
   publiés qui bougent (README, tutoriel, en-têtes d'exemples, tests à bornes) ; corriger
   README § 8 point 6 (explication fausse par `monitor_tol`). Estimé 2 lots. Fait quand :
   test vert, chiffres publiés remesurés, CI verte.
2. **Lot G2 — figures** (audit 3 : L8 sondes tracées comme résidus, L9 bruit d'arrondi à
   pleine échelle dans les profils, L10 cadrage du tuyau axisymétrique, L11 figures de
   champ du tube de Sod). Estimé 1 lot. Un test par point.
3. **Audit 3, partie 2** : chaque accord de validation du README face au maillage et à la
   taille du domaine (cylindre et NACA déjà trouvés ; après G1, les chiffres ayant bougé) ;
   figures restantes (pression, vorticité, ν_t, maillages) ; interface (figures de la page
   Résultats). Convergence apparente trompeuse (C_l de l'O-grid, ordre apparent 6.7) à
   relire. Puis audit 4 : le jalon A est atteint quand un audit complet trouve 0 résultat
   faux silencieux.

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
- Après GCI : l'exemple `cylindre_re20` (champ lointain R = 40) donne C_d 1.5 % au-dessus
  de l'estimation R → ∞ (≈ 2.007, extrapolée de R = 20, 40, 80) ; un domaine plus grand
  changerait les chiffres documentés (README, tutoriel, tests 2.0 < C_d < 2.1) : à mesurer
  (coût, R = 160 ?) avant de décider. Cylindre Re = 100 : même question (limite 7).
- Cavité Re = 100 : u min converge vers −0.21395 ± 1e-4, à 1.4 % de Ghia (−0.21090) ;
  pas de référence plus fine vérifiée depuis l'environnement (recherche du 2026-10-09).
- Exemples NACA (R = 30) : erreur de domaine C_l −1.5 %, C_d +3.5 % (R → ∞ extrapolé,
  2026-10-11). Options : R = 100 ou 300 (+9 % / +19 % de cellules à même 1re maille et
  même progression) ou correction tourbillonnaire du champ lointain ; change les chiffres
  documentés et la polaire : à mesurer (coût) avant de décider.
- Campagnes : c8 écrit des `_audit_*.toml` temporaires dans microrans/examples ; ne pas
  lancer c1 / c6 en même temps que c8 (c6 les a ramassés le 2026-10-11). Scripts à
  corriger (faux constats du 2026-10-13) : c2 ajoute `max_steps` (clé inconnue) et change
  `case_path` sans copier le contour de mesh_naca_multi ; c4 `box_depuis_cylindre` passe
  par `set_combo` (signaux bloqués) ; c14 limite à 900 s (RK3 du cylindre : 16 à 30 min).
- Après convergence, itérations 4 à 10 fois plus chères (résidus à l'arrondi ; plaque
  compressible, cavité cubique ; mesuré sur machine partagée) : remesurer seul, puis voir
  les tolérances des solveurs linéaires.
- NASA TMR (profils, plaque plane) : à faire si l'utilisateur autorise
  turbmodels.larc.nasa.gov (refusé du 10 au 12 octobre 2026).
- Exemples NACA : maillage en C (probablement plus juste : sillage proche résolu ; C_l
  −1.7 %, C_d −3 % par rapport à l'O) et domaine plus grand (R = 30 : C_l −1.5 %, C_d
  +3.5 %) ; les deux effets se compensent en partie sur C_l. Coût mesuré : C 12 288
  cellules, 3 000 it ≈ 180 s et « NON CONVERGÉ » à tol 1e-5, contre O 834 it / 34 s.
  Change les chiffres documentés : décider après mesure complète (C, R = 100).


## Décisions qui appartiennent à l'utilisateur (ne pas trancher seul)

- Publier une release (tag refusé depuis l'environnement).
- Changer de langage / d'architecture (dont MPI ou cœur compilé, jalon C), supprimer une
  fonction, nouvelle dépendance lourde.
- Compressible turbulent (jalon B, point 10) : seulement sur demande.
- Jalon A, côté utilisateur : release, essai de l'exécutable Windows sur un vrai PC.
