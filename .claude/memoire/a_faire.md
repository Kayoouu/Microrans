# Liste de travail (ordre = priorité ; mise à jour le 2026-10-19 ; jalons A, B, C)

Méthode pour chaque tâche : mesurer d'abord (profil, chiffres de départ), changer, mesurer
côte à côte dans les mêmes conditions, tests, docs (README, chiffres), commit, CI verte,
mémoire à jour, compte rendu. Une tâche à la fois. Si une tâche s'avère impossible ou
inutile (mesure à l'appui), l'écrire dans le journal et passer à la suivante.

## En cours

Rien (audit 3 partie 2, `monitor_tol` mesuré et documenté, le 2026-10-19). Suivant :
point 1 ci-dessous (reste de l'audit 3, partie 2).

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

1. **Audit 3, partie 2** : chaque accord de validation du README face au maillage et à la
   taille du domaine (cylindre et NACA déjà trouvés ; après G1, les chiffres ayant bougé) ;
   figures restantes (pression, vorticité, ν_t, maillages) ; interface (figures de la page
   Résultats). Convergence apparente trompeuse (C_l de l'O-grid, ordre apparent 5.8) à
   relire. Pseudo-temps (`pseudo_cfl`) sur maillage non orthogonal : lent sans
   limiteur (NACA 2 063 → > 3 000 it), lent aussi avec (résidus U_y 1e-4) : à comprendre ;
   étude O / C (naca_o_c.py, ψ = 0.5 fixé) à refaire sans pseudo_cfl (le README § 8.6 cite
   encore ses chiffres, dont « ordre apparent 5.8 » en O). `monitor_tol` : fait le 19
   (mesuré, pas de défaut d'extremum, documenté). Compressible : pas de garde-fou sur
   l'arrêt par les résidus, et son monitor_tol (autre critère) n'est pas mesuré. Puis audit 4 : le jalon A est atteint quand un audit complet trouve 0 résultat
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

- Docstring de microrans/fv2d/case.py (sorties [output]) : ajouter « tracées dans
  probes.png » aux sondes (G2 ; au prochain commit qui touche microrans/).
- Profils sur un plan ou un axe de symétrie (mesuré au lot G2) : U_z de la cavité cubique
  4.3e-4 |U| sur ses lignes médianes (z = 0.5, sur des faces de mailles), U_y de l'éolienne
  3e-4 sur l'axe, du filtre poreux 1.5e-3 sur l'axe ; nuls en théorie. Interpolation
  (cellule d'un seul côté + gradient) ou solution non symétrique ? Mesurer aux centres de
  part et d'autre.
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
