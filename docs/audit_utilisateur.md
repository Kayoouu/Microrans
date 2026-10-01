# Audit utilisateur (débutant et expert) — constats et suivi

Campagnes automatiques, toutes reproductibles, menées comme le ferait un utilisateur :

| Campagne | Contenu | Résultat |
|---|---|---|
| 1. Exemples, ligne de commande | les 22 exemples, `microrans run2d <exemple>`, en entier, figures comprises | 20 terminés et convergés ; 2 échecs (exemples de maillage seul) |
| 2. Exemples, interface | ouverture, maillage, calcul court, tracé de chaque champ, distributions pariétales, profil, historique, enregistrement, réouverture | 18 / 22 sans erreur ; 4 plantages (sondes) |
| 3. Entrées invalides, ligne de commande | 54 cas : fautes de frappe, valeurs impossibles, types faux, fichiers absents, TOML mal formé, `--set` | 13 traces Python brutes ; ~15 erreurs acceptées sans avertissement |
| 4. Combinaisons, interface | 115 : 6 types de maillage × 5 types de corps, 6 modèles × 2 traitements pariétaux, 3 algorithmes, 9 schémas en temps, 6 conditions limites, thermique, axisymétrique, 6 lois non newtoniennes, options, TOML, reprise, arrêt, export, balayage, 1D | 75 sans erreur ; 40 échecs (7 causes) |
| 5. Relecture | README, textes de l'interface (captures de chaque page), aide de la ligne de commande, saisie clavier en français | voir U et D ci-dessous |

Après les lots 1 et 2 (mêmes campagnes, relancées) : exemples en ligne de commande 20 / 20
calculs terminés, résultats identiques au bit près, les 2 exemples de maillage seul arrêtés en
3 s avec le bon message ; exemples dans l'interface 22 / 22 ; entrées invalides : 0 trace
brute (13 avant), 4 acceptées (17 avant), toutes légitimes, 44 erreurs claires, 5
avertissements ; combinaisons de l'interface : 26 / 115 en erreur (41 avant), toutes des
fautes volontaires désormais expliquées.

Lot 3 (interface) : U1 à U10 corrigés, plus trois défauts trouvés en chemin (U11 :
l'interface écrivait ν = 0.01 dans un cas non newtonien qui n'en donnait pas ; U12 :
valeurs arrondies à 6 chiffres réécrites dans le cas ; U13 : section [sweep] ignorée) et
un défaut du mailleur hybride (corps superposés ou trop proches : maillage faux sans
erreur, voir « Corps superposés »). Vérifié sur les 22 exemples : chargés puis
enregistrés par l'interface, ils ne reçoivent plus que des valeurs égales aux défauts du
solveur. Durées des exemples remesurées un calcul à la fois (les campagnes 1 et 2
lançaient 4 calculs en parallèle sur 4 cœurs : durées 2 à 6 fois trop longues).
Campagnes de l'interface relancées après le lot 3 : exemples 20 / 22 sans erreur (les 2
exemples de maillage seul, maillés puis lancés, partaient avant avec un ν inventé par
l'interface ; ils sont refusés avec l'explication) ; combinaisons : mêmes erreurs
volontaires qu'au lot 2, sauf les corps ajoutés en maillage hybride (le cercle ajouté
passe ; un corps ajouté puis changé de type passe aussi, au lieu d'un maillage faux ou
d'un refus « non manifold ») ; aucun plantage ; suite de tests : 375 réussis, 1 ignoré.

L'interface est pilotée hors écran en simulant les actions de l'utilisateur (choix dans
les listes, frappe dans les champs) ; les boîtes de dialogue sont interceptées et
enregistrées.

Statut : **à faire** / **corrigé** (avec le test qui le vérifie).

## C. Résultats faux ou fonctions inutilisables

| # | Constat | Reproduction | Statut |
|---|---|---|---|
| C1 | **Saisie décimale avec Windows en français : le point est supprimé en silence.** « 0.5 » devient 5, « -0.25 » devient −25, « 2.5E-4 » devient 0.0025. L'interface affiche pourtant ses valeurs avec un point. | frappe clavier simulée, langue française | **corrigé** : point et virgule acceptés quelle que soit la langue (`test_numbers_typed_with_dot_or_comma_in_french_locale`) |
| C2 | **Sondes : 4 exemples plantent au lancement dans l'interface** (cylindre instationnaire, tube de Sod, filtre poreux, sang). Le champ « Sondes » transforme la liste de points en texte que le solveur ne relit pas. | ouvrir l'exemple, Lancer | **corrigé** : saisie « x y ; x y », relue comme liste de points (`test_probe_points_roundtrip`) |
| C3 | **Maillages « Triangles » et « Hybride » inutilisables depuis l'interface** (`KeyError : 'type'`) dès que le cas n'a pas de section `[domain]` (nouveau cas et 20 exemples sur 22). | ouvrir un exemple rectangle, passer en Triangles, Générer | **corrigé** : domaine sans type = rectangle (`test_switch_to_triangles_from_rectangle_case_meshes`) |
| C4 | **Clés et sections mal orthographiées ignorées sans avertissement** (ligne de commande et onglet TOML). `max_iters = 5` → calcul complet ; `[solveur]` ignoré. | fichier modifié | **corrigé** : avertissement avec suggestion (« vouliez-vous dire « max_iter » ? », noms français reconnus), clé hors section, clé sans effet pour le solveur ou le type de maillage, `[bodies]` au lieu de `[[bodies]]` ; ligne de commande, interface (onglet TOML, avant maillage, calcul et balayage) ; aucun avertissement sur les 22 exemples ni sur les cas des tests (`tests/test_validate.py`, `test_toml_typo_reported`) |
| C5 | **Valeurs impossibles acceptées sans avertissement** : ν < 0, x1 < x0, rayon < 0, tolérance < 0, Mach < 0, vitesse à 3 composantes, condition pour une frontière inexistante, `nu` et `reynolds` donnés ensemble (l'un est ignoré). | fichier modifié | **corrigé** : erreurs claires avant le maillage, toutes signalées en une fois (ν ≤ 0, x1 ≤ x0, nx non entier ou < 1, tolérance < 0, relaxation hors de ]0, 1], Mach < 0, dt ≤ 0, vitesse ≠ 2 composantes, `nu` et `reynolds` ensemble, algorithme ou schéma de convection inconnus — `linearupwind` en minuscules dégradait en silence au 1er ordre) ; condition pour une frontière absente : avertissement avec suggestion (`tests/test_validate.py`) |
| C6 | **Balayage du nombre de Reynolds sans effet sur un cas donné en ν** (trouvé en écrivant la vérification des clés) : `nu` étant prioritaire, tous les points étaient calculés avec le même ν (cavité, Re = 100 et 1000 : Cd identiques à 6 chiffres), sans message. Touchait l'interface (page Balayage) et `microrans sweep`. | balayage `physics.reynolds` de `cavite_re100` | **corrigé** : balayer `reynolds` retire `nu` du point (et inversement) (`test_reynolds_sweep_on_case_given_in_nu`) |

## M. Messages d'erreur

| # | Constat | Statut |
|---|---|---|
| M1 | 13 erreurs affichent une **trace Python brute** : `nx` non entier ou texte, `max_iter` texte, `relax_U = 0`, `mode = "unsteady"`, instationnaire sans `dt`, `dt = 0`, `dt < 0` (après 9 s de calcul), fichier de maillage absent, corps sans `type`, `--set` à 3 niveaux, `--set` sans `=`, `--set` non numérique. | **corrigé** : plus aucune trace brute sur les 54 entrées de la campagne ; fichier de maillage absent → « [mesh] path : fichier de maillage introuvable » ; toute erreur imprévue → « erreur interne inattendue… à signaler », trace complète avec `microrans --debug` (`tests/test_cli.py`) |
| M2 | **Messages trompeurs** : `nx = 0` → « Trop de segments de progression » ; `ν = 0` → « Factor is exactly singular » ; `[mesh]` absente → « Conditions aux limites manquantes pour bottom, inlet… » ; fichier vide → « donner nu ou reynolds » ; incidence « dix » → message Python en anglais ; `U = "vite"` → « nom inconnu v ». | **corrigé** : `nx = 0` → « nx = 0 : doit être ≥ 1 » ; `ν = 0` → « nu = 0.0 : doit être > 0 » ; `[mesh]` absente, fichier vide, incidence « dix », `U = "vite"` : message direct (`tests/test_validate.py`) |
| M3 | L'interface affiche l'exception brute (« KeyError : 'type' ») sans explication ni renvoi au journal. | **corrigé** : fenêtre « Le calcul n'a pas pu aboutir », message en clair, renvoi au journal ; erreur interne signalée comme telle ; titre selon la tâche (« Le maillage n'a pas pu aboutir »…) (`test_error_dialog_text_is_readable`, `test_mesh_failure_dialog`). Trouvé en corrigeant : un corps sans `type` saisi dans l'onglet TOML faisait planter la liste des corps ; corrigé |
| M4 | Exemples de maillage seul (`mesh_cylindre_hybride`, `mesh_naca_multi`) lancés comme un calcul : « donner nu ou reynolds » **après 89 à 116 s de maillage**, sans dire qu'il s'agit d'un exemple de maillage. | **corrigé** : message immédiat (« ce fichier décrit seulement un maillage… microrans mesh ») en ligne de commande ; dans l'interface, « Lancer » sans maillage ni conditions aux limites explique quoi faire (`test_mesh_only_example_run_explained`) |
| M5 | Maillage de 4 millions de cellules : aucun retour pendant 3 minutes, aucune estimation de mémoire ou de durée. | **en partie** : maillage structuré (rectangle, O, multi-blocs) de plus de 500 000 cellules → avertissement immédiat avec mémoire et durée par itération estimées (mesurées : ~1 Ko et ~15 µs par cellule) ; interface : confirmation avant de lancer. Pas d'estimation pour les triangles (raffinements : une formule simple se tromperait d'un facteur 10) (`test_huge_structured_mesh_announced`) |
| M6 | Cartes graphiques proposées dans l'exécutable alors qu'elles ne peuvent pas y fonctionner ; le message conseille `pip install`. | **corrigé** : choix grisés dans l'exécutable avec explication ; backend inconnu refusé avec la liste (`test_gpu_backends_explained_in_executable`) |
| M7 | Maillage « fichier » sans chemin (« Format non supporté : (msh, su2) », puis `KeyError 'path'`), maillage « multi-blocs » vide (`KeyError 'vertices'`), corps « contour importé » sans fichier (`IsADirectoryError`). | **corrigé** : maillage « fichier » sans chemin, « multi-blocs » vide, contour importé sans fichier : message clair avant le maillage |
| M8 | `"strouhal": NaN` dans `summary.json` pour un calcul instationnaire court (JSON invalide pour d'autres outils). | **corrigé** : NaN et infinis écrits `null` dans summary.json et balayage.json (`test_summary_is_strict_json_when_strouhal_undefined`) |
| M9 | **Divergence signalée par « Factor is exactly singular »** (trouvé en relançant les campagnes : exemple `mesh_naca_multi` lancé avec les conditions devinées par l'interface). Les vitesses atteignaient ~1e50 (finies, donc non détectées) avant que la matrice de pression ne devienne singulière. | **corrigé** : arrêt dès que la vitesse dépasse 10⁶ × max(U_ref, 1) ou n'est plus finie, juste après la quantité de mouvement ; message « Le calcul a divergé à l'itération N » avec des pistes, dans l'ordre d'efficacité constaté sur ce cas (upwind : stable ; relax_U = 0.5 et pseudo_cfl = 5 : divergent encore) (`test_divergence_reported_with_tips`) |

Messages déjà bons (à garder comme modèle) : modèle, condition, schéma, matériel ou type
de maillage inconnus (liste des choix) ; expression dangereuse refusée ; condition limite
manquante ; axisymétrique hors du demi-plan y ≥ 0 ; lois de paroi avec k-ε ou transition ;
reprise introuvable ; TOML mal formé (ligne et colonne).

## U. Interface

| # | Constat | Statut |
|---|---|---|
| U1 | **Tableau des conditions limites illisible** : colonnes écrasées (« cyli… », type « Pa », « Ch »). | corrigé : colonnes sans objet pour le cas masquées (T et q sans thermique, débit sans entrée en vitesse, scalaires, u_θ) ; cases sans effet pour le type de la ligne grisées et non modifiables ; colonne Type à la largeur du libellé ; bulles d'aide sur les en-têtes ; virgule décimale acceptée (« 1,5 ») ; « Ω=abc » ne bloque plus la mise à jour, q, omega et débit non numériques refusés avant calcul |
| U2 | **Colonne centrale trop étroite** : champs et notes coupés, défilement horizontal sur toutes les pages. | corrigé : largeur minimale des pages ramenée sous 480 px (avant : 852, 785, 505, 492 px), libellé au-dessus du champ si la place manque, listes déroulantes compactes ; panneau de réglages élargi (640 px à l'ouverture) |
| U3 | Mode « Nombre de Reynolds » : le champ ν grisé affiche une valeur périmée (0.01 au lieu de 0.05 pour Re = 20). | corrigé : le champ grisé affiche la valeur déduite (ν = U L / Re, ou Re = U L / ν dans l'autre mode) et suit les modifications ; changer de mode garde la même viscosité |
| U4 | Non newtonien : 9 paramètres affichés « défaut » alors qu'aucun défaut n'existe → erreur « paramètres manquants » au lancement. | corrigé : seuls les paramètres de la loi choisie sont affichés, marqués « obligatoire » (ν min / ν max : « auto », défaut expliqué en bulle d'aide), formule de la loi sous la liste ; paramètres manquants ou négatifs et loi mal orthographiée signalés avant le lancement (aussi en ligne de commande) |
| U5 | Maillage en O : « Ajouter » crée un 2e corps interdit ; tout corps ajouté est posé exactement sur le premier. | corrigé : « Ajouter » grisé en maillage en O (bulle d'aide : choisir non structuré ou hybride) ; nouveau cercle placé à droite des corps existants, de taille comparable (maillage hybride réel vérifié : 3 corps, frontières propres) ; changement de type d'un corps (cercle → rectangle, ellipse, NACA) : même place et même taille (trouvé par la campagne de l'interface : la nouvelle forme reprenait les valeurs par défaut, sur le premier corps ; les maillages hybrides « réussis » de ces scénarios au lot 2 étaient en fait faux) ; vérification avant maillage, en interface et en ligne de commande, voir « Corps superposés » ci-dessous |
| U6 | Accueil : 22 exemples non classés (noms de fichiers, descriptions coupées), exemples de maillage seul mélangés aux calculs, rien pour dire par où commencer. | corrigé : exemples classés (« Commencer ici », laminaire, turbulence, thermique, fluides particuliers, compressible, maillage seul), titre lisible et durée mesurée (un calcul à la fois) ; description complète de l'exemple choisi ; bouton « Ouvrir l'exemple » (le double-clic seul n'était pas découvrable) ; même classement dans `microrans examples` ; un test signale tout nouvel exemple non classé |
| U7 | Résumé des résultats en JSON brut (16 chiffres). | corrigé : résumé en phrases (convergence, Re, tableau Cd / Cl / Cm / y⁺ par paroi, part pression / frottement, Strouhal, Nusselt, viscosité, zones poreuses, disques, scalaires, sondes), 4 chiffres, bruit d'arrondi affiché 0 ; toute clé non prévue reste listée ; le temps simulé est maintenant écrit dans summary.json en instationnaire |
| U8 | Textes périmés : accueil et « À propos » (« RANS/URANS », 4 modèles ; ni compressible, ni couplé, ni transition) ; « Stationnaire (RANS) » même en laminaire. | corrigé : accueil, titre de fenêtre et « À propos » décrivent l'ensemble (modèles, algorithmes, thermique, scalaires, non newtonien, poreux, disques, axisymétrique, compressible Euler / Navier-Stokes laminaire) ; « Stationnaire » / « Instationnaire » (RANS / URANS expliqués en bulle d'aide) |
| U9 | Lois de paroi proposées avec k-ε et transition : refus seulement au lancement. | corrigé : option grisée avec k-ε et transition (raison en bulle d'aide) ; en laminaire, réglages de turbulence grisés ; combinaison refusée avant le lancement (aussi en ligne de commande) ; lois de paroi en laminaire signalées (sans objet : écart mesuré 4.5·10⁻⁵ sur la cavité) |
| U10 | Balayage : valeurs mal saisies (« 0 à 4 ») → message Python en anglais (« could not convert string to float »). | corrigé : message en français avec la syntaxe attendue (plage, pas nul, pas de mauvais signe, liste vide) ; virgules décimales permises (« 0:1:0,5 », « 0,5; 1; 1,5 ») |
| U11 | *(trouvé pendant le lot 3)* Fluide non newtonien sans ν dans le cas (ν de référence tiré de la loi, ex. sang ≈ 3.6·10⁻⁵ m²/s) : l'interface affichait ν = 0.01 et l'**écrivait dans le cas** à la première modification, ~280 × la valeur de la loi (Re affiché, diffusivité thermique α = ν/Pr). | corrigé : champ ν vide permis, « auto : ν de la loi à γ̇ = U/L » ; vérifié sur les 22 exemples que l'interface n'écrit plus que des valeurs égales aux défauts du solveur. Conséquence vue par la campagne de l'interface : un exemple de maillage seul, maillé puis lancé, partait avec ce ν inventé ; il est maintenant refusé avec l'explication (physique et conditions limites à renseigner, ou ouvrir un exemple de calcul) |
| U12 | *(trouvé pendant le lot 3)* Champs numériques arrondis à 6 chiffres (ν = 1/550 → 0.00181818, t_end = 0.000632455532 → 0.000632456) puis réécrits dans le cas. Écart ~10⁻⁶, sans effet visible, mais valeur modifiée sans action de l'utilisateur. | corrigé : affichage court quand il est exact (3000, 0.1), sinon tous les chiffres |
| U13 | *(trouvé pendant le lot 3)* Section [sweep] du cas ignorée par l'interface : l'exemple de polaire ouvert puis « Lancer » ne calculait qu'un point ; le balayage gardait ses valeurs par défaut (-4:12:2 au lieu de -4:14:2). | corrigé : [sweep] chargé dans la page Calcul (paramètre, valeurs, continuation, calculs en parallèle) avec une note « Lancer le balayage » ; réglages du balayage réécrits dans le cas au lancement |

## L. Ligne de commande

| # | Constat | Statut |
|---|---|---|
| L1 | Description « RANS/URANS » incomplète ; pas de `--version` ; aide mêlant anglais et français. | **corrigé** : `microrans --version` ; description complète et exemple de premier calcul ; aide entièrement en français (« utilisation », « arguments », « afficher cette aide »), toutes les options expliquées (12 n'avaient aucune aide : `--no-plot`, `-q`, `-v`, `--max-iter`…) ; erreurs de saisie en français et courtes (« commande « runn » inconnue — vouliez-vous dire « run2d » ? », « --alpha : 3 valeurs attendues », « « a » n'est pas un nombre ») au lieu de l'usage complet suivi d'un message anglais (`test_help_in_french`, `test_usage_errors_in_french`) |
| L2 | Fin de calcul : bloc JSON de 40 lignes au lieu d'un résumé lisible. | **corrigé** : même résumé en phrases que l'interface (U7), incompressible et compressible ; « mode steady » → « stationnaire » ; summary.json toujours écrit en entier (`test_end_of_run_summary_readable`) |
| L3 | `--set` : deux niveaux seulement (`boundary.lid.U=…` plante). | **corrigé** : `--set` à plusieurs niveaux (`boundary.lid.U=[2,0]`), indice de liste (`bodies.0.radius=0.3`), virgule décimale ; `--set` sans `=` ou clé incomplète → message clair (`tests/test_cli.py`) |

## D. Documentation

| # | Constat | Statut |
|---|---|---|
| D1 | Pas de tutoriel pas à pas pour un premier calcul. | à faire |
| D2 | Pas de section dépannage (divergence, « non convergé », y⁺, qualité du maillage). | à faire |
| D3 | Pas de glossaire (RANS, SIMPLE, y⁺, Cp, Cf, patch, O-grid…). | à faire |
| D4 | Pas de référence complète des clés du fichier de cas ; l'extrait du README oublie `sst_gamma`. | à faire |

## Corps superposés ou trop proches (mesuré pendant le lot 3, corrigé)

Mesures sur deux cercles de rayon 0.5, couches de paroi de 0.054 d'épaisseur
(n = 4, première maille 0.01, raison 1.2) :

| Disposition | Hybride (avant) | Non structuré (avant) |
|---|---|---|
| séparés (écart ≥ 2 × épaisseur des couches) | correct | correct |
| identiques | refus « non manifold » après 6 s (20 s sur le cas de l'audit) | réussi, 2e corps **disparu en silence** |
| l'un dans l'autre | « réussi » mais **faux** : couches de cellules dans le solide, frontières parasites `_b1_top`, `_b2_top` | 2e corps disparu (géométriquement juste) |
| écart 0.08 (< 2 × 0.054) | « réussi » mais faux (mêmes frontières parasites) | correct |
| corps à cheval sur le bord du domaine | « réussi » mais faux (frontière parasite) | corps coupé par le bord (usage légitime : bosse sur une paroi) |

Maintenant, avant de mailler (vérification du cas, < 0,1 s) : en hybride, erreur pour les
corps qui se recouvrent, se touchent, sont plus proches que deux épaisseurs de couches, ou
sont à moins d'une épaisseur du bord du domaine (avec l'écart mesuré et l'épaisseur) ; en
non structuré, avertissement pour les corps qui se recouvrent (obstacle unique) et pour un
corps hors du domaine ; en maillage en O, erreur si le nombre de corps n'est pas 1. Corps
décrits par un fichier de contour : non vérifiés.

## Observations à approfondir (non traitées)

- **Performance** : sur la cavité (laminaire, SIMPLEC), le temps par itération passe de
  0.036 s (10 000 cellules) à 0.79 s (40 000 cellules), soit ×5 par cellule, au seuil de
  20 000 cellules où le solveur de pression passe de la factorisation directe (LU) à l'AMG.
  Le seuil ou le réglage de l'AMG est peut-être mal placé ; à mesurer avant de conclure.
- **Robustesse sur maillage déformé** : sur le maillage de `mesh_naca_multi` (asymétrie max
  19, non-orthogonalité 63°, fente entre les deux éléments), le cas laminaire Re = 100 diverge
  vers l'itération 28 avec `convection_U = "linearUpwind"` (défaut), reste stable en `upwind`.
  Une limitation de la correction linearUpwind dans les cellules très asymétriques (comme le
  `cellLimited` d'OpenFOAM) rendrait le défaut plus robuste ; non fait.

## Ce qui fonctionne (vérifié)

Tous les modèles de turbulence (laminaire, SA, k-ε, k-ω, SST, SST-γ), les trois
algorithmes stationnaires, les neuf schémas en temps, le pas adaptatif, les six types de
conditions limites sur la cavité, la thermique, les scalaires, les zones poreuses, les
disques actuateurs, le démarrage multigrille, le pseudo-transitoire, les sondes (saisies
dans l'interface), l'arrêt sur efforts, la reprise, l'arrêt en cours de calcul, l'export du
maillage, le balayage, le canal 1D (RANS et URANS), les 20 exemples de calcul en ligne de
commande. Aucune valeur NaN ou aberrante dans les 102 résumés produits (hors M8).
