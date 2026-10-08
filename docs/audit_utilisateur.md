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

Lot 4 (ligne de commande et documentation) : L1, L2 et D4 d'abord (aide et erreurs en
français, résumé lisible en fin de calcul, référence complète des clés), puis tutoriel,
guide de dépannage et glossaire (D1 à D3). Ces guides ont été écrits en lançant chaque
commande et en essayant chaque piste. Ce travail a fait trouver sept défauts : une traînée
4 fois trop faible sans signal quand y⁺ ne convient pas au traitement de paroi (C8),
l'interface qui remplaçait en silence les valeurs absentes de ses listes (C9), un résultat
stationnaire « convergé » mais faux pour un écoulement instationnaire (C10, documenté),
la qualité du maillage annoncée mais pas affichée (M10), des pistes compressibles absentes
ou inefficaces (M11), un refus tardif (M12) et l'absence de copie d'un exemple en ligne de
commande (L4). Suite de tests : 393 réussis, 1 ignoré.

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
| C7 | **`linearUpwindLimited` sans effet pour la vitesse, la température et u_θ** (trouvé en écrivant la référence des clés, lot 4) : accepté pour `convection_U` et `convection_T`, mais le gradient limité n'était transmis que pour la turbulence et les scalaires ; vitesse identique au bit près à `linearUpwind` (cavité 24 × 24 : écart 0). L'option qu'on essaie justement sur un maillage déformé ne faisait rien. | écart mesuré entre les deux schémas | **corrigé** : gradient limité transmis pour U, T et u_θ (écart désormais 1.2·10⁻², entre linearUpwind et upwind) ; schéma par défaut inchangé (`test_limited_scheme_acts_on_velocity_and_temperature`). Sur `mesh_naca_multi` il retarde la divergence (itération 40 au lieu de 28) sans l'empêcher : voir « Observations » |
| C8 | *(trouvé en écrivant le dépannage, lot 4)* **y⁺ trop grand en traitement résolu : traînée fausse sans aucun signal.** Plaque plane SST sur le maillage de l'exemple à lois de paroi (1re maille à y⁺ ≈ 50), `wall_treatment = "resolved"` : Cd = 0.00144 au lieu de 0.0055 (maillage fin : 0.00552 ; lois de paroi : 0.00547), affiché « Convergé ». | `plaque_plane_loi_de_paroi`, `--set solver.wall_treatment="resolved"` | **corrigé** : avec un modèle de turbulence en traitement résolu, y⁺ max > 5 sur une paroi → « ATTENTION : y⁺ max = 50 sur « plate »… frottement et traînée sous-estimés » dans le résumé (ligne de commande et interface, mention dans la ligne d'état) et dans summary.json (`warnings`) ; aucun exemple fourni ne la déclenche (y⁺ max ≤ 1.8 en résolu) (`test_yplus_too_high_for_resolved_wall_flagged`) |
| C9 | *(trouvé en vérifiant le tutoriel)* **L'interface remplaçait en silence une valeur absente de ses listes** : un cas en `convection_U = "linearUpwindLimited"` ouvert puis lancé depuis l'interface partait en `linearUpwind` (option absente de la liste) ; même chose pour toute valeur hors liste. | ouvrir le cas, Lancer | **corrigé** : `linearUpwindLimited` ajouté à la liste ; toute valeur hors liste est gardée (option « … (valeur du fichier) »), une faute est donc signalée par la vérification au lieu d'être remplacée (`test_out_of_list_values_kept`) |
| C10 | *(mesuré pour le dépannage)* **Écoulement instationnaire calculé en stationnaire : « Convergé », résultat faux.** Cylindre Re = 100 en `mode = "steady"` : convergé en 127 itérations, Cd = 1.104 (solution symétrique instable) contre 1.33–1.35 en moyenne pour l'écoulement réel, à lâcher de tourbillons. | `cylindre_re100_urans --set solver.mode="steady"` | **documenté** (`docs/depannage.md` § 4) ; pas de détection automatique (un critère général serait peu fiable) |

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
| M10 | *(lot 4)* Le message de divergence renvoyait à la qualité du maillage « affichée après le maillage », mais `microrans run2d` ne l'affichait pas ; l'interface montrait les valeurs sans alerte au-delà des seuils. | **corrigé** : ligne « Maillage : non-orthogonalité max …, asymétrie max … » en tête de chaque calcul (incompressible et compressible), alerte au-delà de 70° / 4 (seuils d'OpenFOAM), même alerte dans la page Maillage ; sur `mesh_naca_multi` : « ATTENTION : asymétrie max 19.25 > 4 », juste avant la divergence ; le message de divergence dit aussi comment repasser en linearUpwind (`--continue`, « Continuer le calcul précédent ») (`test_mesh_quality_shown_at_run_start`, `test_mesh_quality_alert_shown`) |
| M11 | *(lot 4)* Compressible : « Divergence à l'itération N. » sans aucune piste ; « État non physique » proposait de démarrer à l'ordre 1, sans effet mesuré sur un CFL trop grand, et pas le schéma implicite. | **corrigé** : pistes dans l'ordre mesuré sur la rampe Mach 2 (RK3 : cfl 2 converge, 4 et 8 s'arrêtent, même avec `first_order_iter = 300` ; implicite à cfl 8 : même solution à 5·10⁻⁹ près, 1.9 s au lieu de 12 s) (`test_unphysical_state_gives_measured_tips`) |
| M12 | *(lot 4)* Fluide non newtonien avec un modèle de turbulence : refus seulement au lancement (pas à la vérification de l'interface). | **corrigé** : refusé à la vérification, en interface et en ligne de commande (`tests/test_validate.py`) |

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
| L1 | Description « RANS/URANS » incomplète ; pas de `--version` ; aide mêlant anglais et français. | **corrigé** : `microrans --version` ; description complète et exemple de premier calcul ; aide entièrement en français (« utilisation », « arguments », « afficher cette aide »), toutes les options expliquées (21, sur 8 commandes, n'avaient aucune aide : `--no-plot`, `-q`, `-v`, `--max-iter`, `--range`…) ; erreurs de saisie en français et courtes (« commande « runn » inconnue — vouliez-vous dire « run2d » ? », « --alpha : 3 valeurs attendues », « « a » n'est pas un nombre ») au lieu de l'usage complet suivi d'un message anglais (`test_help_in_french`, `test_usage_errors_in_french`) |
| L2 | Fin de calcul : bloc JSON de 40 lignes au lieu d'un résumé lisible. | **corrigé** : même résumé en phrases que l'interface (U7), incompressible et compressible ; « mode steady » → « stationnaire » ; summary.json toujours écrit en entier (`test_end_of_run_summary_readable`) |
| L3 | `--set` : deux niveaux seulement (`boundary.lid.U=…` plante). | **corrigé** : `--set` à plusieurs niveaux (`boundary.lid.U=[2,0]`), indice de liste (`bodies.0.radius=0.3`), virgule décimale ; `--set` sans `=` ou clé incomplète → message clair (`tests/test_cli.py`) |
| L4 | *(trouvé en écrivant le tutoriel)* Aucun moyen simple d'obtenir une copie modifiable d'un exemple : le fichier est dans le dossier d'installation (dossier temporaire pour l'exécutable). | **corrigé** : `microrans examples NOM [-o DOSSIER]` copie l'exemple et les fichiers qu'il lit (contour du profil), sans rien écraser, modifiable même si l'installation est en lecture seule (`test_example_copied_for_editing`) |

## D. Documentation

| # | Constat | Statut |
|---|---|---|
| D1 | Pas de tutoriel pas à pas pour un premier calcul. | **corrigé** : [`docs/tutoriel.md`](tutoriel.md), interface et ligne de commande : cavité (comparaison à Ghia, écart max 0.004), changement de Reynolds, cylindre Re = 20 (Cd 2.037 contre 2.05), contrôle du maillage (Cd 2.033 avec 4 fois plus de cellules) ; toutes les commandes et sorties obtenues en suivant le texte, parcours de l'interface rejoué hors écran |
| D2 | Pas de section dépannage (divergence, « non convergé », y⁺, qualité du maillage). | **corrigé** : [`docs/depannage.md`](depannage.md) : messages réels, et pour chaque piste le cas où elle a été essayée et le résultat, y compris celles qui n'ont pas marché ; a fait trouver C8, C10, M10, M11, M12 |
| D3 | Pas de glossaire (RANS, SIMPLE, y⁺, Cp, Cf, patch, O-grid…). | **corrigé** : [`docs/glossaire.md`](glossaire.md), chaque terme relié à la clé ou au fichier où on le rencontre ; guides accessibles depuis le menu Aide de l'interface et le README (`test_help_menu_links_the_guides`) |
| D4 | Pas de référence complète des clés du fichier de cas ; l'extrait du README oublie `sst_gamma`. | **corrigé** : [`docs/reference_cas.md`](reference_cas.md), générée à partir des clés que le logiciel vérifie (toutes les sections, signification, solveur et types de maillage concernés ; défauts de `[solver]` lus dans le code) ; un test échoue si elle n'est pas régénérée après un changement de clé ; défauts écrits à la main vérifiés un par un contre le code (tous exacts ; la liste des schémas des scalaires omettait `linearUpwindLimited`, qui est le défaut) ; README : `sst_gamma` ajouté et lien vers la référence |

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
  à l'itération 28 avec `convection_U = "linearUpwind"` (défaut), reste stable en `upwind`
  (non convergé à 10⁻⁵ en 600 itérations). Le gradient limité (`linearUpwindLimited`, actif
  pour la vitesse depuis C7) retarde seulement la divergence (itération 40) : la cause n'est
  donc pas (seulement) le dépassement de la reconstruction ; à chercher (correction non
  orthogonale à 63°, couches de paroi très aplaties). Non fait. La piste « démarrer en
  upwind puis repasser en linearUpwind » est vérifiée sur ce cas, en ligne de commande
  (`--continue`) : plus de divergence, résidus en baisse régulière (3.5·10⁻⁵ à 800
  itérations) ; Cl passe de 0.665 (upwind seul) à 0.755, le résultat upwind n'est donc
  pas final.
- **Interface plus lente que la ligne de commande** : cavité, 351 itérations, 8.6 s dans
  l'interface (hors écran) contre 3.5 s en ligne de commande. Suivi en direct des résidus
  probablement en cause ; à mesurer.
- **Type des frontières avant le calcul** : la page Maillage affiche « lid (patch) »,
  « walls (patch) » pour la cavité, alors que ce sont des parois : le type n'est fixé que
  par les conditions limites, au lancement. Peut surprendre un débutant.
- **Rampe Mach 2** : l'exemple est réglé en Runge-Kutta (12 s) ; en implicite à CFL 8 il
  donne la même solution en 1.9 s. Exemple non modifié.
- **Maillage d'un million de cellules** : après l'avertissement de taille, rien ne s'affiche
  pendant plus de 2 minutes (maillage et préparation du solveur, sans indication de
  progression).

## Ce qui fonctionne (vérifié)

Tous les modèles de turbulence (laminaire, SA, k-ε, k-ω, SST, SST-γ), les trois
algorithmes stationnaires, les neuf schémas en temps, le pas adaptatif, les six types de
conditions limites sur la cavité, la thermique, les scalaires, les zones poreuses, les
disques actuateurs, le démarrage multigrille, le pseudo-transitoire, les sondes (saisies
dans l'interface), l'arrêt sur efforts, la reprise, l'arrêt en cours de calcul, l'export du
maillage, le balayage, le canal 1D (RANS et URANS), les 20 exemples de calcul en ligne de
commande. Aucune valeur NaN ou aberrante dans les 102 résumés produits (hors M8).

---

# Audit 2 (2026-10-04) — après la 3D (lots D4, E à E4) et le lot #30

Même démarche que l'audit 1, centrée sur ce qui a changé depuis : 3D (pavé, extrusion,
interface, coupes x / y / z, VTK binaire), `trailing_edge`, `venkat_k`. Les campagnes sont
maintenant **versionnées** dans `tools/audit/` pour être relancées après les corrections :

```bash
python tools/audit/c1_exemples_cli.py SORTIE [exemple …]          # campagne 1
python tools/audit/c2_exemples_interface.py SORTIE [exemple …]    # campagne 2
python tools/audit/c3_entrees_invalides.py SORTIE [variante …]    # campagne 3
python tools/audit/c3b_cles_numeriques.py SORTIE [contexte]       # campagnes 3b / 3c
python tools/audit/c4_combinaisons_3d_interface.py SORTIE [scénario …]  # campagne 4
```

| Campagne | Contenu | Résultat |
|---|---|---|
| 1. Exemples, ligne de commande | la commande écrite en tête de chacun des 25 exemples, copiée telle quelle, lancée dans un dossier vide, un exemple à la fois | 24 / 25 ; **la polaire s'arrête au 2e point** (C11) ; aucune valeur non finie dans les résumés ; durées conformes au catalogue (canal 3D : 23.6 s au premier lancement, 10.6 s ensuite : cache des polices de Matplotlib construit une fois) |
| 2. Exemples, interface | les 25 : ouverture, enregistrement immédiat, maillage, calcul court, tous les champs, 5 grandeurs pariétales, profil, vue du maillage, coupes x / y / z (3D), « Continuer » | aucun plantage ; **canal_turbulent_3d calculé faux** (C13) ; exemples enregistrés : seules des valeurs égales aux défauts ajoutées, sauf ce cas |
| 3. Entrées invalides, ligne de commande | 44 variantes (pavé, extrusion, vecteurs à 2 composantes en 3D, options 2D seulement, commandes polar / sweep / mesh sur des cas 3D, sorties, `trailing_edge`, `venkat_k`) ; 3b / 3c : 116 valeurs (texte, puis −1) sur 46 clés numériques, puis 12 clés reprises dans un cas où elles servent | 36 / 44 corrects (29 refus clairs en moins d'une seconde, 7 cas valides acceptés ou avertis à juste titre) ; 0 trace Python brute ; 8 à reprendre (M13 à M16, M18, M19, garde mémoire) ; clés « expert » non vérifiées (M18) |
| 4. Combinaisons 3D, interface | 38 scénarios + 5 complémentaires : extrusion de 13 exemples 2D par la case à cocher, pavé × 6 modèles + loi de paroi, 4 schémas instationnaires, thermique 3D, 4 types d'extrémités en z, passage 2D ↔ 3D après maillage, changement de type de maillage, sondes, balayage, Continuer, export, plans de coupe au bord et hors domaine | 26 sans aucun message ; 6 refus attendus et clairs (axisymétrique, poreux, animation, maillage seul, export .msh / .su2) ; 1 scénario faussé par le script (refait à part) ; 5 défauts : **maillage périmé utilisé sans avertissement** (U14), sondes / profils non convertis en 3D (U15), saisie des sondes (U16) ; à part : C14 côté interface |
| 5. Relecture | aide de la ligne de commande, « À propos », README (§ 1, 3, 8, 9), tutoriel, glossaire, figures produites, 4 captures de l'interface en 3D | **vecteurs de bruit numérique affichés comme un écoulement** (U17) ; textes périmés (T1 à T3) ; documentation 3D (D5 à D7) |

Les boîtes de dialogue sont interceptées (titre et texte enregistrés). Piège corrigé en
cours de campagne : `set_combo` bloque les signaux de Qt ; les changements de liste doivent
passer par `setCurrentIndex` pour imiter un clic (deux faux constats évités ainsi).

Non testé ici : les exécutables eux-mêmes (téléchargement des artefacts refusé par le
proxy de l'environnement ; l'auto-test de l'interface tourne sur chaque exécutable en CI),
Windows, écran réel ou haute densité, cartes graphiques, calculs de plusieurs heures.

Statut : **à faire** (lot prévu), puis **corrigé** avec le test qui le vérifie.

## C. Résultats faux, calculs perdus

| # | Constat | Reproduction | Statut |
|---|---|---|---|
| C11 | **Polaire, balayage (série et parallèle) et démarrage multigrille cassés hors du dossier du cas.** Le dossier de sortie relatif est pris dans le dossier courant, la reprise relative à côté du fichier de cas. `microrans polar naca0012_polaire --alpha -4 14 2` (commande du README et de l'en-tête) : arrêt au 2e point après 32 s, « Fichier de reprise introuvable : …/microrans/examples/results/… » (reprise jamais demandée). Ne marche que lancé depuis le dossier du cas. L'interface n'est pas touchée (dossier de sortie absolu). | `microrans polar naca0012_polaire --alpha 0 2 2` ; `microrans sweep cavite_re100 … -j 2` ; `microrans run2d cavite_re100 --set solver.fmg_levels=1` | corrigé (F1) : dossiers de sortie rendus absolus (cas, compressible, balayage) ; polaire 3 / 3 points depuis un autre dossier, multigrille idem ; test |
| C13 | **Interface : canal_turbulent_3d calculé faux sans message.** Les faces d'extrusion y sont renommées (`names = { back = "z0", front = "z1" }`) ; l'interface ne reconnaît la périodicité que pour les noms `back` / `front` : `periodic = [["z0", "z1"]]` est effacé à l'enregistrement et au lancement, puis z0 et z1 sont devinées « paroi ». Canal → conduite fermée : C_d du fond 0.00791 au lieu de 0.008889 (−11 %). Tout cas à faces d'extrusion renommées est touché. | ouvrir l'exemple, Lancer | corrigé (F1) : faces z comparées à leurs noms (`names`) ; canal_turbulent_3d relu et réenregistré à l'identique (périodicité z0 / z1 gardée) ; interface = ligne de commande au bit près (campagne 6 : 0 écart hors chemins, avant 36) ; test |
| C14 | **Frontière inexistante dans `[output] forces` : tout le calcul est fait, puis « Erreur interne inattendue (KeyError) » à la fin**, sans résumé. Pour un long calcul, tout est perdu. Interface : même cause quand on change le type de maillage d'un exemple (cylindre_re20 → Rectangle : `forces = ["cylinder"]` reste → « Paramètre manquant dans le cas : 'cylinder' »). | `microrans run2d cavite_re100 --set 'output.forces=["lidd"]'` | corrigé (F1) : frontières de `forces` vérifiées juste après le maillage, avant les itérations (refus en 0.3 s, liste des frontières existantes) ; test |
| U14 | **Interface : maillage périmé utilisé sans avertissement.** Maillage généré, `nx` changé de 64 à 16 dans le formulaire, « Lancer » : calcul sur l'ancien maillage (4 096 cellules), le cas enregistré dit nx = 16. Même cause, case « Extruder » cochée ou décochée après le maillage : messages incompréhensibles (« [boundary.lid] U = [1.0, 0.0] : 3 composantes attendues (maillage 3D) » juste après être repassé en 2D ; « body_force = [np.float64(0.0), …] »). | voir campagne 4, scénarios `decoche_apres_maillage`, `coche_apres_maillage` | corrigé (F1) : réglages [mesh] / [domain] / [[bodies]] mémorisés au maillage ; s'ils ont changé, le calcul (ou le balayage) remaille et le journal le dit ; test (avant : 36 cellules au lieu de 64) |
| U17 | **Coupes 3D : bruit numérique dessiné comme un écoulement.** Conduite carrée laminaire, coupe x = cte, « Vecteurs vitesse » : grandes flèches désordonnées alors que max \|U_y\|, \|U_z\| = 7·10⁻¹⁶ (U_x = 7.3) : l'échelle automatique des flèches agrandit le bruit d'un facteur ~10¹⁶ et fait croire à un écoulement secondaire qui n'existe pas. | capture 3 de la campagne 5 | corrigé (F1) : composantes dans le plan < 0.1 % de |U| : pas de flèches, titre explicite ; écoulement secondaire réel dessiné agrandi, taille en % de |U| dans le titre (cavité cubique, plan x = 0.5 : 11 %) ; test |

## M. Messages

| # | Constat | Statut |
|---|---|---|
| M13 | `periodic = [["inlet", "outlett"]]` (pavé ou extrusion) → « Erreur : 'outlett' is not in list » (message Python brut, ni liste des frontières, ni suggestion). | corrigé (F2) : « periodic = ["walls", "lidd"] : frontière « lidd » inexistante. Frontières du maillage : walls, lid — vouliez-vous dire « lid » ? » (2D et 3D) ; test |
| M14 | Clés de `[solver]` compressible données en texte → « Erreur interne inattendue » (voir M18). | corrigé (F2) : clés vérifiées avant le calcul (voir M18) ; test |
| M15 | Type de frontière mal écrit (`patch_types = { z0 = "symetrie" }`) → « Type de patch inconnu 'symetrie' (z0) », sans liste des choix ni suggestion. | corrigé (F2) : « Type de frontière inconnu « symetrie » (frontière walls). Choix : wall, patch, symmetry, empty — vouliez-vous dire « symmetry » ? » (noms français reconnus) ; test |
| M16 | Clé dans la mauvaise section (`[physics] moment_center`, qui va dans `[output]`) → « clé inconnue, ignorée », sans indiquer la bonne section. | corrigé (F2) : « clé inconnue, ignorée — clé de [output] : la déplacer » (toutes les sections qui la connaissent) ; test |
| M17 | « Conditions aux limites manquantes pour les patches ['front'] » : liste Python, mot anglais. | corrigé (F2) : « … manquantes pour : back, front — ajouter [boundary.back] et [boundary.front] avec type = un de : … ; faces d'extrusion : type = "symmetry" ou [mesh.extrude] periodic » ; test |
| M18 | **Clés « expert » de `[solver]` non vérifiées** (campagnes 3b / 3c, chaque clé dans un cas où elle sert). Texte → « Erreur interne inattendue » pour `venkat_k`, `limiter_freeze`, `entropy_fix`, `cn_theta`, `ddt_phi_coeff`, `cfl_growth`, `viscous_factor` ; message Python en anglais pour `linear_iter` (« invalid literal for int() »), `linear_tol` (« could not convert string to float ») ; −1 accepté en silence pour `venkat_k`, `limiter_freeze`, `entropy_fix`, `ddt_phi_coeff`, `nonorth_limit`, `cfl_growth`, `cfl_cuts`, `linear_*`, `viscous_factor` ; `cn_theta = -1` → « Le calcul a divergé » au lieu d'un refus. Bien vérifiées : `relax_*`, `pseudo_*`, `max_co`, `max_dt`, `n_outer`, `n_corr`, `n_nonorth`, `max_iter`, `tol`, `monitor_*`, `dt`, `t_end`, `fmg_levels`, `cfl`, `cfl_max`, `order`, `[flow]`. | corrigé (F2) : toutes les clés listées vérifiées avant le calcul (nombre, entier, bornes, choix, true / false), erreurs groupées ; cn_theta entre 0.5 et 1 ; dt et t_end vérifiés aussi en stationnaire ; campagnes 3b / 3c refaites : 116 cas (92 + 24), tous refusés clairement avant le calcul ; les 23 exemples passent ; test |
| M19 | Mineurs : `trailing_edge` sur un corps qui n'est pas un NACA, ignoré sans avertissement ; `microrans urans --steps 50` (abréviation acceptée de `--steps-per-period`) → « steps_per_period doit être un multiple de n_phases » (noms internes) ; commande d'en-tête de `mesh_naca_multi` (`--type unstructured`) → « ATTENTION : [mesh] layers : sans effet ». | corrigé (F2) : trailing_edge hors NACA → ATTENTION ; « --steps-per-period = 50 : doit être un multiple de 8 …, par exemple 48 ou 56 » (abréviations d'options laissées à argparse) ; en-tête de mesh_naca_multi : les deux commandes, l'avertissement expliqué ; test |

## U. Interface

| # | Constat | Statut |
|---|---|---|
| U15 | Extrusion : vitesses, force volumique, gravité, U initiale passent à 3 composantes, mais pas les sondes ni les profils `[[output.lines]]` → refus au lancement, après le maillage (`melange_deux_courants` : « [[output.lines]] « sortie » : start et end à 3 composantes attendus »). Case « Extruder » cochable sur un cas axisymétrique, poreux ou animé : refus seulement au lancement (message clair). Extrusion absente pour le maillage multi-blocs (permise dans le fichier de cas). | corrigé (F3) : sondes et lignes passent à 3 composantes (z = milieu du domaine) et reviennent à 2 ; « Extruder » refusé tout de suite (message, case décochée) si le cas a des options 2D seulement ; extrusion proposée pour le maillage multi-blocs |
| U16 | Sondes : virgule décimale refusée (« 0,25 0,75 » lu comme 4 nombres, la virgule séparant x et y), alors qu'elle est acceptée dans les autres champs (C1) ; refus seulement au lancement. | corrigé (F2) : « 0,25 0,75 » accepté (virgule décimale quand les coordonnées sont séparées par des espaces) ; sondes illisibles refusées avant le lancement ; test |
| U18 | « Enregistrer sous » dans un autre dossier : les fichiers relatifs du cas (contour `profil_volet.dat` de `mesh_naca_multi`, fichier de maillage) ne sont plus trouvés (message clair, mais rien pour les suivre ; la ligne de commande les copie avec `microrans examples`). | corrigé (F3) : contours `[[bodies]] type = "file"` et maillage importé copiés à côté du nouveau cas (sans écraser un fichier différent : chemin absolu alors) ; fichier de reprise : chemin absolu |
| U19 | Petits défauts 3D : profil par défaut de (0, 0, 0) à (1, 0, 0) même hors du domaine (conduite : x ∈ [0, 0.5]) ; « Ouvrir l'animation » actif en 3D ; « Zoom sur les corps » coché mais sans effet en coupe x / y (documenté) ; vue 3D du maillage : légende sur le dessin, étiquette z collée à la figure voisine ; types de frontière « patch » / « cyclic » affichés avant les conditions limites (observation de l'audit 1, toujours là). | corrigé (F3) : profil par défaut selon x par le centre du domaine (ligne saisie gardée) ; « Ouvrir l'animation » actif seulement s'il y a une animation ; « Zoom sur les corps » grisé en coupe x / y ; légende de la vue 3D sous la vue, écart avec la coupe ; types de frontière en clair (paroi, symétrie, périodique ; « patch » : sans type) |

## L. Ligne de commande et figures

| # | Constat | Statut |
|---|---|---|
| L5 | Figures 3D de la ligne de commande toujours dans le plan z médian : pour la conduite carrée, bande de 2 mailles le long de l'axe, sans intérêt (la section x = cte est la bonne figure) ; pas de noms d'axes ; grande marge blanche. | corrigé (F3) : `[output] slice_axis` (x, y, z) et `slice_value` (vérifiés avant le calcul ; plan hors du domaine : avertissement, figures omises) ; conduite carrée : `slice_axis = "x"` ; noms d'axes sur toutes les figures de champs ; domaine haut : figure moins large (conduite, plan z : 840 × 1350 px au lieu de 1500 × 1350) |
| L6 | Résumé de fin de calcul 3D : « Vitesse moyenne dans le domaine : (17.64, -4.447e-16, -2.495e-39) » : bruit d'arrondi affiché (U7 l'avait supprimé ailleurs). | corrigé (F2) : composantes < 10⁻⁹ |U| affichées 0 ; test |

## T / D. Textes et documentation

| # | Constat | Statut |
|---|---|---|
| T1 | La 3D manque dans la description générale : `microrans --help` (« Écoulements 2D en volumes finis »), « À propos » (« écoulements 2D en volumes finis »), description du paquet (« RANS/URANS 1D et 2D »). | corrigé (F4) : « 2D et 3D » dans `microrans --help`, « À propos » (paragraphe 3D), description du paquet, description du fichier Windows et docstring du paquet |
| T2 | Interface, ouverture d'un cas 3D : « figures dans un plan z = constante » (périmé depuis E4 : x, y ou z). | corrigé (F4) : « figures dans un plan x, y ou z = constante (page Résultats, « Plan de coupe ») » |
| D5 | README § 3 « Cas 3D » : `periodic = … # ou patch_types = { back = "symmetry", … }` laisse croire que `patch_types` suffit ; il faut aussi `[boundary.back] type = "symmetry"` (sinon « Conditions aux limites manquantes pour les patches ['back', 'front'] »). | corrigé (F4) : extrait complet (cylindre extrudé, calculable) et variante symétrie avec `[boundary.back]` / `[boundary.front]` ; vérifié : `patch_types` seul refusé, les conditions seules suffisent ; les deux variantes calculées par un test |
| D6 | Tutoriel : la 3D tient en un paragraphe ; pas de pas-à-pas (pavé ou extrusion, conditions en z, coupes, ouverture de `fields.vtk` dans ParaView). | corrigé (F4) : tutoriel § 7 « Un cas 3D » (interface et ligne de commande, conduite carrée : +0.38 % sur 32², +1.50 % sur 16² ; extrusion du cylindre sur une couche : C_d 2.037 comme en 2D ; ParaView) ; commandes exécutées par un test |
| D7 | Glossaire : aucun terme 3D (extrusion, hexaèdre / prisme, faces périodiques, plan de coupe). Dépannage : pas d'entrée pour les nouveaux messages 3D de M13 à M17. | corrigé (F4) : glossaire (pavé / extrusion, hexaèdre / prisme, frontières périodiques, plan de coupe) ; dépannage : type de frontière inconnu, periodic sans effet, slice_axis en 2D, figures omises, extrusion impossible ; extraits confrontés aux vrais messages par un test |

## Observations (non classées en défauts)

- `mesh_naca_multi` : 14 233 triangles en 39 s avec un pic de mémoire de 1.2 Go
  (`mesh_cylindre_hybride` : 21 393 cellules, 25 s, 0.57 Go) : cause à mesurer.
- Maillage 3D de 8·10⁶ cellules : l'avertissement annonce ~23.5 Go (machine de test :
  15 Go) puis le maillage commence quand même en ligne de commande (l'interface demande
  confirmation). Pas de garde par rapport à la mémoire disponible.
- `microrans mesh` écrit toujours dans `results/mesh` : deux maillages successifs
  s'écrasent (`run2d` utilise le nom du cas).
- Premier lancement : ~13 s de plus (cache des polices de Matplotlib), une seule fois.

## Ce qui fonctionne (vérifié)

Les 25 exemples en ligne de commande (sauf la polaire hors du dossier du cas) et dans
l'interface ; en 3D : pavé, extrusion calculée des maillages rectangle et O (cylindre,
profil), maillage hybride extrudé (42 786 cellules) ; les 6
modèles de turbulence et la loi de paroi ; SIMPLE, SIMPLEC ; Euler, BDF2, Crank-Nicolson,
RK3 ; thermique (gravité à 3 composantes) ; extrémités périodiques, symétrie, parois, à
régler ; sondes et profils 3D ; balayage ; « Continuer » (reprise exacte) ; export VTK
(.msh et .su2 refusés avec explication) ; coupes x / y / z, plan sur le bord, plan hors du
domaine expliqué ; refus clairs et immédiats des options 2D seulement et des valeurs
impossibles du pavé et de l'extrusion ; profil d'entrée parabolique refusé en 3D (0.6 s,
message clair) ; vecteurs à 2 composantes en 3D signalés ;
sonde hors domaine signalée ; canal 1D (RANS, 4 modèles en 3.9 s).

## Lots de correction

Remplacés par la liste mise à jour à la fin de la partie « Audit 2 approfondi ».

---

# Audit 2 approfondi (2026-10-04) — cohérence des résultats, fichiers, arrêts, exécutable

Demande : « n'hésite pas à approfondir l'audit ». Au lieu de chercher ce qui plante, ces
campagnes cherchent **ce qui donne un résultat faux** par des contrôles de cohérence qui ne
demandent pas de référence extérieure (deux chemins qui doivent donner la même chose), puis
les fichiers et les usages réels (poste Windows, arrêt d'un calcul, commandes de la
documentation, exécutable).

```bash
python tools/audit/c6_interface_egale_cli.py SORTIE      # interface = ligne de commande
python tools/audit/c7_reprise_exacte.py SORTIE           # N + N itérations = 2N
python tools/audit/c8_2d_egale_3d_une_couche.py SORTIE [N=…]   # 2D = 3D une couche
python tools/audit/c9_export_maillages.py SORTIE         # export / relecture, OpenFOAM
python tools/audit/c9b_maillage_importe.py SORTIE        # calcul sur maillage réimporté
python tools/audit/c10_fichiers_windows.py SORTIE        # BOM, CRLF, Latin-1, chemins, CSV
python tools/audit/c11_ctrl_c.py SORTIE                 # Ctrl-C pendant un calcul
python tools/audit/c12_cas_degeneres.py SORTIE           # cas limites
python tools/audit/c13_efforts_csv.py SORTIE            # efforts = intégration des CSV
python tools/audit/c14_commandes_doc.py SORTIE           # commandes du README et des guides
python tools/audit/c15_arrets_interface.py SORTIE        # « Arrêter », mémoire, durée
MICRORANS_EXE=dist/microrans/microrans python tools/audit/c1_exemples_cli.py SORTIE  # exécutable
```

| Campagne | Contrôle | Résultat |
|---|---|---|
| 6. Interface = ligne de commande | 22 exemples (Sod ignoré : pas de pas de temps fixe), 20 itérations (ou 3 pas), résumés comparés | identiques au bit près (hors chemins de fichiers), sauf **canal_turbulent_3d** (C13 : U moyen 18.07 → 14.11) et **conduite_carree_3d** (C15 : C_d × 9) |
| 7. Reprise exacte | 2N itérations d'un coup contre N + `--continue` N, champs de fields.vtk et efforts, 18 exemples | **exacte au bit près** pour tout l'incompressible (2D, 3D, SA, SST, SST-γ, lois de paroi, thermique, scalaires, non newtonien, poreux, disque, axisymétrique, instationnaire) et le compressible RK3 ; **pas exacte en compressible implicite et Navier-Stokes** (C17) ; refaite après F1b : 18 exemples sur 18 exacts au bit près (tube de Sod non découpable : pas de temps par CFL) |
| 8. 2D = 3D une couche entre deux symétries | 9 exemples, à 20 itérations puis convergés | convergé : écart ≤ 2·10⁻⁶ (convection, plaques SA / laminaire / lois de paroi, NACA SA) : l'affirmation du README tient ; **la perturbation de sillage plante en 3D** (C16, corrigé en F1b) |
| 9. Maillages exportés | 6 types (rectangle, périodique, O, blocs, triangles, hybride) en .msh / .su2 / .vtk / OpenFOAM | relecture exacte (aire à 2·10⁻¹⁶) ; polyMesh OpenFOAM correct (propriétaire < voisin, ordre triangulaire supérieur, normales, cellules fermées, frontières contiguës) ; types de frontière et périodicité non transportés par .msh / .su2 (formats) |
| 9b. Calcul sur maillage réimporté | 4 cas, sans `patch_types` | même résultat convergé (2·10⁻¹¹) |
| 10. Fichiers « Windows » | BOM, CRLF, Latin-1, espaces et accents, sortie occupée, contour CSV français, Gmsh 4.1, JSON | CRLF, accents, JSON, Gmsh 4.1 ASCII : bons ; **BOM refusé** (M20) ; **contour CSV français : 14 Go puis processus tué** (C18) |
| 11. Ctrl-C | cavité, cylindre, NACA compressible | **rien n'est écrit** (C19, corrigé en F1b) |
| 12. Cas limites | vitesse nulle, 1 maille, Re 10⁸, plages vides / inversées / pas nul, dt minuscule, balayages | plages et valeurs bien refusées ; **balayage de `mesh.nx` : erreur interne** (C20) ; **balayage d'une clé mal écrite accepté en silence** (C21) |
| 13. Post-traitement | efforts recalculés à partir de wall_*.csv et des faces (cylindre Re 20, NACA SA, conduite 3D) | pression : identique au résumé (C_d de pression du cylindre 1.224690, NACA 0.005859) ; frottement : identique sur la conduite (0.026460), 0.19 % plus bas sur le cylindre (L7) |
| 14. Documentation rejouée | 25 commandes du README, du tutoriel, du dépannage, dossier vide | 17 justes ; 3 du README en échec (C11) ; 3 exemples fictifs ; RK3 du cylindre : 1 061 s sur machine libre, au-dessus des « 2 à 16 min » annoncés (D8) |
| 15. Interface : arrêts, durée, mémoire | « Arrêter » pendant maillage / calcul / balayage ; cavité ; 400 pas instationnaires | arrêt du calcul en 2.4 s, reprise exacte ; balayage 0.8 s ; pas de fuite mémoire ; **arrêt sans effet pendant le maillage** (U20) ; interface 2 fois plus lente (P2) |
| Exécutable | construit localement (même spec que la CI), 25 exemples, côte à côte avec Python | 24 / 25 (polaire : C11) ; même vitesse que Python ; **fils BLAS : deux calculs simultanés jusqu'à 13 fois plus lents** (P1) |
| Relecture du code | hypothèses 2D dans les chemins 3D, exceptions avalées | **polaire 3D : plantage** (C22), **multigrille sur extrusion fine : plantage** (C23) (corrigés en F1b) ; exceptions avalées : toutes volontaires |

## C. Résultats faux, calculs perdus (suite)

| # | Constat | Reproduction | Statut |
|---|---|---|---|
| C15 | **Vitesse de référence par défaut : trois conventions.** Documentation (« défaut 1 ») et ν = U_ref L / Re : 1 ; coefficients et Re affichés par la ligne de commande : vitesse imposée maximale (ou vitesse initiale) ; interface : écrit `reference_velocity = 1` dans le cas. Même fichier `conduite_carree_3d` : C_d 0.02397 en ligne de commande (U_ref = 3), 0.2157 dans l'interface (U_ref = 1), facteur 9. Cas `reynolds = 100` avec une entrée à U = 2 : calcul et résumé à « Re = U L / ν = 200 », sans avertissement. | campagne 6 ; cas rectangle, entrée U = 2, `reynolds = 100` | corrigé : une seule règle (`choose_reference_velocity`) pour ν = U L / Re, les coefficients, γ̇_ref et C_T : `reference_velocity` donnée, sinon vitesse d'entrée (moyenne sur la frontière, la plus grande avec avertissement), sinon paroi mobile, sinon 1 ; jamais la vitesse initiale ; interface sans 1 forcé ; U_ref et son origine affichés. reynolds = 100 + entrée U = 2 : ν = 0.02 (Re = 100) ; conduite : U_ref = 1 partout ; tests |
| C16 | `[initial] perturbation` (perturbation du sillage) plante en 3D : « operands could not be broadcast together with shapes (6144,3) (2,) » ; le champ est proposé par l'interface en 3D. | cylindre_re100_urans extrudé | corrigé (F1b) : en 3D, centre [x, y] (tube selon z) ou [x, y, z] ; test |
| C17 | Reprise compressible implicite et Navier-Stokes **pas exacte** alors que le résumé dit « mode exact » et le README « identique au bit près » : après 30 itérations, C_d 0.0401 d'un coup contre 0.0629 en 15 + 15 (NACA transsonique). Convergé : même résultat (plaque laminaire C_d 0.00423611 / 0.00423610) mais 304 itérations au lieu de 266 (+14 %) : l'état du pilotage (CFL) n'est pas repris. | campagne 7 | corrigé (F1b) : pilotage du CFL implicite (CFL, plafond, stagnation), comptage de la normalisation des résidus, limiteur gelé, circulation du tourbillon et derniers relevés d'efforts enregistrés dans `checkpoint.npz` ; relevés d'efforts calés sur l'itération globale. N + N = 2N au bit près : NACA transsonique 120 + 120 et 90 + 90 (fenêtre 50), arrêt sur efforts à l'itération 240 retrouvé en reprenant à 205, plaque laminaire avec `limiter_freeze`, rampe. Calcul continu identique au bit près à l'ancien code (5 variantes). Réglages du solveur ou écoulement amont changés : pilotage neuf, indiqué dans `summary.json` (`restart.controller`). Test (échoue sur l'ancien code) |
| C18 | **Contour CSV « à la française » (x;y, virgule décimale, export Excel français) mal lu sans message** : « 1,000000;-0,000000 » est découpé sur les virgules et les points-virgules → points absurdes ((0, 998379)) → le mailleur monte à 14 Go et le système tue le processus ; avec une limite de 4 Go : « mémoire insuffisante : réduire le nombre de cellules » (trompeur). Sur un poste de 16 Go : gel ou disparition de l'exécutable. | campagne 10 | corrigé (F1) : « ; » ou tabulation sans point → virgule décimale ; mêmes points que le CSV « x,y » ; test |
| C19 | **Ctrl-C en ligne de commande : tout est perdu** (« Interrompu. », aucun fichier : ni checkpoint, ni résumé, ni champs) ; le README promet un checkpoint « à l'arrêt demandé » (vrai pour le bouton de l'interface seulement) ; aucun autre moyen d'arrêter proprement un calcul en ligne de commande. | campagne 11 | corrigé (F1b) : 1er Ctrl-C : fin de l'itération en cours puis tous les fichiers écrits, code 130 ; 2e Ctrl-C : arrêt immédiat. Campagne 11 refaite : 3 cas sur 3 avec checkpoint, résumé, champs, historique, CSV de paroi, reprise correcte ; balayage parallèle arrêté en 0.4 s avec le tableau des 2 points finis. Test (Linux ; Windows non vérifié) |
| C20 | Balayage du maillage (`--param mesh.nx --values 8 16`, étude de convergence en maillage) : « Erreur interne inattendue (TypeError : 'float' object cannot be interpreted as an integer) ». | campagne 12 | corrigé (F1, avancé) : valeurs entières des clés de maillage converties ; test |
| C21 | **Balayage d'une clé mal orthographiée accepté en silence** : `--param physics.nuu` ou `solver.max_iterr` → points identiques étiquetés de valeurs différentes (une faute dans `physics.reynold` ferait conclure que Re est sans effet). | campagne 12 | corrigé (F1) : clé balayée inconnue ou sans effet pour le cas refusée avant tout calcul, avec la clé proche (« reynolds ? ») ; test |
| C22 | Polaire d'un cas 3D (profil extrudé) : plantage au 2e point (« matmul: Input operand 1 has a mismatch … (size 3 is different from 2) ») : la rotation de U∞ est écrite en 2D. | profil extrudé, `polar --alpha 0 4 2` | corrigé (F1b) : rotation des deux premières composantes seulement ; test |
| C23 | Démarrage multigrille (`fmg_levels`) d'un cas extrudé sur peu de couches (canal_turbulent_3d) : la triangulation 3D de l'interpolation échoue (centres coplanaires au niveau grossier) et affiche ~25 lignes d'aide de Qhull en anglais. Pavé : correct. | `run2d canal_turbulent_3d --set solver.fmg_levels=1` | corrigé (F1b) : axes sans épaisseur retirés avant la triangulation ; canal_turbulent_3d, `fmg_levels = 1`, `tol = 1e-10` : écart 3.7·10⁻⁸ avec le calcul sans multigrille, 2 629 itérations (+ 102 grossières) contre 3 321 ; test |

## P. Performances

| # | Constat | Statut |
|---|---|---|
| P1 | **Fils BLAS par défaut.** Un calcul « sur un cœur » occupe les 4 cœurs sans rien gagner : plaque compressible 10.2 s (CPU 39 s) contre 9.3 s (CPU 9 s) avec un fil BLAS. **Deux calculs simultanés s'effondrent** : 137.2 s au lieu de 9.7 s (×13) ; cavité cubique 3D 104.0 s au lieu de 19.5 s (×5). Les balayages `-j` sont déjà protégés (un fil par processus), pas un calcul seul (ligne de commande, interface, exécutable). Hypothèse à vérifier : le « Numba multi-fil plus lent » du lot D1 et une part des variations « ±50 % selon le jour » en viennent peut-être. | corrigé (F1) : un fil BLAS fixé avant l'import de NumPy (valeur de l'utilisateur conservée) ; aucun cas plus lent ; deux calculs simultanés 9.1 / 9.5 s au lieu de 137.2 s (`docs/multicoeur.md` § 3) |
| P2 | Interface 2 fois plus lente que la ligne de commande (cavité, 351 itérations : 10.2 s contre 4.9 s ; audit 1 : 8.6 / 3.5 s) ; cause non mesurée. | corrigé (F3). Cause mesurée : la courbe de convergence, retracée en entier à chaque progrès (toutes les 0.25 s) dans le fil de l'interface, qui garde le verrou de Python : 0.15 à 0.5 s par tracé, fil principal 30 s de processeur, calcul 34 à 41 s au lieu de 6.6 s (2026-10-07, 4 cœurs). Trouvé en mesurant : un résidu nul (Uy à la 1re itération) tracé à 1e-300, axe de 1e-314 à 1e14, courbes écrasées (interface **et** convergence.png de la ligne de commande). Correction : valeurs ≤ 0 non tracées ; retracé quand le temps écoulé atteint 10 fois le coût processeur du tracé (au plus 2 s), courbes mises à jour sans reconstruire la figure, ≤ 2000 points. Après : interface 8.9 à 9.6 s, sans courbe 8.7 à 8.8 s, ligne de commande 9.6 s (démarrage et figures compris) |

## M / U / L. Messages, interface, sorties (suite)

| # | Constat | Statut |
|---|---|---|
| M20 | Fichier de cas en UTF-8 avec BOM (Bloc-notes, PowerShell) refusé : « Erreur : Invalid statement (at line 1, column 1) » ; Latin-1 : « 'utf-8' codec can't decode byte 0xe9 in position 7: invalid continuation byte » ; sortie = fichier existant : « [Errno 17] File exists: '…' » ; Gmsh binaire : indicateur ASCII / binaire non lu (non vérifié : Gmsh absent). | corrigé (F2) : BOM accepté ; Latin-1 / Windows-1252 lu avec ATTENTION ; syntaxe TOML en français (ligne, colonne, rappels) ; « … existe déjà et n'est pas un dossier : choisir un autre dossier de sortie (-o) » ; Gmsh binaire (en fait : erreur interne KeyError) → « fichier Gmsh binaire, non lu. Le réexporter en texte » ; test |
| M21 | Une seule maille (2D 1 × 1, 3D 1 × 1 × 1) → « Erreur : Factor is exactly singular ». | corrigé (F2) : « Maillage de 1 cellule(s) sans face intérieure … au moins 2 cellules » avant le calcul ; « Factor is exactly singular » restant traduit avec des pistes ; test |
| M22 | Un dossier du dossier courant portant le nom d'un exemple masque l'exemple : `microrans run2d cavite_cubique_re100_3d` → « un fichier est attendu, pas un dossier » (il suffit d'un `-o cavite_cubique_re100_3d` précédent). | corrigé (F2) : l'exemple est pris, avec ATTENTION « … est aussi un dossier ici » ; dossier sans exemple : « est un dossier, pas un fichier de cas » ; test |
| U20 | « Arrêter » sans effet pendant le maillage (hybride : mené à son terme 21.8 s après la demande). | corrigé (F3) : le maillage non structuré teste l'arrêt à chaque itération (`microrans/stop.py`) ; hybride arrêté 0.07 à 1.5 s après la demande (au lieu de 21.8 à 39 s), maillage précédent conservé, pas de message d'échec |
| L7 | Les CSV pariétaux n'ont ni aire ni normale des faces : l'utilisateur ne peut pas refaire l'intégration des efforts (tableur, comparaison avec un autre code). De plus, le frottement du résumé intègre ν (U_P − U_paroi) / d **complet** (composante normale comprise, `solver.forces`), le CSV la seule composante tangentielle : à convergence, 0.195 % du C_d de frottement du cylindre Re 20 (≈ 0.08 % du C_d total), 0.006 % sur le NACA. Écart de discrétisation (la contrainte normale visqueuse est nulle à une paroi sans glissement) : projeter sur la tangente dans les efforts, ou l'écrire. | corrigé (F3) : colonnes p, area, nx, ny (, nz), tau_x, tau_y (, tau_z) dans tous les CSV pariétaux (incompressible 2D, axisymétrique — aire sur 360° —, 3D, compressible) ; Σ (p n + τ) × aire redonne C_d de pression et de frottement du résumé à 1e-10 près (test, 4 cas). Efforts inchangés (choix : ne pas changer les C_d publiés pour un écart de discrétisation de 0.08 %) ; tau_w = composante tangentielle de τ, écrit dans le README |
| D8 | README : « un calcul = un cœur » (faux en CPU consommé, P1) ; « cylindre Re = 100 : 2 à 16 min » : la commande RK3 + `adjust_dt` du § 3 a pris 1 061 s (17.7 min) sur machine libre (tableau : 949 s) ; « reprise … identique au bit près » (faux en compressible implicite, C17 ; vrai depuis F1b) ; un maillage importé ne peut pas avoir de frontières périodiques (`periodic` refusé pour `type = "file"`) : non documenté. | corrigé (F4) : « un calcul = un cœur » vrai depuis P1, chiffre ajouté (cavité : 6.6 s de processeur pour 6.7 s écoulées) ; durée du cylindre remesurée : 1 780 s le 2026-10-08 (1 061 s à l'audit), code d'avant F1 et actuel identiques côte à côte (160 s chacun, t_end = 10) : variation de la machine, README « 16 à 30 min » avec les trois mesures ; reprise exacte : déjà juste depuis F1b ; `periodic` des maillages importés (ignoré avec ATTENTION) : README § 8, limite 6 |

## Ce qui est vérifié juste (en plus de l'audit 2)

Interface = ligne de commande au bit près sur 20 des 22 exemples comparés ; reprise exacte au bit près
sur tout l'incompressible ; 2D = 3D une couche à la convergence ; export OpenFOAM
structurellement correct ; aller-retour .msh / .su2 exact ; maillage réimporté → même
résultat ; efforts de pression du résumé = intégration des CSV pariétaux ; CRLF, chemins accentués, JSON,
Gmsh 4.1 ; exécutable = Python (résultats et vitesse) ; pas de fuite mémoire de l'interface ;
plages de balayage mal formées refusées clairement.

## Lots de correction, mis à jour (dans cet ordre, avant l'import Gmsh 3D)

- **F1 — résultats faux, calculs perdus** : C11, C13, C14, C15 (convention à arrêter), C18,
  C21, U14, U17, P1.
- **F1b — plantages et reprises** : C16, C17, C19, C20, C22, C23 (fait, 2026-10-05 ; C20 avancé en F1).
- **F2 — messages** : M13 à M22, U16, L6 (fait, 2026-10-06).
- **F3 — interface 3D, figures, sorties** : U15, U18, U19, U20, L5, L7, P2 (fait, 2026-10-07 ; trouvé en mesurant P2 : axe des résidus sur 330 décades, corrigé).
- **F4 — textes et documentation** : T1, T2, D5 à D8 (fait, 2026-10-08).
