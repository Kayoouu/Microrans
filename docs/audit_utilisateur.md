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
| U5 | Maillage en O : « Ajouter » crée un 2e corps interdit ; tout corps ajouté est posé exactement sur le premier. | à faire |
| U6 | Accueil : 22 exemples non classés (noms de fichiers, descriptions coupées), exemples de maillage seul mélangés aux calculs, rien pour dire par où commencer. | à faire |
| U7 | Résumé des résultats en JSON brut (16 chiffres). | à faire |
| U8 | Textes périmés : accueil et « À propos » (« RANS/URANS », 4 modèles ; ni compressible, ni couplé, ni transition) ; « Stationnaire (RANS) » même en laminaire. | à faire |
| U9 | Lois de paroi proposées avec k-ε et transition : refus seulement au lancement. | à faire |
| U10 | Balayage : valeurs mal saisies (« 0 à 4 ») → message Python en anglais (« could not convert string to float »). | à faire |
| U11 | *(trouvé pendant le lot 3)* Fluide non newtonien sans ν dans le cas (ν de référence tiré de la loi, ex. sang ≈ 3.6·10⁻⁵ m²/s) : l'interface affichait ν = 0.01 et l'**écrivait dans le cas** à la première modification, ~280 × la valeur de la loi (Re affiché, diffusivité thermique α = ν/Pr). | corrigé : champ ν vide permis, « auto : ν de la loi à γ̇ = U/L » ; vérifié sur les 22 exemples que l'interface n'écrit plus que des valeurs égales aux défauts du solveur |
| U12 | *(trouvé pendant le lot 3)* Champs numériques arrondis à 6 chiffres (ν = 1/550 → 0.00181818, t_end = 0.000632455532 → 0.000632456) puis réécrits dans le cas. Écart ~10⁻⁶, sans effet visible, mais valeur modifiée sans action de l'utilisateur. | corrigé : affichage court quand il est exact (3000, 0.1), sinon tous les chiffres |

## L. Ligne de commande

| # | Constat | Statut |
|---|---|---|
| L1 | Description « RANS/URANS » incomplète ; pas de `--version` ; aide mêlant anglais et français. | **en partie** : `microrans --version` ajouté ; description et langue de l'aide : lot 4 |
| L2 | Fin de calcul : bloc JSON de 40 lignes au lieu d'un résumé lisible. | à faire |
| L3 | `--set` : deux niveaux seulement (`boundary.lid.U=…` plante). | **corrigé** : `--set` à plusieurs niveaux (`boundary.lid.U=[2,0]`), indice de liste (`bodies.0.radius=0.3`), virgule décimale ; `--set` sans `=` ou clé incomplète → message clair (`tests/test_cli.py`) |

## D. Documentation

| # | Constat | Statut |
|---|---|---|
| D1 | Pas de tutoriel pas à pas pour un premier calcul. | à faire |
| D2 | Pas de section dépannage (divergence, « non convergé », y⁺, qualité du maillage). | à faire |
| D3 | Pas de glossaire (RANS, SIMPLE, y⁺, Cp, Cf, patch, O-grid…). | à faire |
| D4 | Pas de référence complète des clés du fichier de cas ; l'extrait du README oublie `sst_gamma`. | à faire |

## Observations à approfondir (non traitées)

- **Performance** : sur la cavité (laminaire, SIMPLEC), le temps par itération passe de
  0.036 s (10 000 cellules) à 0.79 s (40 000 cellules), soit ×5 par cellule, au seuil de
  20 000 cellules où le solveur de pression passe de la factorisation directe (LU) à l'AMG.
  Le seuil ou le réglage de l'AMG est peut-être mal placé ; à mesurer avant de conclure.
- **Corps superposés** (maillage hybride) : « Maillage non manifold : une arête est partagée
  par plus de 2 cellules » après ~20 s de maillage. Le chevauchement devrait être détecté
  avant de mailler (lié à U5).
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
