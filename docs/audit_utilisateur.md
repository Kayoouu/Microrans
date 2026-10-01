# Audit utilisateur (débutant et expert) — constats et suivi

Campagnes automatiques, toutes reproductibles, menées comme le ferait un utilisateur :

| Campagne | Contenu | Résultat |
|---|---|---|
| 1. Exemples, ligne de commande | les 22 exemples, `microrans run2d <exemple>`, en entier, figures comprises | 20 terminés et convergés ; 2 échecs (exemples de maillage seul) |
| 2. Exemples, interface | ouverture, maillage, calcul court, tracé de chaque champ, distributions pariétales, profil, historique, enregistrement, réouverture | 18 / 22 sans erreur ; 4 plantages (sondes) |
| 3. Entrées invalides, ligne de commande | 54 cas : fautes de frappe, valeurs impossibles, types faux, fichiers absents, TOML mal formé, `--set` | 13 traces Python brutes ; ~15 erreurs acceptées sans avertissement |
| 4. Combinaisons, interface | 115 : 6 types de maillage × 5 types de corps, 6 modèles × 2 traitements pariétaux, 3 algorithmes, 9 schémas en temps, 6 conditions limites, thermique, axisymétrique, 6 lois non newtoniennes, options, TOML, reprise, arrêt, export, balayage, 1D | 75 sans erreur ; 40 échecs (7 causes) |
| 5. Relecture | README, textes de l'interface (captures de chaque page), aide de la ligne de commande, saisie clavier en français | voir U et D ci-dessous |

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
| C5 | **Valeurs impossibles acceptées sans avertissement** : ν < 0, x1 < x0, rayon < 0, tolérance < 0, Mach < 0, vitesse à 3 composantes, condition pour une frontière inexistante, `nu` et `reynolds` donnés ensemble (l'un est ignoré). | fichier modifié | **en partie** : corps et `[domain]` vérifiés (rayon, demi-axes, corde > 0, x0 < x1, y0 < y1 ; `test_shape_errors_are_explained`) ; le reste à faire (lot 2) |
| C6 | **Balayage du nombre de Reynolds sans effet sur un cas donné en ν** (trouvé en écrivant la vérification des clés) : `nu` étant prioritaire, tous les points étaient calculés avec le même ν (cavité, Re = 100 et 1000 : Cd identiques à 6 chiffres), sans message. Touchait l'interface (page Balayage) et `microrans sweep`. | balayage `physics.reynolds` de `cavite_re100` | **corrigé** : balayer `reynolds` retire `nu` du point (et inversement) (`test_reynolds_sweep_on_case_given_in_nu`) |

## M. Messages d'erreur

| # | Constat | Statut |
|---|---|---|
| M1 | 13 erreurs affichent une **trace Python brute** : `nx` non entier ou texte, `max_iter` texte, `relax_U = 0`, `mode = "unsteady"`, instationnaire sans `dt`, `dt = 0`, `dt < 0` (après 9 s de calcul), fichier de maillage absent, corps sans `type`, `--set` à 3 niveaux, `--set` sans `=`, `--set` non numérique. | **en partie** : corps sans `type` → message clair ; les 12 autres à faire (lot 2) |
| M2 | **Messages trompeurs** : `nx = 0` → « Trop de segments de progression » ; `ν = 0` → « Factor is exactly singular » ; `[mesh]` absente → « Conditions aux limites manquantes pour bottom, inlet… » ; fichier vide → « donner nu ou reynolds » ; incidence « dix » → message Python en anglais ; `U = "vite"` → « nom inconnu v ». | à faire |
| M3 | L'interface affiche l'exception brute (« KeyError : 'type' ») sans explication ni renvoi au journal. | **corrigé** : fenêtre « Le calcul n'a pas pu aboutir », message en clair, renvoi au journal ; erreur interne signalée comme telle ; titre selon la tâche (« Le maillage n'a pas pu aboutir »…) (`test_error_dialog_text_is_readable`, `test_mesh_failure_dialog`). Trouvé en corrigeant : un corps sans `type` saisi dans l'onglet TOML faisait planter la liste des corps ; corrigé |
| M4 | Exemples de maillage seul (`mesh_cylindre_hybride`, `mesh_naca_multi`) lancés comme un calcul : « donner nu ou reynolds » **après 89 à 116 s de maillage**, sans dire qu'il s'agit d'un exemple de maillage. | à faire |
| M5 | Maillage de 4 millions de cellules : aucun retour pendant 3 minutes, aucune estimation de mémoire ou de durée. | à faire |
| M6 | Cartes graphiques proposées dans l'exécutable alors qu'elles ne peuvent pas y fonctionner ; le message conseille `pip install`. | **corrigé** : choix grisés dans l'exécutable avec explication ; backend inconnu refusé avec la liste (`test_gpu_backends_explained_in_executable`) |
| M7 | Maillage « fichier » sans chemin (« Format non supporté : (msh, su2) », puis `KeyError 'path'`), maillage « multi-blocs » vide (`KeyError 'vertices'`), corps « contour importé » sans fichier (`IsADirectoryError`). | **en partie** : contour sans fichier ou introuvable → message clair ; « fichier » et « multi-blocs » à faire (lot 2) |
| M8 | `"strouhal": NaN` dans `summary.json` pour un calcul instationnaire court (JSON invalide pour d'autres outils). | à faire |

Messages déjà bons (à garder comme modèle) : modèle, condition, schéma, matériel ou type
de maillage inconnus (liste des choix) ; expression dangereuse refusée ; condition limite
manquante ; axisymétrique hors du demi-plan y ≥ 0 ; lois de paroi avec k-ε ou transition ;
reprise introuvable ; TOML mal formé (ligne et colonne).

## U. Interface

| # | Constat | Statut |
|---|---|---|
| U1 | **Tableau des conditions limites illisible** : colonnes écrasées (« cyli… », type « Pa », « Ch »). | à faire |
| U2 | **Colonne centrale trop étroite** : champs et notes coupés, défilement horizontal sur toutes les pages. | à faire |
| U3 | Mode « Nombre de Reynolds » : le champ ν grisé affiche une valeur périmée (0.01 au lieu de 0.05 pour Re = 20). | à faire |
| U4 | Non newtonien : 9 paramètres affichés « défaut » alors qu'aucun défaut n'existe → erreur « paramètres manquants » au lancement. | à faire |
| U5 | Maillage en O : « Ajouter » crée un 2e corps interdit ; tout corps ajouté est posé exactement sur le premier. | à faire |
| U6 | Accueil : 22 exemples non classés (noms de fichiers, descriptions coupées), exemples de maillage seul mélangés aux calculs, rien pour dire par où commencer. | à faire |
| U7 | Résumé des résultats en JSON brut (16 chiffres). | à faire |
| U8 | Textes périmés : accueil et « À propos » (« RANS/URANS », 4 modèles ; ni compressible, ni couplé, ni transition) ; « Stationnaire (RANS) » même en laminaire. | à faire |
| U9 | Lois de paroi proposées avec k-ε et transition : refus seulement au lancement. | à faire |

## L. Ligne de commande

| # | Constat | Statut |
|---|---|---|
| L1 | Description « RANS/URANS » incomplète ; pas de `--version` ; aide mêlant anglais et français. | à faire |
| L2 | Fin de calcul : bloc JSON de 40 lignes au lieu d'un résumé lisible. | à faire |
| L3 | `--set` : deux niveaux seulement (`boundary.lid.U=…` plante). | à faire |

## D. Documentation

| # | Constat | Statut |
|---|---|---|
| D1 | Pas de tutoriel pas à pas pour un premier calcul. | à faire |
| D2 | Pas de section dépannage (divergence, « non convergé », y⁺, qualité du maillage). | à faire |
| D3 | Pas de glossaire (RANS, SIMPLE, y⁺, Cp, Cf, patch, O-grid…). | à faire |
| D4 | Pas de référence complète des clés du fichier de cas ; l'extrait du README oublie `sst_gamma`. | à faire |

## Ce qui fonctionne (vérifié)

Tous les modèles de turbulence (laminaire, SA, k-ε, k-ω, SST, SST-γ), les trois
algorithmes stationnaires, les neuf schémas en temps, le pas adaptatif, les six types de
conditions limites sur la cavité, la thermique, les scalaires, les zones poreuses, les
disques actuateurs, le démarrage multigrille, le pseudo-transitoire, les sondes (saisies
dans l'interface), l'arrêt sur efforts, la reprise, l'arrêt en cours de calcul, l'export du
maillage, le balayage, le canal 1D (RANS et URANS), les 20 exemples de calcul en ligne de
commande. Aucune valeur NaN ou aberrante dans les 102 résumés produits (hors M8).
