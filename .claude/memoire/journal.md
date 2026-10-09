# Journal des lots (le plus récent en haut)

Détail des chiffres : README (§ 6 validation, § 7 performances, § 8 limites) et messages de
commit. Ce journal sert à retrouver ce qui a été fait, pourquoi, et ce qui a été constaté.

## 2026-10-09 — Lot GCI : convergence en maillage (jalon A, point 1)

- Relance quotidienne (03:53 UTC). Vendredi, campagnes de moins de 7 jours : pas de
  campagne. CI verte au départ (cea2d8a).
- `microrans/gci.py` : procédure de Celik et al. (2008) (ordre apparent par point fixe,
  Richardson, GCI facteur 1.25, convergence oscillante signalée) ; `tests/test_gci.py` :
  exemple chiffré de l'article (p 1.53, φ_ext 6.1685, GCI 2.2 %), données exactes d'ordre 2
  en 2D / 3D, entrées refusées. Pas de preuve « échoue avant » possible (module nouveau).
- `tools/validation/gci_maillage.py` (12 calculs, cavité 128² la plus longue : 99 s) :
  cylindre Re = 20 C_d 2.0573 / 2.0374 / 2.0327, p 2.08, extrapolé 2.0312, GCI 0.09 % ;
  cavité u min −0.21151 / −0.21351 / −0.21387, p 2.47, extrapolé −0.21395, GCI 0.05 % ;
  conduite 3D 3.5670 / 3.5276 / 3.5177, p 1.99, extrapolé 3.514401 contre 3.514425 exact.
- Erreurs de ma part corrigées en cours : conduite d'abord avec r = (N1/N2)^(1/3) alors que
  seuls y et z sont raffinés (p 2.98 faux → dim = 2, p 1.99) ; tol 1e-11 jamais atteinte en
  3D (plancher ~1e-9) → tol 1e-9, erreur d'itération mesurée 2e-7.
- Trouvé : l'accord du README « C_d 2.037 contre 2.045, 0.6 % » dépend du rayon du domaine :
  R = 20 / 40 / 80 → 2.0801 / 2.0374 / 2.0195 ; extrapolé R → ∞ ≈ 2.007 (p 1.25,
  estimation). Écrit dans README § 6 (ligne du cylindre, nouvelle section « Incertitude de
  maillage »), limite 7, tutoriel § 6 (« Le domaine est-il assez grand ? », commandes
  vérifiées : 2.08 et 2.02), en-tête de l'exemple, glossaire (GCI). Cavité : u min s'éloigne
  de Ghia en raffinant (au point y = 0.4531 : −0.21095 / −0.21335 / −0.21383 contre −0.21090) ;
  impossible de trancher sans référence plus fine vérifiée (recherche en ligne : Fornberg
  1980 et valeurs spectrales de la cavité non vérifiables depuis l'environnement → non cités).

## 2026-10-08 — Lot F4 : textes et documentation (T1, T2, D5 à D8)

- Relance quotidienne (03:53 UTC). Jeudi : pas de campagnes. CI verte au départ (995af3e).
- T1 : « 2D et 3D » dans `microrans --help`, « À propos » (texte sorti dans
  `gui.app.about_text()`, paragraphe 3D), pyproject, description du fichier Windows
  (`packaging/version_info.txt`), docstring du paquet. T2 : message d'ouverture d'un cas 3D.
- D5 : vérifié d'abord : faces z en symétrie = `[boundary.back]` / `[boundary.front]`
  `type = "symmetry"` suffit ; `patch_types` seul → « Conditions aux limites manquantes ».
  README § 3 : extrait complet (cylindre extrudé) + variante symétrie, calculés par un test.
- D6 : tutoriel § 7 « Un cas 3D » (interface, ligne de commande, extrusion du cylindre).
  Mesuré : conduite 3.528 (+0.38 %, 406 it, 7.3 s), 16² : 3.567 (+1.50 %) ; cylindre
  extrudé 1 couche entre symétries : C_d 2.037 en 183 it (2D : 182 it) ; `fields.vtk` lu
  et coupé (x = 0.25) par la bibliothèque VTK 9.7 (U_x max 7.36) — ParaView lui-même non
  essayé (absent de la machine).
- D7 : glossaire (pavé / extrusion, hexaèdre / prisme, périodiques, plan de coupe) ;
  dépannage : 3 lignes (type inconnu, periodic sans effet, slice en 2D) + figures omises,
  extrusion impossible ; test : extraits cités = vrais messages.
- D8 : « un calcul = un cœur » vrai depuis P1 (cavité : 6.6 s CPU / 6.7 s mur). Cylindre
  RK3 + adjust_dt du § 3 : **1 780 s** aujourd'hui (audit 2 : 1 061 s, tableau : 949 s) ;
  ab.py temps e5deb38 (avant F1) contre l'arbre, t_end = 10 : 160.2 / 160.5 s, CPU = mur
  dans les deux cas (BLAS non parallèle sur ce cas, P1 sans effet) → pas de régression,
  variation de la machine (la cavité, elle, était plus rapide qu'hier). README : « 16 à
  30 min » avec les trois mesures. periodic d'un maillage importé : limite 6.
- 5 tests (`tests/test_docs_f4.py`), 5 « échoue avant : OK ».

## 2026-10-07 — Lot F3 : interface 3D, figures, sorties (P2, L7, L5, U15, U18, U19, U20)

- Relance quotidienne (03:53 UTC). Mercredi : pas de campagnes. CI verte au départ (53acaf5).
- P2 mesuré d'abord (cavité, défaut, 351 it, 4 cœurs, A/B alternés) : interface 34 à 46 s,
  sans la courbe de convergence 8.1 à 9.1 s, ligne de commande 9.6 s (démarrage et figures
  compris ; solveur 5.9 à 6.2 s). Fil de l'interface : 30 s de processeur contre 0.2 s sans
  courbe → la courbe (figure reconstruite toutes les 0.25 s, 0.15 à 0.5 s par tracé, verrou
  de Python tenu) bloquait le calcul. Trouvé en mesurant : un résidu nul tracé à 1e-300,
  axe de 1e-314 à 1e14 (convergence.png de la ligne de commande aussi : figure de la cavité
  illisible depuis toujours). Corrigé : `postprocess.log_values` (≤ 0 → non tracé, 4
  endroits) ; retracé quand le temps écoulé atteint 10 × le coût processeur du tracé (≤ 2 s),
  courbes mises à jour sans reconstruire la figure, ≤ 2000 points. Après : 8.9 à 9.6 s avec
  courbe contre 8.7 à 8.8 s sans.
- L7 : colonnes p, area, nx, ny (, nz), tau_x, tau_y (, tau_z) dans tous les CSV
  pariétaux (aire sur 360° en axisymétrique) ; Σ (p n + τ) A redonne C_d de pression et de
  frottement du résumé à 1e-10 (cylindre, sphère axisymétrique, conduite 3D, plaque
  compressible). Efforts non changés (choix : 0.08 % du C_d du cylindre, écart de
  discrétisation qui tend vers 0 ; ne pas changer les chiffres publiés pour cela).
- L5 : `[output] slice_axis` / `slice_value` (vérifiés avant calcul ; hors domaine :
  ATTENTION et figures omises) ; conduite carrée : section x ; noms d'axes ; figure moins
  large pour un domaine haut (1500 → 840 px).
- U20 : `microrans/stop.py` (fonction de test installée par la tâche de l'interface,
  `check_stop()` dans les boucles de `triangulate`) ; arrêt du maillage hybride 0.07 à
  1.5 s après la demande (avant 21.8 à 39 s), signal `stopped` (pas de message d'échec).
- U15 : sondes / lignes à 3 composantes (z = milieu) à l'extrusion et retour ; extrusion
  refusée tout de suite avec des options 2D seulement ; extrusion proposée en multi-blocs.
  U19 : profil par défaut à travers le domaine, bouton animation, zoom grisé en coupe x / y,
  légende de la vue 3D sous la vue, types de frontière en clair (`patch_label`). U18 :
  fichiers relatifs copiés à l'« Enregistrer sous » ailleurs, reprise en chemin absolu.
- Outil : `ab.py` compare les colonnes communes quand des colonnes sont ajoutées (avant :
  « DIFFÉRENT » sans détail) ; test ajouté.
- Vérifié : 8 tests (`tests/test_sorties_f3.py`), 8 « échoue avant : OK » ; `ab.py egalite`
  sur 6 exemples (cylindre, sphère axisymétrique, conduite 3D, plaque compressible, cavité,
  NACA SA) : tout identique au bit près sauf les colonnes ajoutées.

## 2026-10-06 — Lot F2 : messages (M13 à M22, U16, L6)

- Relance quotidienne (03:53 UTC), premier lot suivant le skill `lot`. Mardi : pas de
  campagnes de non-régression.
- Reproduction d'abord : `tools/audit/c16_messages_f2.py` (campagne 16, 17 cas). Constats
  en plus de l'audit : Gmsh binaire → « Erreur interne (KeyError 'Nodes') » (l'audit
  supposait un message obscur) ; clés compressibles en texte → erreur interne
  (UFuncTypeError).
- Corrections : `tomlio.read_text` (BOM, repli Windows-1252 avec ATTENTION) et `loads`
  (syntaxe TOML en français), utilisés par `load_config` et l'interface ;
  `mesh2d.mesh.check_patch_names` / `unknown_patch_type` (2D et 3D) ;
  `validate.missing_bc_message` (incompressible et compressible ; faces d'extrusion :
  symmetry ou [mesh.extrude] periodic, vérifié) ; `_Check._home` (clé mal placée) ;
  clés expert de [solver] vérifiées (nombres, entiers, bornes, choix, booléens) ;
  `build_solver` refuse 0 face intérieure ; `_error_text` (fichier existant, pas un
  dossier, « Factor is exactly singular ») ; `resolve_example` (dossier homonyme) ;
  `parse_points` (virgule décimale, U16) et sondes vérifiées avant le lancement ;
  résumé : bruit d'arrondi de U moyen → 0 (L6) ; trailing_edge hors NACA → ATTENTION ;
  `--steps-per-period` dans le message de la phase ; en-tête de mesh_naca_multi.
- Non fait (choix) : abréviations d'options (argparse) laissées actives ; seul le message
  de M19 change.
- Vérifié : 12 tests (`tests/test_messages.py`), 12 « échoue avant : OK » (echoue_avant,
  aucune preuve faible) ; `ab.py egalite` identique au bit près sur cavité, rampe M = 2,
  conduite 3D, filtre poreux, cylindre URANS (sondes) ; les 23 exemples de calcul passent
  la vérification ; campagne 16 après : chaque cas a un message en français avec choix
  ou suggestion.
- Campagnes relancées : 3c (12 clés « expert » dans un cas où elles servent, texte puis
  −1) : 24 / 24 erreurs claires avant calcul (audit : erreurs internes et acceptations
  silencieuses) ; 3b (toutes les clés numériques) : 88 erreurs claires, 4 acceptées en
  silence (dt, t_end invalides dans un cas stationnaire, où ils ne servent pas) →
  désormais vérifiés aussi en stationnaire (second commit, test complété) ; 3b refaite :
  92 / 92 erreurs claires.
- CI de a063d20 verte (tests 37412579973, exécutables 37412583026) ; suite après le second
  commit : 489 réussis, 2 ignorés.

## 2026-10-05 — Outils et automatismes de travail (demande de l'utilisateur)

- Demande : « améliore autonomement ton workflow », avant les lots quotidiens.
- Constats tirés des derniers lots : preuve « le test échoue avant » faite à la main
  (`git stash`, risque de perdre du travail) ; comparaisons avant / après et mesures
  (getrusage) réécrites à chaque fois ; suite lancée pendant l'édition de la doc (course
  avec test_reference_document_up_to_date) ; attente de la CI par boucles de sleep + MCP ;
  règles de CLAUDE.md (ruff, suite, pas de force, pas de nom de modèle, tomllib) vérifiées
  de mémoire ; aucune mesure estimé / réel alors que des délais sont annoncés.
- Fait : `tools/dev/` (suite.py, echoue_avant.py, ab.py egalite / temps, ci.py,
  avant_push.py, _commun.py) ; hook PreToolUse (Bash) avant poussée ; hook Stop étendu
  (travail non commité / non poussé) ; état de la CI au démarrage (SessionStart) ; skill
  `lot` (procédure) ; `bilan_lots.md` (estimé / réel, indicateur du jalon A : résultats
  faux silencieux par audit 8 → 3 → 2, campagnes de non-régression) ; routine mise à jour
  (suit le skill ; campagnes le lundi). Accès API : add_repo Kayoouu/Microrans (push), car
  `gh api` refusait le dépôt (ancien nom → chemin numérique refusé par le proxy).
- Essais : echoue_avant sur les tests de F1b contre 898cffa~1 : 3 « échoue avant : OK »,
  C17 « preuve faible » (KeyError avant l'égalité des champs : juste), test FMG préexistant
  « ne prouve pas » (juste). ab egalite : cavité identique ; conduite 3D avant / après C15 :
  coefficients ×9, champs identiques (déjà mesuré à la main). ab temps : cavité 200 it,
  B / A = 1.007. ci.py : état et artefacts du run 37263377558 relus.
- Défauts trouvés en essayant : ab avalait --ref avec la commande, et deux exécutions en
  erreur donnaient « identique » (corrigés) ; le hook prenait « git push » cité dans un
  texte pour une poussée et bloquait sa propre correction (motif limité à la position de
  commande ; correction faite avec l'outil Edit). Tests : tests/test_outils_dev.py (4).

## 2026-10-05 — Jalons A, B, C inscrits dans a_faire.md

- Question de l'utilisateur (2026-10-04) : délai pour un usage académique / pré-industriel.
  Réponse : A (2D académique crédible) fin octobre à mi-novembre 2026, B (3D académique,
  petite géométrie) janvier–février 2027, C (pré-industriel) pas atteignable avec
  l'architecture actuelle. Estimations, pas des mesures. Accord pour les inscrire :
  « oui écris ».
- Réordonnancement qui en découle : GCI et maillage en C + NASA TMR (jalon A) passent
  avant l'import Gmsh 3D et le multi-cœur (jalon B) ; ajout de l'audit 3 (critère du
  jalon A : plus de résultat faux silencieux) et d'un cas 3D turbulent comparé à des
  données publiées (jalon B). Compressible turbulent : seulement sur demande.
- Restes de l'estimation A après C15 et F1b (faits en 2 lots, estimés 2 à 3) : 9 à 16 lots.

## 2026-10-05 — Lot F1b : plantages et reprises

- Relance quotidienne (03:53 UTC). Un test par point, chacun vérifié en échec sur l'ancien
  code (`git stash push <fichiers>`).
- C16 : `perturbation_center` en 3D : [x, y] = tube selon z, ou [x, y, z] ; autre taille
  refusée (case.py, doc de la clé, reference_cas.md régénéré).
- C22 : `_freestream_change.rot` (sweep.py) ne tourne que les 2 premières composantes.
- C23 : `restart._interpolator` retire les axes d'épaisseur ≤ 1e-9 × taille avant
  Delaunay (couche unique au niveau grossier). Canal 3D, fmg_levels = 1, tol 1e-10 :
  écart 3.7e-8 avec le calcul sans FMG ; 2 629 it (+102 grossières) contre 3 321.
- C19 : cli.py `_ctrl_c_stops_cleanly` : 1er Ctrl-C = drapeau d'arrêt (fin d'itération,
  tous les fichiers, code 130), 2e = KeyboardInterrupt ; workers de balayage ignorent
  SIGINT. Campagne 11 : 3/3 cas avec checkpoint, résumé, champs, historique, CSV paroi,
  reprise ok ; balayage parallèle arrêté en 0.4 s, 2 points écrits. Windows non vérifié
  (test ignoré sous win32).
- C17 : état du pilotage gardé par le solveur (`_ctrl` : cfl, cfl_cap, prev, best, since,
  cuts, path, q_start, derniers relevés d'efforts), `_res0_count`, `_psi_frozen`, Γ du
  tourbillon que l'itération suivante lirait (calculé avant le résidu final, ou à la
  sauvegarde automatique) → `restart_state()` / `set_restart_state()` dans le checkpoint
  (méta JSON + tableaux `state_*`). Relevés monitor / monitor_tol sur `gi % 10` (avant
  `it % 10`, décalés à la reprise). Restauré seulement si stationnaire → stationnaire et
  signature (réglages hors max_iter/tol/monitor_tol/log_every, gaz, amont) inchangée ;
  sinon pilotage neuf (`restart.controller` dans summary.json). Vérifié : N + N = 2N au
  bit près (NACA 120+120, 90+90 fenêtre 50, arrêt sur efforts à 240 repris à 205, plaque
  limiter_freeze 20, rampe) ; calcul continu = ancien code au bit près (5 variantes) ;
  campagne 7 : 18/18 exemples exacts (Sod non découpable).
- Suite : 474 réussis, 1 ignoré (6 min 30 s). README (reprise, Ctrl-C), dépannage § 5,
  audit (statuts) mis à jour.

## 2026-10-04 — C15 : une seule vitesse de référence

- Demande : « La meilleure solution pour C15 » (choix délégué). Mesuré d'abord : 4 usages,
  2 conventions (1 pour ν et γ̇_ref ; vitesse imposée max OU vitesse initiale pour les
  coefficients) ; interface : 1 écrit. Seuls 2 exemples sans `reference_velocity` : cavité
  cubique (auto 1 = 1) et conduite carrée (auto 3 = **vitesse initiale**, arbitraire).
- Règle retenue (`choose_reference_velocity`, solver.py, utilisée par le solveur et par
  build_solver avant ν) : donnée > vitesse d'entrée (inlet U ou débit, farfield ;
  moyenne pondérée par l'aire, ×r en axisymétrique ; la plus grande + avertissement) >
  paroi mobile > 1. Entrée prioritaire sur la paroi : cylindre tournant → Re sur U∞.
  Reynolds sans vitesse imposée : avertissement. U_ref et origine dans le journal et
  summary.json (`reference_velocity_source`). Interface : champ vide = auto ; affichage
  de ν / Re dérivé avec la même règle.
- Vérifié : reynolds = 100 + entrée U = 2 → ν = 0.02 (avant 0.01, Re 200) ; campagne 6 sur
  22 exemples : interface = ligne de commande partout (hors chemins) ; ligne de commande
  avant / après : seuls les coefficients de la conduite changent (×1/9 ; champs identiques
  au bit près). T3A : C_l (≈ 5e-4) change de 1.7e-6 relatif entre l'audit 2 et maintenant :
  dû à P1 (ordre des sommes BLAS), identique au bit près avant / après C15 à BLAS égal.
- Tests : 3 nouveaux ; suite 468 réussis, 1 ignoré (5 min 35 s).

## 2026-10-04 — Lot F1 : résultats faux, calculs perdus

- Relance quotidienne (03:53 UTC). C15 laissé de côté : décision de convention posée à
  l'utilisateur, pas de réponse (le silence ne vaut pas accord).
- P1 : `microrans/__init__.py` fixe OMP / OPENBLAS / MKL / VECLIB à 1 avant l'import de
  NumPy (valeur utilisateur conservée, NUMBA_NUM_THREADS non touché : Numba garde 4 fils).
  Mesuré A/B/B/A sur 6 cas (cavité 2D, plaque compressible, cavité 32³, NACA SA 200 it,
  cylindre URANS t = 10, cavité 64³ 10 it) : aucun cas plus lent avec 1 fil ; CPU ÷ 4 sur
  plaque et cavité 32³. Deux calculs simultanés : plaque 9.1 / 9.5 s (avant 137.2 s),
  cavité 32³ 17.9 / 18.6 s (avant 104.0 s).
- C11 : dossiers de sortie `.resolve()` (case.py, compressible_case.py, sweep.py).
- C13 : `_load_extrude` compare aux noms `names.back/front`.
- C14 : `check_force_patches` avant les itérations (refus en 0.3 s).
- C18 : `_numbers` : « ; » (ou tabulation sans point) → virgule décimale. Première version
  (découpage sur les seules tabulations) écartée : une ligne .dat « 1.0 \t 2.0  3.0 »
  aurait été ignorée.
- C20 (avancé de F1b) : valeurs entières pour les clés de maillage balayées.
- C21 : `check_sweep_key` (nouvelle clé comparée aux avertissements de check_case) :
  refuse aussi une clé sans effet (mesh.nx sur un maillage en O, solver.cfl en
  incompressible) ; `bodies.0.radius` : refus clair (set_key ne traverse pas les listes).
- U14 : signature [mesh] / [domain] / [[bodies]] mémorisée au maillage ; remaillage annoncé.
- U17 : première version (flèches à l'échelle de |U| complet) écartée après avoir regardé
  la figure : l'écoulement secondaire réel de la cavité cubique (plan x = 0.5 : 11 % de
  |U|) devenait invisible. Retenu : rien sous 0.1 % de |U|, sinon échelle automatique +
  « flèches agrandies (max 11 % de |U|) » dans le titre.
- Tests : 11 nouveaux, chacun vérifié en échec sur l'ancien code ; suite 465 réussis,
  1 ignoré (5 min 33 s). Campagne 1 : 25 / 25 (avant 24 / 25).
- Campagnes 2 et 4 relancées : identiques à l'audit 2, sauf « extrusion cochée / décochée
  après maillage » (1 erreur → 0, effet de U14). Campagne 6 : canal interface = ligne de
  commande au bit près ; conduite toujours ×9 (C15). Exécutable local : 1 fil BLAS
  effectif. CI : tests 37177299915, executables 37177304448 (verts).

## 2026-10-04 — Audit 2 approfondi : cohérence des résultats, fichiers, arrêts

- Demande : « Faire ce qui a déjà était fait est quand même une bonne ideee mais n'hésite
  pas à approfondir l'audit ».
- Méthode : chercher les résultats faux par des contrôles sans référence extérieure (deux
  chemins qui doivent donner la même chose), puis les usages réels. Campagnes 6 à 15 dans
  `tools/audit/` (c6 à c15 ; c1 accepte `MICRORANS_EXE` pour l'exécutable).
- Constats les plus graves : C15 (vitesse de référence : C_d × 9 entre interface et ligne
  de commande sur conduite_carree_3d ; `reynolds = 100` + entrée U = 2 calculé à Re = 200
  sans avertissement), C18 (contour CSV « x;y » virgule décimale → 14 Go, processus tué),
  C21 (balayage d'une clé mal écrite accepté), C19 (Ctrl-C : rien d'écrit), P1 (BLAS
  4 fils : 2 calculs simultanés 137.2 s au lieu de 9.7 s).
- Vérifié juste : interface = ligne de commande au bit près (21 / 23) ; reprise exacte au bit
  près en incompressible et RK3 (pas en compressible implicite / NS : C17, convergé
  identique mais +14 % d'itérations) ; 2D = 3D une couche à 2e-6 près ; OpenFOAM
  structurellement correct ; .msh / .su2 aller-retour exact ; pression = intégration des
  CSV (frottement : 0.195 % d'écart sur le cylindre, composante normale incluse dans le
  résumé, ajouté à L7) ;
  exécutable = Python (24 / 25, même vitesse, pic mémoire 113 contre 158 Mo).
- README : la commande RK3 + adjust_dt du cylindre Re 100 a pris 1 061 s (17.7 min) sur
  machine libre, au-dessus des « 2 à 16 min » écrits (ajouté à D8).
- Pièges : un essai de contour CSV français sans limite mémoire a fait tuer un processus à
  14 Go (les scripts limitent maintenant à 4 Go) ; deux calculs simultanés faussent tous
  les chronométrages (P1) : un seul calcul à la fois pour mesurer.
- Lots F1 (avec P1, C15, C18, C21 en tête), F1b (plantages, reprises) ajoutés ; F2 à F4
  complétés ; la mesure « Numba multi-fil plus lent » est à refaire après P1.

## 2026-10-04 — Audit 2 : fouille bugs et expérience utilisateur

- Demande : « Avant de continuer le lot 3D, je pense que ce serait bien de faire une session
  fouille pour ses bugs, pour l'expérience utilisateur, ect ect ».
- 5 campagnes (scripts versionnés dans `tools/audit/`, sorties en scratchpad) :
  1. 25 exemples, commande d'en-tête copiée, dossier vide : 24 / 25 ; la polaire échoue au
     2e point (C11) ; durées conformes (canal 3D 23.6 s au 1er lancement = cache des
     polices de Matplotlib, 10.6 s ensuite).
  2. 25 exemples dans l'interface : aucun plantage ; canal_turbulent_3d calculé faux (C13 :
     C_d du fond 0.00791 au lieu de 0.008889, faces z devenues parois).
  3. 44 variantes invalides (36 correctes, 8 à reprendre) + 116 valeurs sur les clés
     numériques (M18 : clés expert non vérifiées, « erreur interne »).
  4. 38 scénarios 3D dans l'interface (+ 5 à part) : maillage périmé utilisé (U14), sondes
     / profils non convertis en 3D (U15), virgule décimale des sondes (U16).
  5. Relecture + 4 captures : vecteurs de bruit (|U_y|, |U_z| ≈ 7e-16) dessinés comme un
     écoulement secondaire (U17) ; textes périmés (T1, T2) ; doc 3D (D5 à D7).
- Pièges de la campagne : `set_combo` bloque les signaux Qt (deux faux constats évités en
  refaisant avec `setCurrentIndex`) ; une variante de ma part était fausse
  (`[physics] moment_center` : la clé est dans [output]) → devenue M16 (message sans la
  bonne section).
- Non testé : exécutables eux-mêmes (artefacts non téléchargeables), Windows, écran réel.
- Lots F1 à F4 placés avant l'import Gmsh 3D dans a_faire.md.

## 2026-10-03 — Lot #30 : écart transsonique NACA 0012 expliqué

- Demande : « Envoie la purée pour un lot entier ». Budget : une journée.
- Références : accès aux articles bloqué par le proxy (Stanford, scispace, arXiv, NASA
  GRC, su2code.github.io) ; non contourné. Valeurs trouvées seulement dans des résumés de
  moteur de recherche : Dolejší & Roskovec (DG hp-adaptatif, arXiv 2007.06840)
  C_l = 0.333, C_d = 0.02135 ; Vassberg & Jameson 2010 (profil pointu) « C_L ≈ 0.347 »,
  « C_D = 0.022453440 » ; workshop High-Order CFD : profil fermé −0.1036 (= le nôtre),
  chocs ≈ 0.6 et 0.35.
- Mesures (192 × 64 sauf mention, résidus < 1e-8) :
  - bord de fuite (K = 0.3) : closed 0.3337 / sharp 0.3343 / open 0.3460 ; 384 × 128 :
    0.3333 / 0.3340 / 0.3477. « ≈ 0.35 » = profil ouvert (AGARD).
  - limiteur : K = 0.05 → implicite 0.3139, RK3 0.3162 (96 × 32) : deux solutions ; K ≥ 0.1
    ou sans limiteur → une seule ; sans limiteur 0.3308 / 0.3332 / 0.3330 (96/192/384) ;
    K = 0.3 0.3308 / 0.3337 / 0.3333 ; dépassement C_p au choc 0.04 (K 0.05) → 0.23 (0.3)
    → 0.29 (sans) ; K = 0.1 : 0.10, C_l 0.3267 (96) / 0.3348 (192).
  - J'avais d'abord conclu à un effet géométrique (sharp −2.2 %) : faux, artefact de
    K = 0.05 (avec 0.3 : +0.2 %). Corrigé avant d'écrire la doc.
  - sans effet notable : 128 mailles radiales (−0.04 %), valeur de paroi des gradients
    extrapolée (+0.25 %, essai non retenu).
  - entropie parasite dans la 1re maille de paroi : 6.3e-3 (K 0.05), 3.8e-3 (sans limiteur),
    née au bord d'attaque ; amont < 3e-8. C_d +2 % (+0.0004) non expliqué (traînée parasite
    subsonique 384 × 128 : 0.00013, un tiers).
  - M = 0.5 : K = 0.05 abaisse C_l de 0.2-0.5 % et augmente la traînée parasite de 31-57 %.
  - méthode des panneaux (Hess-Smith, vérifiée Joukowski / Kármán-Trefftz) : C_l
    incompressible 0.15093 ; corrections de compressibilité trop dispersées (0.174-0.182 à
    M = 0.5) pour un contrôle fin.
- Fait : `trailing_edge = closed | open | sharp` (NACA, défaut inchangé, testé), exemple
  `venkat_k = 0.3`, test `test_transonic_example_limiter_threshold_matches_unlimited_solution`
  (22.6 s), figure régénérée, docs (compressible § 2, 3.5, 4, 5 ; README), E731 corrigée.
- Chronométré machine libre : exemple 192 × 64 50 s (452 it) contre 112 s (1 011 it) avec
  K = 0.05 ; 384 × 128 409 s (868 it).
- Commit 688df2a ; suite locale 454 réussis, 1 ignoré (6 min 37 s). CI tests verte
  (run 37145806547), exécutables verts (run 37145810874, Windows 11281804767, Linux
  11282345490).
- Environnement : `gh api` refusé pour ce dépôt (ancien nom → redirection par identifiant
  numérique refusée ; nouveau nom « non activé pour la session ») ; les outils MCP
  (owner Kayoouu, repo Claude-test) marchent : les utiliser pour suivre la CI.

## 2026-10-03 — Lot E4 : coupes x / y dans l'interface 3D

- Relance quotidienne (routine) ; CI verte au départ.
- Fait : `PlaneSlice` (`mesh3d/slice.py`) : cellules coupées repérées par le signe de
  x − c aux sommets, intersection des arêtes (tableaux par type : hexa 12 arêtes, prisme 9,
  pyramide 8, tétraèdre 6), points en double fusionnés (plan par un sommet), polygones
  d'aire nulle écartés, sommets rangés par angle autour du barycentre ; plan sur une face :
  cellule du dessus (comme ZSlice), face extrême haute : cellule du bord. `slice_mesh` :
  ZSlice pour z (inchangé), PlaneSlice pour x / y.
- Interface : « Plan de coupe (3D) » = axe (z / x / y) + cote ; vecteurs = composantes
  dans le plan ; vorticité = composante normale au plan (ω_x, ω_y, ω_z) ; axes nommés ;
  zoom sur les corps seulement en coupe z ; vue du maillage : plan dessiné en perspective.
  `_selftest` : coupe x de la conduite (titre « plan x = 0.25 », axe horizontal y).
- Vérifié : somme des aires = aire exacte de la section à 1e-15 près (pavé gradué, plans par
  les nœuds et les faces extrêmes ; extrusion autour d'un cylindre : épaisseur × longueur de
  la droite dans le maillage 2D ; 1 296 tétraèdres de Kuhn déformés ; 6 pyramides, plan par
  le sommet commun) ; coupe z PlaneSlice = ZSlice (mêmes cellules, mêmes aires). Coût 10⁶
  hexaèdres : 0.35 à 0.41 s (ZSlice 0.65 s).
- Captures regardées : conduite carrée x = cte (section, U_x max au centre), cavité cubique
  x = 0.5 (écoulement secondaire), y = 0.5 ω_y (antisymétrique par rapport à z = 0.5,
  attendu), vue du maillage avec le plan (plan gris discret mais visible).
- Défauts trouvés en route : numérotation des sommets cumulée par ligne au lieu de
  globalement (aire 4.13 au lieu de 3 : vue par le contrôle des aires) ; variable de boucle
  `w` qui écrasait la page Résultats (Qt : « QVBoxLayout already deleted », vue par les
  tests). Corrigés avant commit.
- Commit c1747e2 ; suite 452 réussis / 1 ignoré ; CI verte (tests 3.10 / 3.12 ; exécutables
  Windows et Linux par workflow_dispatch, run 37095760981, auto-test avec coupe x).

## 2026-10-02 — Lot E3 : VTK binaire

- Fait : `write_vtk(..., binary=True)` (VTK legacy BINARY, gros-boutiste, `>f8` points et
  champs, `>i4` connectivité ; en-têtes UTF-8 comme le texte), ancien texte gardé dans
  `_write_vtk_ascii` (identique octet pour octet à l'ancien code sur 5 maillages, polygones
  compris), `read_vtk` (texte et binaire), `[output] vtk_format = binary | ascii` (valeur
  inconnue refusée avant calcul), `write_mesh(binary=)`. Interface : clé conservée (le
  formulaire n'écrit que ses champs).
- Mesures côte à côte, mêmes champs (cavité cubique, U + 4 scalaires) : 10⁶ cellules
  14.0–15.1 s / 173 Mo → 0.31–0.35 s / 121 Mo ; 64³ 1.9 s / 48 Mo → 0.07 s / 32 Mo ; 32³
  5.6 → 4.0 Mo. Gain en taille limité (~30 %) : double précision gardée (exactitude).
- Vérifié avec le lecteur officiel VTK 9.7 (vtkUnstructuredGridReader, installé dans le
  conteneur, pas une dépendance) : triangles, quads, polygones (type 7), hexaèdres,
  prismes + hexa, vecteurs 2D complétés, noms accentués ; valeurs et points au bit près ;
  volumes VTK / solveur 1.7e-14 (2e-8 en texte).
- Défaut trouvé à la relecture du diff : nom de champ accentué (scalaire passif nommé par
  l'utilisateur) → UnicodeEncodeError en binaire (en-tête encodé ASCII) ; corrigé (UTF-8)
  et testé.
- Constat annexe : hexaèdres à faces gauches → volumes VTK ≠ volumes du solveur par cellule
  (proportionnel au gauchissement, total égal) ; nos maillages actuels ont des faces planes ;
  noté pour la tâche « import 3D ».
- Tests : `tests/test_vtk.py` (8), `test_case3d` / `test_mesh2d` relisent par `read_vtk` ;
  CI exécutables : en-tête « BINARY » vérifié sur le calcul 3D.
- Commit a4693bd ; CI verte (tests 3.10 / 3.12, exécutables Windows et Linux, run
  37076055233). Titre du commit « 40 fois » : vrai à 10⁶ cellules (40 à 49), 27 fois à 64³
  (corps du message exact ; historique non réécrit).

## 2026-10-02 — Lot E2 : distance à la paroi accélérée

- Profil (64³) : 91 % des cellules « incertaines » (test best ≤ d_k − R avec R = plus grand
  rayon de face), ~43 faces calculées chacune ; 28 s sur 37 dans le point le plus proche
  sur triangle. Pire : paroi à faces de tailles très différentes (cylindre extrudé sur un
  fond maillé) : 112 s pour 24 576 cellules, 207 s pour 64 179 (82 M paires).
- Fait : borne inférieure exacte par face (`_face_lower_bound` : écart au plan + écart dans
  le plan, marge 1e-9), arbre de boîtes englobantes commun 2D / 3D (`mesh2d/bvh.py`) à la
  place de query_ball_point, recherche k-d sur tous les cœurs (`workers=-1` : 18.8 s →
  8.3 s à 10⁶, même résultat). 2D : k-d + arbre de boîtes au lieu de tous les segments,
  égalités départagées comme avant (plus petit numéro).
- Essayé et rejeté : supprimer la recherche k-d (première estimation par descente vers la
  boîte la plus proche) → 64³ : 228 s (boîtes contenant le point, distance 0) et vecteurs
  différents aux égalités (plan diagonal de la cavité).
- Mesures côte à côte, même jour : distance 32³ 2.72 → 0.41 s, 64³ 34 → 3.7 s, 100³
  206 → 25.5 s ; calcul complet 10⁶ cellules 10 it : 375 → 195 s ; fields.vtk et
  history.csv identiques octet pour octet ; 2D 490 000 cellules 91 → 2.7 s.
- Tests : `tests/test_wall_distance.py` (force brute sur faces gauches, tailles de faces
  ×200, prismes + paroi courbe, 2D étirés, structure de l'arbre) ; ils passent tous par
  l'arbre de boîtes (vérifié). Avertissement « gros maillage » : préparation ~1 à 2 min à
  10⁶ (avant 3 à 5).

## 2026-10-02 — Mémoire persistante

- `CLAUDE.md` (consignes, préférences, pièges), `.claude/memoire/` (état, liste de travail,
  journal), hooks `.claude/settings.json` : SessionStart réinjecte `etat.md` et
  `a_faire.md` ; Stop bloque une fois si du code a été commité après la dernière mise à jour
  de la mémoire.
- `.gitignore` : `.claude/*` reste ignoré sauf `settings.json`, `memoire/`, `hooks/`.
- Autonomie : l'utilisateur a choisi « 1 lot par jour sauf exception » → routine
  `trig_01RCTjquBZw82GuKW83rPuSG`, tous les jours à 03h53 UTC dans cette session.
- Non vérifiable dans le tour : les hooks du projet ne sont chargés qu'au démarrage d'une
  session (ou relance du conteneur) ; scripts testés à la main (pipe-test), JSON validé.

## 2026-10-01 — Lot E : 3D dans l'exécutable (f08fb62, b289a91)

- Demande : « La 3D et le plus gros doit fonctionner coûte que coûte sur le .exe. Quitte à
  changer de format / de langage ». Changer de langage : écarté, mesure à l'appui
  (exe = Python à 1 % près sur 10⁶ cellules).
- Interface : type « Pavé 3D », cadre d'extrusion (case, z0, z1, nz, faces z périodiques /
  symétrie / parois / à régler), vecteurs à 3 composantes (`Vec2.set_dim`, `follow_dim`,
  `z_default` : raffinement z = 1, pas 0 — bug trouvé par le test : cas refusé « rapports
  > 0 »), colonne Uz des conditions limites, page Résultats : plan de coupe z, champ Uz,
  vue du maillage en perspective + coupe, plan hors domaine → message (ZSlice lève
  ValueError), messages du canevas repliés (textwrap).
- Tests : 3 tests d'interface 3D (formulaire pavé, case d'extrusion aller-retour, maillage
  + calcul + tracés) ; `_selftest` : conduite carrée 3D 128 cellules ; build.yml : calcul
  et maillage 3D par l'exe en CI (Windows et Linux verts).
- Mesures : 10⁶ cellules exe 461 s / 13.8 s par it / 3.09 Go ; Python 456 s / 13.6 s /
  3.18 Go ; la valeur publiée 8.5 s / it venait d'un jour où la machine était plus rapide →
  README, dépannage et avertissement donnent la fourchette 8.5 à 14 s. Interface 64³ :
  maillage 11.6 s, 3 it 65.5 s, figures < 1 s, 1.19 Go.
- Limite assumée : coupes z seulement (pour la conduite carrée, la section x = cte serait
  plus parlante) → tâche E4.

## 2026-10-01 — Lot D4 : 3D, périmètre réduit (1701163 → 2f818ee)

- D4.1 maillage 3D (Mesh3D hexa/prismes, pavé, extrusion) ; D4.2 opérateurs VF en
  dimension d ; D4.3 solveur incompressible 3D validé (conduite carrée : débit +0.38 % sur
  la série exacte) ; D4.4 turbulence + PIMPLE 3D (ABC, canal périodique identique au 2D et
  au 1D, conduite carrée turbulente SST : λ à 10 % de Blasius, pas d'écoulement secondaire —
  limite des modèles linéaires) ; D4.5 fichier de cas 3D, VTK 3D, sondes, reprise, efforts
  (A_ref = L_ref × étendue z, Cs, Cm autour de z), exemples, docs.
- Constats : extrusion multi-couches ≠ 2D (Rhie-Chow, Cd 2.17561 au lieu de 2.17575) et
  SIMPLE 2× plus lent ; une couche entre symétries = 2D exact. Résidu U_z plancher ~1e-9.
- Défauts de coût corrigés : distance à la paroi calculée deux fois (cache vidé par
  set_patch_types), géométrie des faces d'un bloc (5.2 Go → 1.6 Go par paquets), listes de
  query_ball_point (~5 Go) par paquets.

## 2026-10-01 — Audit utilisateur et lots 2 à 4 (1c07b07 → 283313e)

- Audit débutant / expert (`docs/audit_utilisateur.md`) : bugs critiques de l'interface,
  fautes de frappe signalées, valeurs impossibles refusées avant calcul, ligne de commande
  en français sans trace brute, interface (U1 à U13), tutoriel, dépannage, glossaire,
  référence des clés générée. Tous les points du tableau sont « corrigé ».

## 2026-09-30 — Lots D1 à D3 et extensions

- D1 noyaux Numba facultatifs (8ac9dbb, ~10 % de gain, multi-fil plus lent sur la VM) ;
  D2 solveur compressible en densité (Roe/HLLC, MUSCL, implicite SGS-GMRES, Sod, rampe,
  plaque, NACA transsonique : C_l ~4.5 % bas non expliqué) ; D3 solveur couplé
  pression-vitesse (1a79008).
- Scalaires passifs, non newtonien, zones poreuses, swirl, disques actuateurs, transition γ
  (Menter 2015), balayages parallèles, backends Intel / AMD, durcissement sécurité avant
  passage public, dépôt renommé Kayoouu/Microrans.

## 2026-09-29 — Socle

- 1D RANS/URANS (SA, k-ε, k-ω, SST), mailleur 2D, solveur 2D SIMPLE/PIMPLE, AMG, schémas en
  temps, interface PySide6, exécutables PyInstaller (macOS retiré ensuite), backend GPU
  (CuPy, non testé sur carte réelle), énergie + Boussinesq, lois de paroi, reprise,
  polaires, axisymétrique, FMG, sondes, débit / pression totale, animations, licence MIT.
