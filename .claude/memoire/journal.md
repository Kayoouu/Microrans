# Journal des lots (le plus récent en haut)

Détail des chiffres : README (§ 6 validation, § 7 performances, § 8 limites) et messages de
commit. Ce journal sert à retrouver ce qui a été fait, pourquoi, et ce qui a été constaté.

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
