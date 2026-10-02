# Journal des lots (le plus récent en haut)

Détail des chiffres : README (§ 6 validation, § 7 performances, § 8 limites) et messages de
commit. Ce journal sert à retrouver ce qui a été fait, pourquoi, et ce qui a été constaté.

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
