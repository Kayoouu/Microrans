# Journal des lots (le plus récent en haut)

Détail des chiffres : README (§ 6 validation, § 7 performances, § 8 limites) et messages de
commit. Ce journal sert à retrouver ce qui a été fait, pourquoi, et ce qui a été constaté.

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
