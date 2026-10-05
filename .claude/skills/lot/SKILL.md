---
name: lot
description: Procédure d'un lot de travail microrans (relance quotidienne, ou lot demandé par l'utilisateur) - départ, preuve que chaque test échoue sans la correction, comparaisons avant / après, suite complète, poussée, CI, mémoire, rétrospective. À suivre pour tout lot de a_faire.md.
---

# Procédure d'un lot

Un lot = une tâche de `.claude/memoire/a_faire.md`, dans l'ordre. Les règles de fond sont
dans `CLAUDE.md` ; cette page dit comment les appliquer avec les outils `tools/dev`.

## 1. Départ (5 min)

1. `git pull origin claude/sweet-knuth-tdnyri`. L'état de la CI du dernier commit poussé est
   affiché au démarrage (hook) ; sinon `python tools/dev/ci.py etat`. **CI rouge : la
   corriger d'abord**, c'est le lot du jour.
2. Lire la tâche dans `a_faire.md` et, pour un point d'audit, sa ligne dans
   `docs/audit_utilisateur.md` (reproduction).
3. Lundi (ou si la dernière campagne de non-régression notée dans `bilan_lots.md` date de
   plus de 7 jours) : relancer d'abord les campagnes c1, c6, c7, c8 de `tools/audit`
   (en arrière-plan) et comparer au dernier relevé.
4. Écrire le lot dans « En cours » de `a_faire.md` avec : points, critère « fait quand »,
   estimation (lots ou heures). Une reprise après compaction repart de là.

## 2. Chaque point

1. **Reproduire d'abord** (commande ou test qui montre le défaut), puis corriger.
2. **Le test prouve la correction** :
   `python tools/dev/echoue_avant.py tests/test_x.py::test_y`. Attendu : « échoue avant :
   OK ». « Preuve faible » (échec par fonction ou clé absente) : ajouter une assertion sur
   le comportement, qui échoue sur l'ancien code.
3. **Résultats inchangés là où ils doivent l'être** : `python tools/dev/ab.py egalite --
   run2d <exemple> -o {out} --no-plot -q --set solver.max_iter=30` sur les exemples
   touchés ; une différence doit être voulue et expliquée.
4. **Tout chiffre écrit est mesuré** : temps et mémoire par
   `python tools/dev/ab.py temps --repet 2 -- …` (ordre A B B A, machine libre : la charge
   est affichée) ; préciser les conditions dans le texte.
5. Relecture adverse du diff (`git diff`) avec la liste des défauts déjà trouvés par les
   audits : chemins 2D / 3D / axisymétrique, reprise (exacte ?), interface = ligne de
   commande, dossier courant ≠ dossier du cas, Windows (chemins, encodage, virgule
   décimale), Python 3.10, exécutable (imports cachés de PyInstaller), messages en
   français, documentation / référence des clés à jour.

## 3. Avant de pousser

1. `python tools/dev/suite.py` **en arrière-plan** (suite complète sur une copie : on peut
   éditer la doc pendant ce temps). Le hook avant poussée la **réclame** si microrans/ ou
   tests/ changent ; il vérifie aussi ruff, branche, force, tags, noms de modèles, tomllib.
2. Docs (README, référence des clés : `python -m microrans.fv2d.validate >
   docs/reference_cas.md` si une clé change), statuts de l'audit.
3. Commit en français : quoi, pourquoi, chiffres mesurés ; lignes d'attribution.
4. `git push -u origin claude/sweet-knuth-tdnyri`.

## 4. Après la poussée

1. `python tools/dev/ci.py attendre --executables` en arrière-plan (`--executables` si
   microrans/ ou packaging/ ont changé) : conclusion, étapes en échec, artefacts et ligne
   prête pour `etat.md`. Échec : `get_job_logs` (MCP) pour le détail, corriger, repousser.
2. Campagnes `tools/audit` concernées par le lot, relancées.

## 5. Fin du lot

1. Mémoire : `etat.md` (dernier lot, liens CI / exécutables), `a_faire.md` (« En cours »
   vidé, tâche retirée), `journal.md` (détail, chiffres), `bilan_lots.md` (ligne : estimé,
   réel = relance → dernier commit, tests ajoutés, rétro).
2. **Rétrospective (2 min)** : qu'est-ce qui a coûté du temps ou causé une erreur ? Si un
   outil, un hook ou une ligne de cette procédure l'évite : le faire tout de suite si
   c'est moins de 15 min, sinon l'ajouter aux « Petits travaux ». Noter le changement et,
   aux lots suivants, son effet constaté (sinon le retirer).
3. Compte rendu en français : fait / ne marche pas / chiffres / liens (CI, exécutables) /
   prochaine tâche. Lot suivant seulement si petit (< 1 h).
