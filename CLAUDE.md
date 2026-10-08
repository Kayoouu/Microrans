# Consignes pour Claude (microrans)

Lu automatiquement au début de chaque session. La mémoire de travail (état, liste de
travail, journal) est dans `.claude/memoire/` : `etat.md` et `a_faire.md` sont réinjectés
automatiquement au démarrage et après chaque compaction (hook SessionStart) ; `journal.md`
est à lire quand il faut retrouver le détail d'un lot passé.

## Préférences de l'utilisateur (textuelles)

- « N'essaie jamais de me plaire dans tes réponses, je veux savoir exactement ce qui est bon
  et mauvais, ce qui est faisable ou non et ce qui est VRAI ou FAUX »
- « tu vas un peu trop vite, doucement » : un lot à la fois, montrer les résultats
  (chiffres mesurés, figures regardées) avant de passer au suivant.
- Répondre en français. Mettre les liens importants (exécutables, exécutions de CI) dans la
  réponse finale.
- « Dans l'idéal je ne devrais même pas intervenir et tu pourrais poursuivre en autonomie » :
  suivre `a_faire.md` sans redemander, dans les limites de la section « Autonomie » ci-dessous.

## Règles de travail

- Dépôt : github.com/Kayoouu/Microrans (outils GitHub MCP : owner `Kayoouu`, repo
  `Claude-test`, ancien nom, toujours accepté). Branche unique : `claude/sweet-knuth-tdnyri`
  (branche par défaut) ; travailler, commiter et pousser seulement là. Pas de pull request
  sauf demande.
- Messages de commit en français, détaillés (quoi, pourquoi, chiffres mesurés), terminés
  par les lignes d'attribution fournies par la session. Aucun nom de modèle dans le dépôt.
- Tags et releases : l'envoi d'un tag est refusé (403) depuis cet environnement ; ne pas
  contourner. La publication d'une release est faite par l'utilisateur.
- Tout chiffre écrit dans le README, la doc ou un message d'avertissement doit avoir été
  mesuré ; préciser les conditions (machine virtuelle 4 cœurs, durées variables d'un jour à
  l'autre, jusqu'à 1.7 fois pour le même code — cylindre RK3 : 1 061 s puis 1 780 s — et pas
  dans le même sens pour tous les cas : mesurer les comparaisons côte à côte, machine libre ;
  une durée publiée se donne en fourchette de mesures datées).
- Avant de pousser : ruff, tests concernés, suite complète si microrans/ ou tests/
  changent (le hook avant poussée le vérifie). Après la poussée : vérifier la CI (workflow
  `tests`, et `executables` si l'exécutable est concerné).

## Procédure et outils (détail : skill `lot`, `.claude/skills/lot/SKILL.md`)

- Chaque lot suit le skill `lot` : départ (CI, tâche), un test qui échoue sans la
  correction, comparaisons avant / après, suite, poussée, CI, mémoire, rétrospective.
- `tools/dev/` (scripts, `-h` pour l'aide) :
  - `suite.py` : suite complète sur une copie de l'arbre (arrière-plan ; on peut éditer
    pendant ce temps) ;
  - `echoue_avant.py` : le test échoue-t-il sur l'ancien code ? (sans `git stash`) ;
  - `ab.py egalite | temps` : mêmes sorties au bit près, ou temps / CPU / mémoire
    A B B A, entre une révision et l'arbre de travail ;
  - `ci.py etat | attendre [--executables]` : CI par l'API (gh api), artefacts ;
  - `avant_push.py` : contrôles avant poussée.
- Hooks (`.claude/settings.json`) :
  - démarrage : mémoire + état de la CI du dernier commit poussé ;
  - avant une poussée (commande Bash) : branche, pas de force ni de tag, ruff, pas de nom
    de modèle ni de `tomllib` dans les tests, suite passée sur ce code exact ;
  - fin de tour : mémoire à jour, rien de non commité ni de non poussé.
- `.claude/memoire/bilan_lots.md` : estimé / réel par lot, rétrospectives, indicateur du
  jalon A (résultats faux silencieux par audit).

## Autonomie

Rythme choisi par l'utilisateur : « 1 lot par jour sauf exception ». Une relance
quotidienne programmée (routine `trig_01RCTjquBZw82GuKW83rPuSG`, 03h53 UTC) réveille la
session ; chaque relance = un lot. Exceptions : lot fini tôt et suivant petit (< 1 h) →
l'enchaîner ; CI rouge → la corriger d'abord ; décision qui revient à l'utilisateur →
s'arrêter et demander. Si l'utilisateur demande d'arrêter : désactiver la routine
(`update_trigger`, enabled = false) et le noter dans `etat.md`.

Sans demander : les tâches de `a_faire.md` dans l'ordre, corrections de bugs, tests,
documentation, mesures, relance de la construction des exécutables (workflow_dispatch).
Réordonner la liste est permis si une mesure le justifie (noter la raison dans le journal).

Demander d'abord (et s'arrêter sur ce point) : changer de langage ou d'architecture,
supprimer une fonction existante, ajouter une dépendance lourde aux exécutables, toute
action sur les releases / tags / historique Git (pas de force-push), toute dépense ou
service extérieur, tout ce que la liste ne couvre pas et qui change le périmètre du projet.

Fin d'un lot : tests verts en local et en CI, `etat.md` / `a_faire.md` / `journal.md` mis à
jour, commit et poussée, compte rendu honnête (ce qui marche, ce qui ne marche pas,
chiffres), puis lot suivant.

## Mémoire : quand la mettre à jour

- À la fin de chaque lot, et dès qu'une mesure, une décision ou une limite change.
- Avant une réponse finale à l'utilisateur.
- Le hook Stop bloque une fois la fin du tour si du code a été commité après la dernière
  modification de `.claude/memoire/` : mettre à jour, commiter, pousser.
- `etat.md` et `a_faire.md` restent courts (ils sont réinjectés à chaque compaction) ; le
  détail va dans `journal.md`.

## Pièges connus de l'environnement

- Pas de `pytest-xdist` (`-n 4` refusé) : suite complète en série, 5 à 8 min
  (`python tools/dev/suite.py`, en arrière-plan).
- Scripts hors du dépôt : `PYTHONPATH=/home/user/Claude-test` ; interface hors écran :
  `QT_QPA_PLATFORM=offscreen` (et `MICRORANS_RESULTS=<dossier>` pour les sorties).
- Pas de `/usr/bin/time` : temps, CPU et pic mémoire par `tools/dev/ab.py temps`.
- `pkill -f motif` tue aussi le shell qui le lance (code 144) : viser un PID.
- `sleep` au premier plan est bloqué : boucles `for i in $(seq ..); do sleep 20; done` dans
  une commande, ou tâche en arrière-plan.
- Le téléchargement des artefacts de CI (blob.core.windows.net) est refusé par le proxy :
  vérifier les exécutables par les codes de retour des étapes de CI.
- Le workflow `executables` ne part tout seul que si `packaging/` ou
  `.github/workflows/build.yml` change ; sinon `python tools/dev/ci.py attendre
  --executables` (ou outil MCP `actions_run_trigger`, workflow `build.yml`).
- API GitHub par `gh api` : seulement avec le nom `Kayoouu/Microrans` (dépôt ajouté à la
  session par add_repo, access push, le 2026-10-05) ; l'ancien nom `Claude-test` redirige
  vers un chemin numérique refusé par le proxy. API refusée : outils MCP (owner
  `Kayoouu`, repo `Claude-test`).
- Le hook avant poussée lit les commandes Bash : une commande qui commence par
  « git push » dans un texte (heredoc) est prise pour une poussée ; éditer alors avec
  l'outil Edit / Write.
- La CI (`tests`) lance `python -m microrans verify` puis `pytest` sous Python 3.10 et 3.12,
  mais pas `ruff` : le lancer en local. `ruff check tests` signale une erreur E731
  préexistante (`tests/test_mesh2d.py:74`). Pas de `import tomllib` dans les tests
  (absent de Python 3.10 : utiliser `microrans.tomlio`).
- Le conteneur est éphémère : tout ce qui doit survivre va dans le dépôt (commit + push).
