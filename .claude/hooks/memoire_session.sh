#!/bin/sh
# SessionStart (démarrage, reprise, après compaction) : réinjecte la mémoire de travail
# dans le contexte (la sortie standard d'un hook SessionStart est ajoutée au contexte).
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}" || exit 0
echo "=== Mémoire de travail du projet (.claude/memoire/, réinjectée automatiquement) ==="
for f in .claude/memoire/etat.md .claude/memoire/a_faire.md; do
  if [ -f "$f" ]; then
    echo
    echo "----- $f -----"
    cat "$f"
  fi
done
echo
echo "Journal détaillé des lots : .claude/memoire/journal.md (à lire si un détail manque)."
echo "Dernier commit : $(git log -1 --format='%h %ad %s' --date=short 2>/dev/null)"
# état de la CI du dernier commit poussé (API GitHub ; rouge → la corriger d'abord)
ci=$(timeout 10 python3 tools/dev/ci.py etat --court 2>/dev/null) || ci=""
echo "${ci:-CI : état inconnu (API inaccessible) ; vérifier avec les outils MCP.}"
echo "Outils : tools/dev (suite, echoue_avant, ab, ci, avant_push) ; procédure : skill lot."
exit 0
