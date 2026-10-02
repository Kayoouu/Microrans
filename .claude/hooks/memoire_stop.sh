#!/bin/sh
# Stop : si du code a été commité après la dernière mise à jour de .claude/memoire/,
# bloque la fin du tour une fois (stop_hook_active évite la boucle) pour faire mettre la
# mémoire à jour.
input=$(cat)
case "$input" in
  *'"stop_hook_active":true'* | *'"stop_hook_active": true'*) exit 0 ;;
esac
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}" || exit 0
code=$(git log -1 --format=%ct -- microrans tests packaging .github README.md docs \
       2>/dev/null)
mem=$(git log -1 --format=%ct -- .claude/memoire 2>/dev/null)
[ -n "$code" ] || exit 0
if [ -z "$mem" ] || [ "$code" -gt "$mem" ]; then
  printf '%s\n' '{"decision": "block", "reason": "Mémoire à mettre à jour : du code a été commité après la dernière modification de .claude/memoire/. Mettre à jour etat.md, a_faire.md et journal.md (résultats, mesures, décisions, prochaine étape), commiter, pousser, puis terminer."}'
fi
exit 0
