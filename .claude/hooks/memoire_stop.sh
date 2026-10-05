#!/bin/sh
# Stop : bloque la fin du tour une fois (stop_hook_active évite la boucle) si
# - du code a été commité après la dernière mise à jour de .claude/memoire/ ;
# - du travail risque d'être perdu (conteneur éphémère) : fichiers suivis modifiés non
#   commités, ou commits non poussés.
input=$(cat)
case "$input" in
  *'"stop_hook_active":true'* | *'"stop_hook_active": true'*) exit 0 ;;
esac
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}" || exit 0
raisons=""
code=$(git log -1 --format=%ct -- microrans tests tools packaging .github README.md docs \
       2>/dev/null)
mem=$(git log -1 --format=%ct -- .claude/memoire 2>/dev/null)
if [ -n "$code" ] && { [ -z "$mem" ] || [ "$code" -gt "$mem" ]; }; then
  raisons="Mémoire à mettre à jour : du code a été commité après la dernière modification de .claude/memoire/ (etat.md, a_faire.md, journal.md : résultats, mesures, décisions, prochaine étape), puis commiter et pousser."
fi
modifs=$(git status --porcelain 2>/dev/null | head -5 | tr '\n' ' ')
nonpousses=$(git log --oneline "@{u}..HEAD" 2>/dev/null | wc -l)
if [ -n "$modifs" ] || [ "${nonpousses:-0}" -gt 0 ]; then
  raisons="$raisons Travail non sauvegardé (conteneur éphémère) : fichiers modifiés ou nouveaux non commités [$modifs] ; $nonpousses commit(s) non poussé(s). Commiter et pousser, ou, si c'est voulu (travail en cours, question posée), le dire en une ligne dans la réponse."
fi
if [ -n "$raisons" ]; then
  python3 -c 'import json, sys; print(json.dumps({"decision": "block", "reason": sys.argv[1].strip()}, ensure_ascii=False))' "$raisons"
fi
exit 0
