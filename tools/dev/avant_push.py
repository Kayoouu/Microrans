"""Contrôles avant « git push » (règles de CLAUDE.md), lancés par le hook PreToolUse.

Refuse la poussée (code 2, raisons sur la sortie d'erreur) si :
- poussée forcée, poussée de tags, ou branche autre que la branche de travail ;
- ruff signale une erreur (microrans, tools, tests) ;
- une ligne ajoutée dans un fichier contient un nom de modèle, ou `import tomllib` dans
  tests/ (absent de Python 3.10, utilisé par la CI) ;
- les commits à pousser changent microrans/ ou tests/ sans que la suite complète
  (tools/dev/suite.py) ait passé sur ce code exact. Exception voulue : créer le fichier
  indiqué dans le message (trace de la décision).

    python tools/dev/avant_push.py            # contrôle manuel
    (hook) python tools/dev/avant_push.py --hook < entrée JSON de Claude Code
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _commun import BRANCHE, RACINE, arbres_code, dossier_travail, git  # noqa: E402

# noms de modèles écrits en morceaux : ce fichier ne doit pas lui-même les contenir
_FAMILLES = "|".join(["o" "pus", "son" "net", "hai" "ku", "fa" "ble"])
_MODELE_ID = re.compile(r"claude-(?:" + _FAMILLES + r")", re.IGNORECASE)
_MODELE_NOM = re.compile(r"\b(?:" + _FAMILLES.title() + r") [0-9]")
# « git push » en position de commande (début de ligne, après ; && || | ( do then else) :
# pas le mot push dans un message de commit, ni « git push » cité dans un texte
_PUSH = re.compile(r"(?:^[ \t]*|[;&|(][ \t]*|\b(?:do|then|else)[ \t]+)"
                   r"git[ \t]+(?:-[Cc][ \t]+\S+[ \t]+)*push\b([^;&|\n]*)", re.MULTILINE)


def problemes_commande(cmd: str) -> list[str]:
    """Règles sur la commande elle-même (sans regarder le dépôt)."""
    out = []
    for m in _PUSH.finditer(cmd):
        reste = m.group(1)
        mots = reste.split()
        if any(x in ("-f", "--force", "--force-with-lease", "--mirror") or
               x.startswith("--force") or x.startswith("+") for x in mots):
            out.append("poussée forcée interdite (CLAUDE.md : pas de force-push ; demander "
                       "à l'utilisateur)")
        if "--tags" in mots or "--follow-tags" in mots or any("refs/tags" in x for x in mots):
            out.append("poussée de tags : refusée par l'environnement (403) ; les releases "
                       "sont faites par l'utilisateur")
        if "--delete" in mots or "-d" in mots or any(x.startswith(":") for x in mots):
            out.append("suppression de branche distante interdite")
        if BRANCHE not in reste:
            out.append(f"écrire la branche explicitement : git push -u origin {BRANCHE}")
    return out


def lignes_interdites(diff: str) -> list[str]:
    """Lignes ajoutées interdites dans un diff unifié (git diff)."""
    out, fichier = [], ""
    for ligne in diff.splitlines():
        if ligne.startswith("+++ "):
            fichier = ligne[6:] if ligne.startswith("+++ b/") else ligne[4:]
            continue
        if not ligne.startswith("+") or ligne.startswith("+++"):
            continue
        txt = ligne[1:]
        if _MODELE_ID.search(txt) or _MODELE_NOM.search(txt):
            out.append(f"{fichier} : nom de modèle dans le dépôt (CLAUDE.md) : "
                       f"{txt.strip()[:100]}")
        if fichier.startswith("tests/") and re.search(r"^\s*(import|from)\s+tomllib\b", txt):
            out.append(f"{fichier} : `tomllib` absent de Python 3.10 (CI) : utiliser "
                       f"microrans.tomlio")
    return out


def _suite_ok(arbres: dict) -> bool:
    try:
        res = json.loads((dossier_travail("suite") / "resultats_suite.json")
                         .read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return any(r.get("arbres") == arbres and r.get("rc") == 0 and r.get("complete")
               for r in res)


def problemes_depot() -> list[str]:
    out = []
    if git("rev-parse", "--abbrev-ref", "HEAD") != BRANCHE:
        out.append(f"branche courante différente de {BRANCHE}")
    ruff = [shutil.which("ruff")] if shutil.which("ruff") else [sys.executable, "-m", "ruff"]
    r = subprocess.run([*ruff, "check", "microrans", "tools", "tests"], cwd=RACINE,
                       capture_output=True, text=True)
    if r.returncode:
        out.append("ruff :\n" + (r.stdout + r.stderr).strip()[-1500:])
    base = f"origin/{BRANCHE}"
    if not git("rev-parse", "--verify", "-q", base, check=False):
        return out
    out += lignes_interdites(git("diff", "--no-color", "-U0", f"{base}..HEAD"))
    changes = git("diff", "--name-only", f"{base}..HEAD").splitlines()
    if any(f.startswith(("microrans/", "tests/")) for f in changes):
        arbres = arbres_code("HEAD")
        if not _suite_ok(arbres):
            exc = dossier_travail("suite") / f"sans_suite_{arbres['microrans'][:12]}"
            if not exc.exists():
                out.append("microrans/ ou tests/ changés, mais la suite complète n'a pas passé "
                           "sur ce code exact : lancer `python tools/dev/suite.py` (arrière-"
                           "plan, 5 à 8 min) puis pousser. Exception voulue (la CI fera la "
                           f"suite) : `touch {exc}` et le dire dans le compte rendu.")
    return out


def main() -> int:
    if "--hook" in sys.argv:
        try:
            entree = json.load(sys.stdin)
        except ValueError:
            return 0
        cmd = (entree.get("tool_input") or {}).get("command", "")
        if not _PUSH.search(cmd):
            return 0
        problemes = problemes_commande(cmd) + problemes_depot()
    else:
        problemes = problemes_depot()
    if problemes:
        print("Poussée refusée par tools/dev/avant_push.py :\n- " + "\n- ".join(problemes),
              file=sys.stderr)
        return 2
    if "--hook" not in sys.argv:
        print("Contrôles avant poussée : OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
