"""Suite complète des tests sur une COPIE de l'arbre de travail.

On peut continuer d'éditer pendant qu'elle tourne (plus de course avec
test_reference_document_up_to_date quand on régénère la doc). Le résultat est enregistré
avec l'identifiant des sous-arbres microrans/ et tests/ testés : le contrôle avant poussée
(avant_push.py) vérifie que le code poussé est exactement celui qui a passé la suite.

    python tools/dev/suite.py                # à lancer en arrière-plan (5 à 8 min)
    python tools/dev/suite.py tests/test_x.py -k reprise   # arguments passés à pytest

Sortie : journal complet dans <tmp>/microrans_dev/suite/suite.log, résumé en fin de sortie.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _commun import (arbres_code, copier_arbre, dossier_travail,  # noqa: E402
                     env_microrans, git)

RESULTATS = "resultats_suite.json"


def resultats_enregistres() -> list:
    f = dossier_travail("suite") / RESULTATS
    try:
        return json.loads(f.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []


def main(args) -> int:
    base = dossier_travail("suite")
    copie = base / "copie"
    shutil.rmtree(copie, ignore_errors=True)
    arbres = arbres_code()
    copier_arbre(copie)
    sha = git("rev-parse", "--short", "HEAD")
    modifie = bool(git("status", "--porcelain", "--untracked-files=no"))
    journal = base / "suite.log"
    cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *args]
    t0 = time.time()
    with open(journal, "w", encoding="utf-8") as fh:
        p = subprocess.Popen(cmd, cwd=copie, env=env_microrans(copie), stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True)
        for ligne in p.stdout:
            fh.write(ligne)
        rc = p.wait()
    duree = time.time() - t0
    lignes = journal.read_text(encoding="utf-8").splitlines()
    resume = next((x for x in reversed(lignes) if " in " in x and
                   ("passed" in x or "failed" in x or "error" in x)), "(pas de résumé)")
    echecs = [x for x in lignes if x.startswith(("FAILED", "ERROR"))]
    complete = not args
    res = {"date": time.strftime("%Y-%m-%d %H:%M:%S"), "commit": sha, "arbre_modifie": modifie,
           "arbres": arbres, "rc": rc, "complete": complete, "resume": resume.strip("= "),
           "duree_s": round(duree)}
    if complete:
        hist = [r for r in resultats_enregistres() if r.get("arbres") != arbres][-9:]
        (base / RESULTATS).write_text(json.dumps(hist + [res], indent=1), encoding="utf-8")
    print("\n".join(echecs[:30]))
    print(f"Suite {'complète' if complete else 'partielle'} sur {sha}"
          f"{' + modifications non commitées' if modifie else ''} : {res['resume']} "
          f"({duree / 60:.1f} min) ; journal : {journal}")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
