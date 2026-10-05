"""Le test prouve-t-il la correction ? Lance des tests sur le code d'une révision (HEAD par
défaut, donc sans les modifications en cours) avec les tests ACTUELS, puis sur l'arbre de
travail. Attendu : échec avant, succès après. L'arbre de travail n'est pas touché
(remplace `git stash push … ; pytest ; git stash pop`).

    python tools/dev/echoue_avant.py tests/test_x.py::test_a tests/test_y.py::test_b
    python tools/dev/echoue_avant.py --ref HEAD~1 tests/test_x.py::test_a

Un échec par ImportError / AttributeError / NameError (fonction absente de l'ancien code)
est signalé « preuve faible » : il montre que l'API manque, pas que le défaut existait.
Code de sortie 0 seulement si tous les tests échouent avant (preuve forte) et passent après.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _commun import (RACINE, copier_arbre, dossier_travail, env_microrans,  # noqa: E402
                     fichiers_arbre, git, instantane_revision)

_FAIBLE = ("ImportError", "ModuleNotFoundError", "AttributeError", "NameError", "KeyError",
           "unexpected keyword argument", "fixture '")


def lancer(racine: Path, noeuds, junit: Path) -> dict:
    """{nom du test : (état, message)} ; état = passé | échec | erreur | ignoré."""
    subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                    f"--junitxml={junit}", *noeuds], cwd=racine, env=env_microrans(racine),
                   capture_output=True, text=True)
    res = {}
    if not junit.is_file():
        return res
    for tc in ET.parse(junit).getroot().iter("testcase"):
        nom = tc.get("name")
        etat, msg = "passé", ""
        for balise, e in (("failure", "échec"), ("error", "erreur"), ("skipped", "ignoré")):
            el = tc.find(balise)
            if el is not None:
                etat = e
                msg = (el.get("message") or el.text or "").strip()
                break
        res[nom] = (etat, msg)
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref", default="HEAD", help="révision « avant » (défaut HEAD)")
    ap.add_argument("noeuds", nargs="+", help="tests pytest (fichier::test)")
    a = ap.parse_args()
    base = dossier_travail("echoue_avant")
    avant = base / "avant"
    shutil.rmtree(avant, ignore_errors=True)
    shutil.copytree(instantane_revision(a.ref), avant)
    copier_arbre(avant, [f for f in fichiers_arbre() if f.startswith("tests/")])
    r_avant = lancer(avant, a.noeuds, base / "avant.xml")
    r_apres = lancer(RACINE, a.noeuds, base / "apres.xml")
    if not r_avant and not r_apres:
        print("Aucun test exécuté : vérifier les noms (fichier::test).")
        return 2
    ok = True
    ref = git("rev-parse", "--short", a.ref)
    for nom in sorted(set(r_avant) | set(r_apres)):
        ea, ma = r_avant.get(nom, ("absent", ""))
        ep, mp = r_apres.get(nom, ("absent", ""))
        if ea in ("échec", "erreur"):
            faible = any(s in ma for s in _FAIBLE)
            verdict = "preuve faible (API absente avant)" if faible else "échoue avant : OK"
        else:
            faible, verdict = False, f"NE prouve PAS la correction ({ea} sur {ref})"
        ok &= ea in ("échec", "erreur") and not faible and ep == "passé"
        print(f"{nom} : avant ({ref}) {ea} ; après {ep} → {verdict}")
        if ea in ("échec", "erreur"):
            print(f"    avant : {ma.splitlines()[0][:160] if ma else ''}")
        if ep != "passé":
            print(f"    après : {mp.splitlines()[0][:160] if mp else ''}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
