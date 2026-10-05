"""CI GitHub du dépôt (workflows « tests » et « executables »), par l'API (gh api).

    python tools/dev/ci.py etat              # état des workflows du dernier commit poussé
    python tools/dev/ci.py attendre          # attend la fin des workflows du commit poussé
    python tools/dev/ci.py attendre --executables   # lance aussi la construction (build.yml)

« attendre » : à lancer en arrière-plan ; échec → jobs et étapes en échec (le détail des
journaux : outil MCP get_job_logs). Exécutables : identifiants des artefacts et une ligne
prête pour etat.md. Codes de sortie : 0 tout vert, 1 échec, 2 délai dépassé ou API
inaccessible (alors : outils MCP actions_list / actions_get, ou add_repo Kayoouu/Microrans
avec access push).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _commun import BRANCHE, DEPOT, git  # noqa: E402


class ApiInaccessible(RuntimeError):
    pass


def api(chemin, methode="GET", champs=(), delai=20):
    cmd = ["gh", "api", "-X", methode, f"repos/{DEPOT}/{chemin}"]
    for c in champs:
        cmd += ["-f", c]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=delai)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise ApiInaccessible(str(e)) from e
    if r.returncode:
        raise ApiInaccessible((r.stderr or r.stdout).strip()[-300:])
    return json.loads(r.stdout) if r.stdout.strip() else {}


def executions(sha):
    return api(f"actions/runs?head_sha={sha}&per_page=30").get("workflow_runs", [])


def derniere_par_workflow(runs):
    """Exécution la plus récente de chaque workflow."""
    out = {}
    for r in sorted(runs, key=lambda r: r["created_at"]):
        out[r["name"]] = r
    return out


def ligne(r):
    etat = r["conclusion"] or r["status"]
    return f"{r['name']} : {etat} — {r['html_url']}"


def details_echec(r):
    jobs = api(f"actions/runs/{r['id']}/jobs?per_page=50").get("jobs", [])
    for j in jobs:
        if j.get("conclusion") not in ("success", "skipped", None):
            etapes = [s["name"] for s in j.get("steps", [])
                      if s.get("conclusion") not in ("success", "skipped", None)]
            print(f"    job « {j['name']} » (id {j['id']}) : {j['conclusion']} ; étapes en "
                  f"échec : {', '.join(etapes) or '?'}")


def artefacts(r):
    arts = api(f"actions/runs/{r['id']}/artifacts").get("artifacts", [])
    for a in arts:
        print(f"    artefact {a['name']} : id {a['id']}, expire le {a['expires_at'][:10]}")
    exe = {a["name"]: a for a in arts if a["name"].startswith("microrans-")}
    if exe:
        w, lx = exe.get("microrans-Windows"), exe.get("microrans-Linux")
        print(f"Pour etat.md : run {r['html_url']} (Windows : artefact "
              f"{w['id'] if w else '?'}, Linux : {lx['id'] if lx else '?'} ; expirent le "
              f"{(w or lx)['expires_at'][:10]})")


def sha_pousse():
    """Dernier commit de la branche distante (après git fetch)."""
    git("fetch", "-q", "origin", BRANCHE, check=False)
    return git("rev-parse", f"origin/{BRANCHE}")


def cmd_etat(court=False) -> int:
    sha = sha_pousse()
    runs = derniere_par_workflow(executions(sha))
    if court:
        etats = ", ".join(f"{n} {r['conclusion'] or r['status']}" for n, r in runs.items())
        print(f"CI du dernier commit poussé ({sha[:7]}) : {etats or 'aucune exécution'}")
        return 0
    print(f"Commit poussé : {sha[:7]} ({git('log', '-1', '--format=%s', sha)})")
    if git("rev-parse", "HEAD") != sha:
        print("  (HEAD local différent : commits non poussés)")
    for r in runs.values():
        print("  " + ligne(r))
    if not runs:
        print("  aucune exécution pour ce commit")
    return 0


def cmd_attendre(executables, delai) -> int:
    sha = sha_pousse()
    if git("rev-parse", "HEAD") != sha:
        print("HEAD n'est pas poussé : pousser d'abord (la CI suit le commit distant).")
        return 2
    # marge d'une minute : horloges locale et GitHub pas tout à fait synchrones
    debut = (datetime.now(timezone.utc) - timedelta(minutes=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
    if executables:
        api("actions/workflows/build.yml/dispatches", "POST", [f"ref={BRANCHE}"])
        print(f"Construction des exécutables lancée sur {sha[:7]}.")
    attendus = {"tests"} | ({"executables"} if executables else set())
    t0 = time.time()
    while True:
        runs = executions(sha)
        if executables:          # l'exécution lancée maintenant, pas une plus ancienne
            runs = [r for r in runs if r["name"] != "executables" or
                    r["created_at"] >= debut]
        der = derniere_par_workflow(runs)
        presents = attendus | set(der)
        finis = all(n in der and der[n]["status"] == "completed" for n in presents)
        if finis:
            break
        if time.time() - t0 > delai:
            print("Délai dépassé ; état :")
            for r in der.values():
                print("  " + ligne(r))
            return 2
        time.sleep(30)
    ok = True
    print(f"CI de {sha[:7]} terminée en {(time.time() - t0) / 60:.1f} min :")
    for r in der.values():
        print("  " + ligne(r))
        if r["conclusion"] != "success":
            ok = False
            details_echec(r)
        elif r["name"] == "executables":
            artefacts(r)
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("action", choices=("etat", "attendre"))
    ap.add_argument("--executables", action="store_true",
                    help="attendre : lancer aussi build.yml (workflow_dispatch)")
    ap.add_argument("--court", action="store_true", help="etat : une seule ligne")
    ap.add_argument("--delai", type=float, default=2400, help="attendre : secondes max")
    a = ap.parse_args()
    try:
        if a.action == "etat":
            return cmd_etat(a.court)
        return cmd_attendre(a.executables, a.delai)
    except ApiInaccessible as e:
        print(f"API GitHub inaccessible ({e}). Utiliser les outils MCP (actions_list, "
              f"actions_get), ou add_repo Kayoouu/Microrans avec access push.")
        return 2


if __name__ == "__main__":
    sys.exit(main())
