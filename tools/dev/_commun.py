"""Fonctions communes aux outils de développement (tools/dev) : racine du dépôt, appels git,
instantanés (révision ou arbre de travail), exécution mesurée (temps, CPU, pic mémoire)."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RACINE = Path(__file__).resolve().parents[2]
BRANCHE = "claude/sweet-knuth-tdnyri"
DEPOT = "Kayoouu/Microrans"


def dossier_travail(nom: str) -> Path:
    """Dossier hors du dépôt pour les instantanés et les sorties (MICRORANS_DEV_TMP, sinon
    <tmp>/microrans_dev)."""
    base = Path(os.environ.get("MICRORANS_DEV_TMP",
                               Path(tempfile.gettempdir()) / "microrans_dev"))
    d = base / nom
    d.mkdir(parents=True, exist_ok=True)
    return d


def git(*args, cwd=RACINE, check=True, env=None) -> str:
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, env=env)
    if check and r.returncode:
        raise RuntimeError(f"git {' '.join(args)} : {r.stderr.strip()}")
    return r.stdout.strip()


def fichiers_arbre() -> list[str]:
    """Fichiers suivis et nouveaux non ignorés de l'arbre de travail (supprimés exclus)."""
    out = git("ls-files", "-co", "--exclude-standard", "-z")
    return [f for f in out.split("\0") if f and (RACINE / f).is_file()]


def instantane_revision(ref: str = "HEAD") -> Path:
    """Copie des fichiers suivis à la révision `ref` (git archive), en cache par commit."""
    sha = git("rev-parse", f"{ref}^{{commit}}")
    d = dossier_travail("revisions") / sha[:12]
    if not (d / ".complet").exists():
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
        arch = subprocess.Popen(["git", "archive", sha], cwd=RACINE, stdout=subprocess.PIPE)
        subprocess.run(["tar", "-x", "-C", str(d)], stdin=arch.stdout, check=True)
        if arch.wait():
            raise RuntimeError(f"git archive {sha} a échoué")
        (d / ".complet").touch()
    return d


def copier_arbre(dest: Path, fichiers=None) -> Path:
    """Copie de l'arbre de travail courant (ou des seuls `fichiers`) dans `dest`."""
    for f in fichiers if fichiers is not None else fichiers_arbre():
        cible = dest / f
        cible.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(RACINE / f, cible)
    return dest


def arbres_code(ref: str | None = None) -> dict:
    """Identifiants git des sous-arbres microrans/ et tests/ : de la révision `ref`, ou de
    l'arbre de travail (index temporaire, l'index réel n'est pas touché) si ref est None."""
    if ref is None:
        with tempfile.TemporaryDirectory() as t:
            env = dict(os.environ, GIT_INDEX_FILE=str(Path(t) / "index"))
            git("add", "-A", "--", "microrans", "tests", env=env)
            arbre = git("write-tree", env=env)
    else:
        arbre = git("rev-parse", f"{ref}^{{tree}}")
    return {d: git("rev-parse", f"{arbre}:{d}", check=False) for d in ("microrans", "tests")}


def env_microrans(racine: Path) -> dict:
    """Environnement pour lancer microrans depuis `racine` sans fenêtre ni affichage."""
    return dict(os.environ, PYTHONPATH=str(racine), MPLBACKEND="Agg",
                QT_QPA_PLATFORM="offscreen")


# Processus intermédiaire : RUSAGE_CHILDREN ne compte alors que la commande mesurée.
_MESURE = (
    "import json, resource, subprocess, sys, time\n"
    "t = time.perf_counter()\n"
    "r = subprocess.run(sys.argv[1:])\n"
    "u = resource.getrusage(resource.RUSAGE_CHILDREN)\n"
    "print('\\n@@mesure ' + json.dumps({'rc': r.returncode, "
    "'mur_s': time.perf_counter() - t, 'cpu_s': u.ru_utime + u.ru_stime, "
    "'pic_mo': u.ru_maxrss / 1024.0}))\n")


def executer_mesure(cmd, cwd, env=None, journal: Path | None = None) -> dict:
    """Lance `cmd` ; renvoie rc, temps mur (s), CPU (s) et pic mémoire (Mo) de la commande
    seule. Sorties de la commande écrites dans `journal` si donné."""
    p = subprocess.run([sys.executable, "-c", _MESURE, *map(str, cmd)], cwd=cwd, env=env,
                       capture_output=True, text=True)
    if journal is not None:
        journal.write_text(p.stdout + "\n--- stderr ---\n" + p.stderr, encoding="utf-8")
    for ligne in reversed(p.stdout.splitlines()):
        if ligne.startswith("@@mesure "):
            return json.loads(ligne[9:])
    raise RuntimeError(f"mesure impossible : {p.stderr.strip()[-500:]}")
