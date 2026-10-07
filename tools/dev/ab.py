"""Comparaisons avant / après : le même calcul microrans avec le code d'une révision (A, HEAD
par défaut) et avec l'arbre de travail (B).

Égalité des résultats (champs VTK, nombres de summary.json, CSV, checkpoint, au bit près) :
    python tools/dev/ab.py egalite -- run2d cavite_re100 -o {out} --no-plot -q \\
        --set solver.max_iter=30
Temps, CPU et pic mémoire, ordre A B B A (répété --repet fois ; machine libre) :
    python tools/dev/ab.py temps --repet 2 -- run2d plaque_plane_sa -o {out} --no-plot -q

{out} est remplacé par un dossier de sortie propre à chaque exécution. --ref change la
révision A ; --ref-b en donne une pour B au lieu de l'arbre de travail.
Code de sortie de « egalite » : 0 si tout est identique, 1 sinon.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import statistics
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(1, str(Path(__file__).resolve().parents[2]))
from _commun import (RACINE, dossier_travail, env_microrans, executer_mesure,  # noqa: E402
                     git, instantane_revision)

# clés de summary.json qui dépendent de la machine ou du dossier, pas du calcul
_IGNORE = ("wall_time_s", "wall_time", "time_mesh_s", "timings", "checkpoint", "backend",
           "out_dir", "output", "file")


def _aplatir(x, prefixe=""):
    if isinstance(x, dict):
        for k, v in x.items():
            if k in _IGNORE or k.endswith("_time_s"):
                continue
            yield from _aplatir(v, f"{prefixe}{k}.")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _aplatir(v, f"{prefixe}{i}.")
    else:
        yield prefixe.rstrip("."), x


def _ecart(a, b):
    try:
        a, b = np.asarray(a, float), np.asarray(b, float)
        if a.shape != b.shape:
            return f"tailles {a.shape} / {b.shape}"
        return f"écart max {float(np.max(np.abs(a - b))):.3g}"
    except (TypeError, ValueError):
        return "valeurs différentes"


def _comparer_json(fa, fb, out_a, out_b):
    da = dict(_aplatir(json.loads(fa.read_text(encoding="utf-8"))))
    db = dict(_aplatir(json.loads(fb.read_text(encoding="utf-8"))))
    diff = []
    for k in sorted(set(da) | set(db)):
        va, vb = da.get(k, "<absent>"), db.get(k, "<absent>")
        if isinstance(va, str) and isinstance(vb, str):
            va, vb = va.replace(str(out_a), "{out}"), vb.replace(str(out_b), "{out}")
        if va != vb:
            diff.append(f"{k} : {va!r} / {vb!r}")
    return diff


def _comparer_csv(fa, fb):
    with open(fa, newline="", encoding="utf-8") as ha, open(fb, newline="",
                                                             encoding="utf-8") as hb:
        la, lb = list(csv.reader(ha)), list(csv.reader(hb))
    if la == lb:
        return []
    if len(la) != len(lb) or not la or not lb:
        return [f"{len(la)} / {len(lb)} lignes, en-têtes {la[:1]} / {lb[:1]}"]
    diff = []
    if la[0] != lb[0]:              # colonnes ajoutées ou retirées : communes comparées
        plus = [c for c in lb[0] if c not in la[0]]
        moins = [c for c in la[0] if c not in lb[0]]
        diff.append("colonnes" + (f" ajoutées en B : {', '.join(plus)}" if plus else "")
                    + (f" retirées en B : {', '.join(moins)}" if moins else "")
                    + ("" if plus or moins else " dans un autre ordre"))
    for nom in [c for c in la[0] if c in lb[0]]:
        i, j = la[0].index(nom), lb[0].index(nom)
        try:
            ca = [float(r[i]) for r in la[1:]]
            cb = [float(r[j]) for r in lb[1:]]
        except (ValueError, IndexError):
            ca, cb = [r[i:i + 1] for r in la[1:]], [r[j:j + 1] for r in lb[1:]]
        if ca != cb:
            diff.append(f"colonne {nom} : {_ecart(ca, cb)}")
    if la[0] != lb[0] and len(diff) == 1:
        diff[0] += " ; colonnes communes identiques"
    return diff


def _comparer_npz(fa, fb):
    with np.load(fa, allow_pickle=False) as za, np.load(fb, allow_pickle=False) as zb:
        diff = [f"{k} absent d'un côté" for k in set(za.files) ^ set(zb.files)]
        for k in sorted(set(za.files) & set(zb.files)):
            a, b = za[k], zb[k]
            if k == "meta":
                ja, jb = (dict(_aplatir(json.loads(str(v)))) for v in (a, b))
                diff += [f"meta.{x}" for x in sorted(set(ja) | set(jb))
                         if ja.get(x) != jb.get(x)]
            elif not np.array_equal(a, b):
                diff.append(f"{k} : {_ecart(a, b)}")
    return diff


def _comparer_vtk(fa, fb):
    from microrans.mesh2d.io import read_vtk
    va, vb = read_vtk(fa), read_vtk(fb)
    diff = []
    for partie in ("points", "cells", "cell_types", "cell_data"):
        a, b = va.get(partie), vb.get(partie)
        if isinstance(a, dict) or isinstance(b, dict):
            a, b = a or {}, b or {}
            diff += [f"{partie}.{k} absent d'un côté" for k in set(a) ^ set(b)]
            diff += [f"{partie}.{k} : {_ecart(a[k], b[k])}" for k in sorted(set(a) & set(b))
                     if not np.array_equal(a[k], b[k])]
        elif a is not None and b is not None and not np.array_equal(a, b):
            diff.append(f"{partie} : {_ecart(a, b)}")
    return diff


def comparer_dossiers(a: Path, b: Path) -> dict:
    """{fichier relatif : liste des différences} ; liste vide = identique ; None = format
    non comparé (images…)."""
    fa = {p.relative_to(a).as_posix() for p in a.rglob("*") if p.is_file()}
    fb = {p.relative_to(b).as_posix() for p in b.rglob("*") if p.is_file()}
    res = {}
    for rel in sorted(fa | fb):
        if rel not in fa or rel not in fb:
            res[rel] = [f"présent seulement dans {'B' if rel in fb else 'A'}"]
            continue
        x, y = a / rel, b / rel
        suf = x.suffix.lower()
        if suf == ".json":
            res[rel] = _comparer_json(x, y, a, b)
        elif suf in (".csv", ".dat"):
            res[rel] = _comparer_csv(x, y)
        elif suf == ".npz":
            res[rel] = _comparer_npz(x, y)
        elif suf == ".vtk":
            res[rel] = _comparer_vtk(x, y)
        else:
            res[rel] = None
    return res


def _cote(nom, racine, args, n):
    out = dossier_travail("ab") / f"{nom}{n}"
    shutil.rmtree(out, ignore_errors=True)
    cmd = [sys.executable, "-m", "microrans", *[x.replace("{out}", str(out)) for x in args]]
    m = executer_mesure(cmd, racine, env_microrans(racine),
                        journal=dossier_travail("ab") / f"{nom}{n}.log")
    return out, m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("egalite", "temps"))
    ap.add_argument("--ref", default="HEAD", help="révision A (défaut HEAD)")
    ap.add_argument("--ref-b", default=None, help="révision B (défaut : arbre de travail)")
    ap.add_argument("--repet", type=int, default=1, help="temps : nombre de cycles A B B A")
    argv = sys.argv[1:]
    if "--" not in argv:
        ap.error("donner la commande microrans après -- (voir l'aide)")
    k = argv.index("--")
    a, args = ap.parse_args(argv[:k]), argv[k + 1:]
    if not args:
        ap.error("donner la commande microrans après --")
    ra = instantane_revision(a.ref)
    rb = instantane_revision(a.ref_b) if a.ref_b else RACINE
    nom_a = git("rev-parse", "--short", a.ref)
    nom_b = git("rev-parse", "--short", a.ref_b) if a.ref_b else "arbre de travail"
    print(f"A = {nom_a}, B = {nom_b}")
    if a.mode == "egalite":
        oa, ma = _cote("A", ra, args, 0)
        ob, mb = _cote("B", rb, args, 0)
        print(f"codes de retour : A {ma['rc']}, B {mb['rc']} ; durées {ma['mur_s']:.1f} s / "
              f"{mb['mur_s']:.1f} s")
        for cote, m in (("A", ma), ("B", mb)):
            if m["rc"] not in (0, 1):            # 1 = non convergé, résultats écrits
                log = (dossier_travail("ab") / f"{cote}0.log").read_text(encoding="utf-8")
                print(f"ERREUR côté {cote} (code {m['rc']}) :\n{log.strip()[-800:]}")
                return 2
        res = comparer_dossiers(oa, ob)
        if not res:
            print("Aucun fichier de sortie : rien à comparer (vérifier -o {out}).")
            return 2
        egal = ma["rc"] == mb["rc"]
        for rel, d in res.items():
            if d is None:
                print(f"  {rel} : non comparé")
            elif d:
                egal = False
                print(f"  {rel} : DIFFÉRENT")
                for x in d[:12]:
                    print(f"      {x}")
            else:
                print(f"  {rel} : identique")
        print("RÉSULTAT : identique au bit près" if egal else "RÉSULTAT : différences")
        return 0 if egal else 1
    charge = os.getloadavg()[0]
    if charge > 0.5:
        print(f"ATTENTION : charge moyenne {charge:.2f} (machine pas libre : mesures faussées)")
    mes = {"A": [], "B": []}
    n = 0
    for _ in range(a.repet):
        for cote in "ABBA":
            _, m = _cote(cote, ra if cote == "A" else rb, args, n)
            n += 1
            mes[cote].append(m)
            print(f"  {cote} : {m['mur_s']:.2f} s mur, {m['cpu_s']:.2f} s CPU, "
                  f"{m['pic_mo']:.0f} Mo, rc {m['rc']}", flush=True)
    med = {c: {k: statistics.median(m[k] for m in v) for k in ("mur_s", "cpu_s", "pic_mo")}
           for c, v in mes.items()}
    for k, u in (("mur_s", "s"), ("cpu_s", "s CPU"), ("pic_mo", "Mo")):
        print(f"médiane {k} : A {med['A'][k]:.2f} {u}, B {med['B'][k]:.2f} {u} "
              f"(B / A = {med['B'][k] / max(med['A'][k], 1e-12):.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
