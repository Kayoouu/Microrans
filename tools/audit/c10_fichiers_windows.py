"""Campagne 10 : fichiers et chemins tels qu'on les trouve sur un poste Windows : cas
enregistré par le Bloc-notes (BOM UTF-8), fins de ligne CRLF, encodage Latin-1 (cp1252),
dossiers avec espaces et accents, sortie en lecture seule ou occupée par un fichier, contour
CSV « à la française » (point-virgule, virgule décimale), maillage Gmsh 4.1 ASCII et
binaire, cas JSON. Pour chaque essai : code de retour et dernière ligne utile."""

import json
import os
import resource
import shutil
import stat
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]).resolve()
if OUT.exists():
    shutil.rmtree(OUT, onerror=lambda f, p, e: (os.chmod(p, stat.S_IWRITE), f(p)))
OUT.mkdir(parents=True)
EX = REPO / "microrans/examples"
env = dict(os.environ, PYTHONPATH=str(REPO), MPLBACKEND="Agg")
CAV = (EX / "cavite_re100.toml").read_text(encoding="utf-8")
SHORT = ["--set", "solver.max_iter=5", "--no-plot"]
res = []


def run(tag, args, cwd=OUT, expect="ok"):
    def cap():                     # 4 Go : un mailleur emballé ne doit pas tuer la machine
        resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))
    p = subprocess.run([sys.executable, "-m", "microrans", *args], cwd=cwd, env=env,
                       capture_output=True, text=True, preexec_fn=cap, timeout=600)
    lines = [ln for ln in (p.stderr + "\n" + p.stdout).splitlines() if ln.strip()]
    msg = next((ln for ln in lines if ln.startswith(("Erreur", "ATTENTION"))),
               lines[-1] if lines else "")
    r = {"essai": tag, "rc": p.returncode, "attendu": expect, "message": msg[:300],
         "trace": "Traceback" in p.stderr}
    res.append(r)
    print(f"{tag:34s} rc={p.returncode} [{expect}] {msg[:160]}", flush=True)
    return p


def write(path, text, enc="utf-8", bom=False, crlf=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if crlf:
        text = text.replace("\n", "\r\n")
    data = text.encode(enc, errors="replace")
    path.write_bytes((b"\xef\xbb\xbf" if bom else b"") + data)
    return path


# 1. fichiers de cas
write(OUT / "bom.toml", CAV, bom=True)
run("BOM UTF-8 (Bloc-notes)", ["run2d", "bom.toml", "-o", "o_bom", *SHORT])
write(OUT / "crlf.toml", CAV, crlf=True)
run("fins de ligne CRLF", ["run2d", "crlf.toml", "-o", "o_crlf", *SHORT])
write(OUT / "latin1.toml", "# Cavité entraînée, réglée à Re = 100\n" + CAV, enc="cp1252")
run("encodage Latin-1 (cp1252)", ["run2d", "latin1.toml", "-o", "o_latin1", *SHORT],
    expect="erreur claire ou ok")
write(OUT / "bom_crlf.toml", CAV, bom=True, crlf=True)
run("BOM + CRLF", ["run2d", "bom_crlf.toml", "-o", "o_bomcrlf", *SHORT])
# 2. chemins
d = OUT / "Mes Documents" / "Été 2026"
write(d / "cas cavité.toml", CAV)
run("cas : espaces et accents", ["run2d", str(d / "cas cavité.toml"), "-o",
                                 str(d / "résultats cavité"), *SHORT])
ro = OUT / "lecture_seule"
ro.mkdir()
os.chmod(ro, stat.S_IREAD | stat.S_IEXEC)
run("sortie en lecture seule", ["run2d", "cavite_re100", "-o", str(ro / "res"), *SHORT],
    expect="erreur claire")
os.chmod(ro, stat.S_IWRITE | stat.S_IREAD | stat.S_IEXEC)
(OUT / "un_fichier").write_text("x")
run("sortie = fichier existant", ["run2d", "cavite_re100", "-o", str(OUT / "un_fichier"),
                                  *SHORT], expect="erreur claire")
run("examples : copie", ["examples", "cavite_re100", "-o", "copies"])
run("examples : copie déjà faite", ["examples", "cavite_re100", "-o", "copies"],
    expect="erreur claire ou refus")
# 3. contour CSV à la française (profil du volet de l'exemple multi-profils)
pts = [ln.split() for ln in (EX / "profil_volet.dat").read_text().splitlines()[1:] if ln.strip()]
fr = "x;y\n" + "\n".join(f"{x.replace('.', ',')};{y.replace('.', ',')}" for x, y in pts) + "\n"
write(OUT / "volet_fr.csv", fr, crlf=True)
us = "x,y\n" + "\n".join(f"{x},{y}" for x, y in pts) + "\n"
write(OUT / "volet_us.csv", us)
for tag, f in (("contour CSV x,y (point)", "volet_us.csv"),
               ("contour CSV x;y virgule déc.", "volet_fr.csv")):
    case = {"mesh": {"type": "unstructured", "h_max": 0.5, "h_surface": 0.02},
            "domain": {"type": "rectangle", "x0": -1, "x1": 3, "y0": -1, "y1": 1},
            "bodies": [{"type": "file", "path": f, "name": "volet"}]}
    write(OUT / f"{f}.json", json.dumps(case))
    p = run(tag, ["mesh", f"{f}.json", "-o", f"m_{f}", "--no-plot"], expect="mêmes cellules")
    res[-1]["cellules"] = next((ln for ln in p.stdout.splitlines() if "cellules" in ln), "")
# 4. maillages Gmsh 4.1 (ASCII, binaire) : carré de 2 triangles, 4 frontières nommées
MSH41 = """$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
1 1 "paroi"
1 2 "couvercle"
2 3 "fluide"
$EndPhysicalNames
$Entities
4 4 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 1 2 1 -2
2 1 0 0 1 1 0 1 1 2 2 -3
3 0 1 0 1 1 0 1 2 2 3 -4
4 0 0 0 0 1 0 1 1 2 4 -1
1 0 0 0 1 1 0 1 3 4 1 2 3 4
$EndEntities
$Nodes
5 4 1 4
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
2 1 0 0
$EndNodes
$Elements
5 6 1 6
1 1 1 1
1 1 2
1 2 1 1
2 2 3
1 3 1 1
3 3 4
1 4 1 1
4 4 1
2 1 2 2
5 1 2 3
6 1 3 4
$EndElements
"""
write(OUT / "carre41.msh", MSH41)
case = {"mesh": {"type": "file", "path": "carre41.msh"}, "physics": {"nu": 0.01},
        "boundary": {"paroi": {"type": "wall"}, "couvercle": {"type": "wall", "U": [1, 0]}},
        "solver": {"max_iter": 3}}
write(OUT / "carre41.json", json.dumps(case))
run("Gmsh 4.1 ASCII", ["run2d", "carre41.json", "-o", "o_msh41", "--no-plot"])
write(OUT / "carre41b.msh", MSH41.replace("4.1 0 8", "4.1 1 8"))
case["mesh"]["path"] = "carre41b.msh"
write(OUT / "carre41b.json", json.dumps(case))
run("Gmsh 4.1 binaire (en-tête)", ["run2d", "carre41b.json", "-o", "o_msh41b", "--no-plot"],
    expect="erreur claire")
# 5. cas JSON
from_toml = subprocess.run([sys.executable, "-c", "import json,sys;sys.path.insert(0,sys.argv[1]);"
                            "from microrans.tomlio import loads;print(json.dumps(loads(open(sys.argv[2],"
                            "encoding='utf-8').read())))", str(REPO), str(EX / "cavite_re100.toml")],
                           capture_output=True, text=True).stdout
write(OUT / "cavite.json", from_toml)
run("cas JSON", ["run2d", "cavite.json", "-o", "o_json", *SHORT])
(OUT / "c10.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
print("fin")
