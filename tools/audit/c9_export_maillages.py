"""Campagne 9 : maillages exportés puis relus. Pour des maillages de chaque type : écriture
.msh / .su2 / .vtk / OpenFOAM, relecture .msh et .su2 (cellules, aire totale, frontières,
nombre de faces par frontière), contrôle de structure du polyMesh OpenFOAM (propriétaire <
voisin, ordre triangulaire supérieur, normales sortantes du propriétaire, cellules fermées,
frontières contiguës couvrant toutes les faces de bord), sans OpenFOAM installé."""

import json
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from microrans.mesh2d.builder import build_mesh  # noqa: E402
from microrans.mesh2d.io import read_mesh, write_mesh, write_openfoam  # noqa: E402
from microrans.tomlio import loads  # noqa: E402

OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
EX = REPO / "microrans/examples"
CASES = {
    "rectangle": ("cavite_re100", {}),
    "periodique": ("canal_turbulent_3d", {"drop_extrude": True}),
    "ogrid": ("cylindre_re20", {}),
    "blocs": ("plaque_plane_sa", {}),
    "triangles": ("cylindre_re20", {"mesh": {"type": "unstructured", "h_max": 2.0,
                                             "h_surface": 0.1}}),
    "hybride": ("cylindre_re20", {"mesh": {"type": "hybrid", "h_max": 2.0, "h_surface": 0.1,
                                           "layers": {"n": 4, "first_height": 0.01,
                                                      "ratio": 1.2}}}),
}


def foam_list(path, kind):
    txt = Path(path).read_text()
    body = txt[txt.index("\n(", txt.index("FoamFile") + 8 if "FoamFile" in txt else 0):]
    m = re.search(r"(\d+)\s*\n\(", txt[txt.index("}") + 1:])
    n = int(m.group(1))
    inner = txt[txt.index("}") + 1:]
    inner = inner[inner.index("(", m.end() - 1) + 1: inner.rindex(")")]
    if kind == "label":
        return np.array(inner.split(), int)[:n]
    if kind == "point":
        return np.array(re.findall(r"\(([^()]*)\)", inner), dtype=object)
    del body
    return re.findall(r"\d+\(([^()]*)\)", inner)


def check_foam(d):
    pts = np.array([[float(v) for v in s.split()] for s in foam_list(d / "points", "point")])
    faces = [np.array(s.split(), int) for s in foam_list(d / "faces", "face")]
    owner = foam_list(d / "owner", "label")
    neigh = foam_list(d / "neighbour", "label")
    nint = len(neigh)
    ncell = int(max(owner.max(), neigh.max() if nint else 0)) + 1
    out = {"faces": len(faces), "internes": nint, "cellules": ncell}
    out["owner<neighbour"] = bool(np.all(owner[:nint] < neigh))
    key = owner[:nint].astype(np.int64) * (ncell + 1) + neigh
    out["ordre_triangulaire_sup"] = bool(np.all(np.diff(key) > 0))
    fc = np.array([pts[f].mean(0) for f in faces])
    area = np.zeros((len(faces), 3))
    for i, f in enumerate(faces):
        p = pts[f]
        area[i] = 0.5 * np.cross(p, np.roll(p, -1, axis=0)).sum(0)
    cc = np.zeros((ncell, 3))
    cnt = np.zeros(ncell)
    for i, f in enumerate(faces):
        cc[owner[i]] += fc[i]
        cnt[owner[i]] += 1
        if i < nint:
            cc[neigh[i]] += fc[i]
            cnt[neigh[i]] += 1
    cc /= cnt[:, None]
    d_on = np.einsum("ij,ij->i", area, fc - cc[owner])
    out["normales_sortantes_proprietaire"] = int(np.sum(d_on <= 0))
    if nint:
        d_nb = np.einsum("ij,ij->i", area[:nint], cc[neigh] - fc[:nint])
        out["normales_vers_voisin_fausses"] = int(np.sum(d_nb <= 0))
    closed = np.zeros((ncell, 3))
    np.add.at(closed, owner, area)
    np.add.at(closed, neigh, -area[:nint])
    scale = np.abs(area).max()
    out["cellules_non_fermees"] = int(np.sum(np.linalg.norm(closed, axis=1) > 1e-9 * scale))
    btxt = (d / "boundary").read_text()
    starts = [int(x) for x in re.findall(r"startFace\s+(\d+);", btxt)]
    nfs = [int(x) for x in re.findall(r"nFaces\s+(\d+);", btxt)]
    ok = starts and starts[0] == nint and all(starts[i] + nfs[i] == starts[i + 1]
                                               for i in range(len(starts) - 1)) \
        and starts[-1] + nfs[-1] == len(faces)
    out["frontieres_contigues"] = bool(ok)
    return out


def summary(m):
    return {"cellules": int(m.n_cells), "aire": float(np.sum(m.cell_volumes)),
            "frontieres": {p.name: (p.type, p.size) for p in m.patches},
            "periodiques": [list(p[:2]) for p in m.periodic_pairs]}


for label, (ex, mod) in CASES.items():
    cfg = loads((EX / f"{ex}.toml").read_text(encoding="utf-8"))
    if mod.get("drop_extrude"):
        cfg["mesh"].pop("extrude", None)
    if "mesh" in mod:
        cfg["mesh"] = {**{k: v for k, v in cfg["mesh"].items() if k in ("n_around",)},
                       **mod["mesh"]}
        cfg.setdefault("domain", {"type": "rectangle", "x0": -10, "x1": 20, "y0": -10, "y1": 10})
    r = {"cas": label}
    try:
        m = build_mesh(cfg, base_dir=EX)
        r["origine"] = summary(m)
        w = OUT / label
        w.mkdir(exist_ok=True)
        for ext in ("msh", "su2"):
            write_mesh(m, w / f"m.{ext}")
            back = read_mesh(w / f"m.{ext}")
            s = summary(back)
            r[f"relu_{ext}"] = {k: s[k] == r["origine"][k] for k in s}
            if not all(r[f"relu_{ext}"].values()):
                r[f"relu_{ext}_detail"] = s
        write_mesh(m, w / "m.vtk")
        write_openfoam(m, w / "foam")
        r["openfoam"] = check_foam(w / "foam" / "constant" / "polyMesh")
    except Exception:
        import traceback
        r["exception"] = traceback.format_exc()[-800:]
    with open(OUT / "c9.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False, default=str) + "\n")
    print(label, json.dumps({k: v for k, v in r.items() if k != "origine"}, ensure_ascii=False,
                            default=str)[:700], flush=True)
print("fin")
