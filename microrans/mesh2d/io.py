"""Entrées/sorties : maillages (Gmsh, SU2, VTK, OpenFOAM) et contours géométriques.

Maillages
  - Gmsh .msh : lecture ASCII v2.2 et v4.1, écriture v2.2 (groupes physiques = patches)
  - SU2 .su2  : lecture/écriture 2D (marqueurs = patches)
  - VTK .vtk  : écriture (ParaView), avec champs aux cellules
  - OpenFOAM  : écriture d'un constant/polyMesh extrudé d'une maille (faces avant/arrière 'empty')
Contours (objets importés)
  - .dat/.txt : profils Selig ou Lednicer, ou colonnes x y
  - .csv      : colonnes x, y (en-tête toléré)
  - .svg      : polygon, polyline, rect, circle, ellipse, path (M L H V C Q Z, absolus/relatifs)
  - .dxf      : LWPOLYLINE, POLYLINE/VERTEX, LINE, ARC, CIRCLE (ASCII), segments chaînés
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .mesh import Mesh2D

# ============================================================================ contours
def read_curve(path) -> np.ndarray:
    """Lit un contour fermé 2D (le plus grand si le fichier en contient plusieurs)."""
    curves = read_curves(path)
    if not curves:
        raise ValueError(f"Aucun contour trouvé dans {path}.")
    from .geometry import signed_area
    return max(curves, key=lambda c: abs(signed_area(c)) if len(c) > 2 else 0.0)


def read_curves(path) -> list[np.ndarray]:
    path = Path(path)
    ext = path.suffix.lower()
    text = path.read_text(encoding="utf-8", errors="replace")
    if ext in (".dat", ".txt", ".xy"):
        return [_read_airfoil_dat(text)]
    if ext == ".csv":
        return [_read_numeric_columns(text)]
    if ext == ".svg":
        return _read_svg(text)
    if ext == ".dxf":
        return _read_dxf(text)
    raise ValueError(f"Format de contour non supporté : {ext} (dat, txt, csv, svg, dxf)")


def _numbers(line):
    try:
        return [float(v) for v in re.split(r"[,\s;]+", line.strip()) if v]
    except ValueError:
        return None


def _read_numeric_columns(text):
    rows = []
    for line in text.splitlines():
        v = _numbers(line)
        if v and len(v) >= 2:
            rows.append(v[:2])
    if len(rows) < 3:
        raise ValueError("Moins de 3 points lus.")
    pts = np.array(rows)
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    return pts


def _read_airfoil_dat(text):
    lines = [ln for ln in text.splitlines() if ln.strip()]
    first = _numbers(lines[0])
    body = lines if (first and len(first) >= 2) else lines[1:]
    second = _numbers(body[0])
    # Lednicer : 1re ligne numérique = nombres de points extrados/intrados (> 1)
    if second and len(second) >= 2 and second[0] > 1.5 and second[1] > 1.5 \
            and float(second[0]).is_integer():
        nu, nl = int(second[0]), int(second[1])
        vals = [_numbers(ln) for ln in body[1:]]
        vals = [v[:2] for v in vals if v and len(v) >= 2]
        upper, lower = np.array(vals[:nu]), np.array(vals[nu:nu + nl])
        pts = np.vstack([upper[::-1], lower[1:]])
    else:
        vals = [v[:2] for v in (_numbers(ln) for ln in body) if v and len(v) >= 2]
        pts = np.array(vals)
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    # suppression des doublons consécutifs
    keep = np.concatenate([[True], np.linalg.norm(np.diff(pts, axis=0), axis=1) > 1e-12])
    return pts[keep]


def _parse_floats(s):
    return [float(v) for v in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", s)]


def _read_svg(text):
    curves = []
    for m in re.finditer(r"<(polygon|polyline)\b[^>]*\bpoints\s*=\s*\"([^\"]*)\"", text):
        v = _parse_floats(m.group(2))
        curves.append(np.array(v).reshape(-1, 2))
    for m in re.finditer(r"<rect\b([^>]*)>", text):
        a = _svg_attrs(m.group(1))
        x, y, w, h = (float(a.get(k, 0)) for k in ("x", "y", "width", "height"))
        curves.append(np.array([[x, y], [x + w, y], [x + w, y + h], [x, y + h]]))
    for m in re.finditer(r"<(circle|ellipse)\b([^>]*)>", text):
        a = _svg_attrs(m.group(2))
        cx, cy = float(a.get("cx", 0)), float(a.get("cy", 0))
        rx = float(a.get("r", a.get("rx", 0)))
        ry = float(a.get("r", a.get("ry", rx)))
        t = np.linspace(0, 2 * np.pi, 128, endpoint=False)
        curves.append(np.column_stack([cx + rx * np.cos(t), cy + ry * np.sin(t)]))
    for m in re.finditer(r"<path\b[^>]*\bd\s*=\s*\"([^\"]*)\"", text):
        curves += _svg_path(m.group(1))
    out = []
    for c in curves:
        c = np.asarray(c, float).copy()
        c[:, 1] = -c[:, 1]              # l'axe y du SVG pointe vers le bas
        if len(c) > 2 and np.allclose(c[0], c[-1]):
            c = c[:-1]
        if len(c) >= 3:
            out.append(c)
    return out


def _svg_attrs(s):
    return dict(re.findall(r"([\w-]+)\s*=\s*\"([^\"]*)\"", s))


def _svg_path(d):
    tokens = re.findall(r"[MmLlHhVvCcQqZz]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", d)
    curves, cur = [], []
    pos = np.zeros(2)
    start = np.zeros(2)
    i, cmd = 0, None
    while i < len(tokens):
        t = tokens[i]
        if re.match(r"[A-Za-z]", t):
            cmd = t
            i += 1
            if cmd in "Zz":
                if cur:
                    curves.append(np.array(cur))
                cur = []
                pos = start.copy()
                continue
        rel = cmd.islower()
        c = cmd.upper()
        need = {"M": 2, "L": 2, "H": 1, "V": 1, "C": 6, "Q": 4}[c]
        vals = np.array([float(v) for v in tokens[i:i + need]])
        i += need
        base = pos if rel else np.zeros(2)
        if c == "M":
            if cur:
                curves.append(np.array(cur))
            pos = base + vals
            start = pos.copy()
            cur = [pos.copy()]
            cmd = "l" if rel else "L"
        elif c == "L":
            pos = base + vals
            cur.append(pos.copy())
        elif c == "H":
            pos = np.array([vals[0] + (pos[0] if rel else 0.0), pos[1]])
            cur.append(pos.copy())
        elif c == "V":
            pos = np.array([pos[0], vals[0] + (pos[1] if rel else 0.0)])
            cur.append(pos.copy())
        elif c in "CQ":
            ctrl = [pos.copy()] + [base + vals[k:k + 2] for k in range(0, need, 2)]
            u = np.linspace(0, 1, 17)[1:, None]
            if c == "C":
                p0, p1, p2, p3 = ctrl
                pts = ((1 - u) ** 3 * p0 + 3 * (1 - u) ** 2 * u * p1 + 3 * (1 - u) * u ** 2 * p2
                       + u ** 3 * p3)
            else:
                p0, p1, p2 = ctrl
                pts = (1 - u) ** 2 * p0 + 2 * (1 - u) * u * p1 + u ** 2 * p2
            cur += list(pts)
            pos = ctrl[-1]
    if cur:
        curves.append(np.array(cur))
    return curves


def _read_dxf(text):
    lines = [ln.strip() for ln in text.splitlines()]
    pairs = list(zip(lines[0::2], lines[1::2]))
    entities, cur = [], None
    in_entities = False
    for code, val in pairs:
        if code == "2" and val == "ENTITIES":
            in_entities = True
        if not in_entities:
            continue
        if code == "0":
            if cur is not None:
                entities.append(cur)
            cur = {"type": val, "codes": []}
            if val == "ENDSEC":
                break
        elif cur is not None:
            cur["codes"].append((code, val))
    closed_curves, segments = [], []
    poly = None
    for e in entities:
        t, codes = e["type"], e["codes"]
        g = {}
        for c, v in codes:
            g.setdefault(c, []).append(v)
        if t == "LWPOLYLINE":
            pts = np.column_stack([np.array(g["10"], float), np.array(g["20"], float)])
            flags = int(g.get("70", ["0"])[0])
            if flags & 1:
                closed_curves.append(pts)
            else:
                segments += [pts[k:k + 2] for k in range(len(pts) - 1)]
        elif t == "POLYLINE":
            poly = {"pts": [], "closed": int(g.get("70", ["0"])[0]) & 1}
        elif t == "VERTEX" and poly is not None:
            poly["pts"].append([float(g["10"][0]), float(g["20"][0])])
        elif t == "SEQEND" and poly is not None:
            pts = np.array(poly["pts"])
            if poly["closed"]:
                closed_curves.append(pts)
            else:
                segments += [pts[k:k + 2] for k in range(len(pts) - 1)]
            poly = None
        elif t == "LINE":
            segments.append(np.array([[float(g["10"][0]), float(g["20"][0])],
                                      [float(g["11"][0]), float(g["21"][0])]]))
        elif t in ("CIRCLE", "ARC"):
            cx, cy, r = float(g["10"][0]), float(g["20"][0]), float(g["40"][0])
            if t == "CIRCLE":
                a = np.linspace(0, 2 * np.pi, 128, endpoint=False)
                closed_curves.append(np.column_stack([cx + r * np.cos(a), cy + r * np.sin(a)]))
            else:
                a0, a1 = np.radians(float(g["50"][0])), np.radians(float(g["51"][0]))
                if a1 <= a0:
                    a1 += 2 * np.pi
                a = np.linspace(a0, a1, max(int(32 * (a1 - a0) / np.pi), 4))
                segments.append(np.column_stack([cx + r * np.cos(a), cy + r * np.sin(a)]))
    closed_curves += _chain_segments(segments)
    return closed_curves


def _chain_segments(segments, tol=1e-8):
    """Assemble des segments/polylignes ouvertes en contours fermés."""
    segs = [np.asarray(s, float) for s in segments]
    loops = []
    while segs:
        chain = segs.pop(0)
        progress = True
        while progress and not np.allclose(chain[0], chain[-1], atol=tol):
            progress = False
            for k, s in enumerate(segs):
                if np.allclose(chain[-1], s[0], atol=tol):
                    chain = np.vstack([chain, s[1:]])
                elif np.allclose(chain[-1], s[-1], atol=tol):
                    chain = np.vstack([chain, s[::-1][1:]])
                else:
                    continue
                segs.pop(k)
                progress = True
                break
        if np.allclose(chain[0], chain[-1], atol=tol) and len(chain) > 3:
            loops.append(chain[:-1])
    return loops


# ============================================================================ maillages
def read_mesh(path, patch_types: dict | None = None) -> Mesh2D:
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".msh":
        return read_gmsh(path, patch_types)
    if ext == ".su2":
        return read_su2(path, patch_types)
    raise ValueError(f"Format de maillage non supporté en lecture : {ext} (msh, su2)")


def write_mesh(mesh: Mesh2D, path, cell_data: dict | None = None):
    """Écrit selon l'extension : .msh, .su2, .vtk ; un dossier (ou 'foam') -> polyMesh OpenFOAM."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".msh":
        write_gmsh(mesh, path)
    elif ext == ".su2":
        write_su2(mesh, path)
    elif ext == ".vtk":
        write_vtk(mesh, path, cell_data)
    elif ext in ("", ".foam"):
        write_openfoam(mesh, path.with_suffix(""))
    else:
        raise ValueError(f"Format non supporté en écriture : {ext} (msh, su2, vtk, dossier OpenFOAM)")


def _check_polygons(mesh, fmt):
    if np.any(mesh.cell_nv > 4):
        raise ValueError(f"{fmt} : seules les cellules triangles/quadrilatères sont exportables.")


# --- Gmsh -----------------------------------------------------------------------------
def write_gmsh(mesh: Mesh2D, path):
    _check_polygons(mesh, "Gmsh")
    patches = mesh.all_boundary_patches()
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$PhysicalNames",
             str(len(patches) + 1)]
    for k, (name, _, _) in enumerate(patches, start=1):
        lines.append(f'1 {k} "{name}"')
    fluid_tag = len(patches) + 1
    lines += [f'2 {fluid_tag} "fluid"', "$EndPhysicalNames", "$Nodes", str(mesh.n_points)]
    lines += [f"{i + 1} {x:.16g} {y:.16g} 0" for i, (x, y) in enumerate(mesh.points)]
    lines += ["$EndNodes", "$Elements"]
    elems = []
    eid = 1
    for k, (_, _, edges) in enumerate(patches, start=1):
        for a, b in edges:
            elems.append(f"{eid} 1 2 {k} {k} {a + 1} {b + 1}")
            eid += 1
    for row, nv in zip(mesh.cell_nodes, mesh.cell_nv):
        etype = 2 if nv == 3 else 3
        elems.append(f"{eid} {etype} 2 {fluid_tag} {fluid_tag} " + " ".join(str(v + 1) for v in row[:nv]))
        eid += 1
    lines += [str(len(elems))] + elems + ["$EndElements"]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_gmsh(path, patch_types: dict | None = None) -> Mesh2D:
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    sec = {}
    for m in re.finditer(r"\$(\w+)\s*\n(.*?)\n\$End\1", text, re.S):
        sec[m.group(1)] = m.group(2).split("\n")
    version = float(sec["MeshFormat"][0].split()[0])
    phys = {}
    if "PhysicalNames" in sec:
        for ln in sec["PhysicalNames"][1:]:
            parts = ln.split(maxsplit=2)
            if len(parts) == 3:
                phys[(int(parts[0]), int(parts[1]))] = parts[2].strip().strip('"')
    if version < 3:
        node_lines = sec["Nodes"][1:]
        tags = np.array([int(ln.split()[0]) for ln in node_lines])
        coords = np.array([[float(v) for v in ln.split()[1:3]] for ln in node_lines])
        elements = []
        for ln in sec["Elements"][1:]:
            v = [int(x) for x in ln.split()]
            etype, ntags = v[1], v[2]
            ptag = v[3] if ntags > 0 else 0
            elements.append((etype, ptag, v[3 + ntags:]))
    else:
        ent_phys = {}
        if "Entities" in sec:
            e = sec["Entities"]
            counts = [int(x) for x in e[0].split()]
            row = 1
            for dim, cnt in enumerate(counts):
                for _ in range(cnt):
                    v = e[row].split()
                    row += 1
                    tag = int(v[0])
                    k = 4 if dim == 0 else 7
                    nph = int(float(v[k]))
                    ent_phys[(dim, tag)] = int(v[k + 1]) if nph > 0 else 0
        nl = sec["Nodes"]
        nblocks = int(nl[0].split()[0])
        row = 1
        tags, coords = [], []
        for _ in range(nblocks):
            _, _, _, n = (int(x) for x in nl[row].split())
            row += 1
            tags += [int(nl[row + k]) for k in range(n)]
            coords += [[float(x) for x in nl[row + n + k].split()[:2]] for k in range(n)]
            row += 2 * n
        tags, coords = np.array(tags), np.array(coords)
        el = sec["Elements"]
        nblocks = int(el[0].split()[0])
        row = 1
        elements = []
        for _ in range(nblocks):
            edim, etag, etype, n = (int(x) for x in el[row].split())
            row += 1
            ptag = ent_phys.get((edim, etag), 0)
            for k in range(n):
                v = [int(x) for x in el[row + k].split()]
                elements.append((etype, ptag, v[1:]))
            row += n
    index = {t: i for i, t in enumerate(tags)}
    cells, boundary = [], {}
    for etype, ptag, nodes in elements:
        nodes = [index[t] for t in nodes]
        if etype == 1:
            name = phys.get((1, ptag), f"patch{ptag}")
            boundary.setdefault(name, []).append(nodes[:2])
        elif etype == 2:
            cells.append(nodes[:3])
        elif etype == 3:
            cells.append(nodes[:4])
    used = np.unique(np.concatenate([np.asarray(c) for c in cells]))
    remap = -np.ones(len(coords), dtype=int)
    remap[used] = np.arange(len(used))
    cells = [remap[np.asarray(c)] for c in cells]
    boundary = {k: remap[np.asarray(v)] for k, v in boundary.items()}
    return Mesh2D(coords[used], cells, boundary, patch_types)


# --- SU2 --------------------------------------------------------------------------------
def write_su2(mesh: Mesh2D, path):
    _check_polygons(mesh, "SU2")
    lines = ["NDIME= 2", f"NELEM= {mesh.n_cells}"]
    for i, (row, nv) in enumerate(zip(mesh.cell_nodes, mesh.cell_nv)):
        et = 5 if nv == 3 else 9
        lines.append(f"{et} " + " ".join(str(v) for v in row[:nv]) + f" {i}")
    lines.append(f"NPOIN= {mesh.n_points}")
    lines += [f"{x:.16g} {y:.16g} {i}" for i, (x, y) in enumerate(mesh.points)]
    patches = mesh.all_boundary_patches()
    lines.append(f"NMARK= {len(patches)}")
    for name, _, edges in patches:
        lines += [f"MARKER_TAG= {name}", f"MARKER_ELEMS= {len(edges)}"]
        lines += [f"3 {a} {b}" for a, b in edges]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_su2(path, patch_types: dict | None = None) -> Mesh2D:
    lines = [ln.split("%")[0].strip() for ln in
             Path(path).read_text(encoding="utf-8", errors="replace").splitlines()]
    lines = [ln for ln in lines if ln]
    i = 0
    cells, pts, boundary = [], [], {}
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("NDIME="):
            if int(ln.split("=")[1]) != 2:
                raise ValueError("Seuls les maillages SU2 2D sont supportés.")
            i += 1
        elif ln.startswith("NELEM="):
            n = int(ln.split("=")[1])
            for k in range(n):
                v = [int(x) for x in lines[i + 1 + k].split()]
                cells.append(v[1:4] if v[0] == 5 else v[1:5])
            i += n + 1
        elif ln.startswith("NPOIN="):
            n = int(ln.split("=")[1].split()[0])
            pts = [[float(x) for x in lines[i + 1 + k].split()[:2]] for k in range(n)]
            i += n + 1
        elif ln.startswith("NMARK="):
            nm = int(ln.split("=")[1])
            i += 1
            for _ in range(nm):
                name = lines[i].split("=")[1].strip()
                ne = int(lines[i + 1].split("=")[1])
                boundary[name] = [[int(x) for x in lines[i + 2 + k].split()[1:3]] for k in range(ne)]
                i += ne + 2
        else:
            i += 1
    return Mesh2D(np.array(pts), cells, boundary, patch_types)


# --- VTK ---------------------------------------------------------------------------------
def write_vtk(mesh: Mesh2D, path, cell_data: dict | None = None, title: str = "microrans"):
    """VTK legacy ASCII (UNSTRUCTURED_GRID), lisible par ParaView."""
    lines = ["# vtk DataFile Version 3.0", title, "ASCII", "DATASET UNSTRUCTURED_GRID",
             f"POINTS {mesh.n_points} double"]
    lines += [f"{x:.10g} {y:.10g} 0" for x, y in mesh.points]
    total = int(np.sum(mesh.cell_nv + 1))
    lines.append(f"CELLS {mesh.n_cells} {total}")
    lines += [f"{nv} " + " ".join(str(v) for v in row[:nv])
              for row, nv in zip(mesh.cell_nodes, mesh.cell_nv)]
    lines.append(f"CELL_TYPES {mesh.n_cells}")
    lines += [str({3: 5, 4: 9}.get(int(nv), 7)) for nv in mesh.cell_nv]
    if cell_data:
        lines.append(f"CELL_DATA {mesh.n_cells}")
        for name, arr in cell_data.items():
            arr = np.asarray(arr, float)
            safe = re.sub(r"\W", "_", name)
            if arr.ndim == 1:
                lines += [f"SCALARS {safe} double 1", "LOOKUP_TABLE default"]
                lines += [f"{v:.10g}" for v in arr]
            else:
                lines.append(f"VECTORS {safe} double")
                lines += [f"{a:.10g} {b:.10g} 0" for a, b in arr[:, :2]]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


# --- OpenFOAM ----------------------------------------------------------------------------
def _foam_header(cls, obj, note=None):
    s = ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
         f"    class       {cls};\n")
    if note:
        s += f'    note        "{note}";\n'
    s += f'    location    "constant/polyMesh";\n    object      {obj};\n}}\n\n'
    return s


def write_openfoam(mesh: Mesh2D, case_dir, thickness: float = 1.0):
    """Écrit constant/polyMesh (maillage 2D extrudé d'une maille en z, faces avant/arrière 'empty').

    Les paires périodiques sont écrites en patches 'cyclic'. Vérifiable avec `checkMesh`.
    """
    case_dir = Path(case_dir)
    d = case_dir / "constant" / "polyMesh"
    d.mkdir(parents=True, exist_ok=True)
    npt = mesh.n_points
    pts3 = np.vstack([np.column_stack([mesh.points, np.zeros(npt)]),
                      np.column_stack([mesh.points, np.full(npt, thickness)])])
    nreg = mesh.n_regular_internal
    faces, owner, neigh = [], [], []
    for (a, b), o, n in zip(mesh.face_nodes[:nreg], mesh.owner[:nreg], mesh.neighbour[:nreg]):
        faces.append((a, b, b + npt, a + npt))
        owner.append(o)
        neigh.append(n)
    boundary = []
    ptypes = {"wall": "wall", "patch": "patch", "symmetry": "symmetry", "empty": "empty"}
    patch_list = [(p.name, ptypes[p.type], mesh.face_nodes[p.faces], mesh.owner[p.faces], None)
                  for p in mesh.patches]
    for pa, pb, t in mesh.periodic_pairs:
        patch_list.append((pa, "cyclic", mesh._periodic_boundary[pa], mesh._periodic_owner[pa], pb))
        patch_list.append((pb, "cyclic", mesh._periodic_boundary[pb], mesh._periodic_owner[pb], pa))
    for name, ptype, fn, own, nb in patch_list:
        start = len(faces)
        for (a, b), o in zip(fn, own):
            faces.append((a, b, b + npt, a + npt))
            owner.append(o)
        boundary.append((name, ptype, start, len(fn), nb))
    start = len(faces)
    for c, (row, nv) in enumerate(zip(mesh.cell_nodes, mesh.cell_nv)):
        faces.append(tuple(int(v) for v in row[:nv][::-1]))          # z = 0, normale -z
        owner.append(c)
        faces.append(tuple(int(v) + npt for v in row[:nv]))          # z = h, normale +z
        owner.append(c)
    boundary.append(("frontAndBack", "empty", start, 2 * mesh.n_cells, None))
    note = f"nPoints:{2 * npt} nCells:{mesh.n_cells} nFaces:{len(faces)} nInternalFaces:{nreg}"
    (d / "points").write_text(_foam_header("vectorField", "points") + f"{len(pts3)}\n(\n"
                              + "\n".join(f"({x:.16g} {y:.16g} {z:.16g})" for x, y, z in pts3)
                              + "\n)\n", encoding="utf-8")
    (d / "faces").write_text(_foam_header("faceList", "faces") + f"{len(faces)}\n(\n"
                             + "\n".join(f"{len(f)}(" + " ".join(map(str, f)) + ")" for f in faces)
                             + "\n)\n", encoding="utf-8")
    (d / "owner").write_text(_foam_header("labelList", "owner", note) + f"{len(owner)}\n(\n"
                             + "\n".join(map(str, owner)) + "\n)\n", encoding="utf-8")
    (d / "neighbour").write_text(_foam_header("labelList", "neighbour", note) + f"{len(neigh)}\n(\n"
                                 + "\n".join(map(str, neigh)) + "\n)\n", encoding="utf-8")
    body = f"{len(boundary)}\n(\n"
    for name, ptype, start, n, nb in boundary:
        body += f"    {name}\n    {{\n        type            {ptype};\n"
        if ptype == "wall":
            body += "        inGroups        List<word> 1(wall);\n"
        if ptype == "cyclic":
            # transformation (translation) laissée au calcul automatique d'OpenFOAM
            body += f"        neighbourPatch  {nb};\n"
        body += f"        nFaces          {n};\n        startFace       {start};\n    }}\n"
    body += ")\n"
    (d / "boundary").write_text(_foam_header("polyBoundaryMesh", "boundary") + body, encoding="utf-8")
    return d
