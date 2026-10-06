"""Entrées/sorties : maillages (Gmsh, SU2, VTK, OpenFOAM) et contours géométriques.

Maillages
  - Gmsh .msh : lecture ASCII v2.2 et v4.1, écriture v2.2 (groupes physiques = patches)
  - SU2 .su2  : lecture/écriture 2D (marqueurs = patches)
  - VTK .vtk  : écriture (ParaView), avec champs aux cellules
  - OpenFOAM  : écriture d'un constant/polyMesh extrudé d'une maille (faces avant/arrière 'empty')
Contours (objets importés)
  - .dat/.txt : profils Selig ou Lednicer, ou colonnes x y
  - .csv      : colonnes x, y (en-tête toléré ; « x;y » à virgule décimale accepté)
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
    """Nombres d'une ligne : séparateurs virgule, espaces, tabulation ou point-virgule.
    Export d'un tableur français (« 0,5;1,25 », ou tabulations sans aucun point) : la
    virgule est décimale. Avant, « 1,000;-0,000 » donnait les points (1, 0), (0, 0)… sans
    message, et un contour de 100 points pouvait demander 14 Go au mailleur."""
    s = line.strip()
    if ";" in s or ("\t" in s and "," in s and "." not in s):
        parts = [v.replace(",", ".") for v in re.split(r"[;\s]+", s)]
    else:
        parts = re.split(r"[,\s]+", s)
    try:
        return [float(v) for v in parts if v]
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


def write_mesh(mesh: Mesh2D, path, cell_data: dict | None = None, binary: bool = True):
    """Écrit selon l'extension : .msh, .su2, .vtk ; un dossier (ou 'foam') -> polyMesh OpenFOAM.
    Maillage 3D : .vtk seulement. binary : .vtk binaire (défaut) ou texte."""
    path = Path(path)
    ext = path.suffix.lower()
    if getattr(mesh, "dim", 2) == 3 and ext != ".vtk":
        raise ValueError(f"Maillage 3D : export {ext or 'OpenFOAM'} non disponible, seulement "
                         ".vtk (ParaView).")
    if ext == ".msh":
        write_gmsh(mesh, path)
    elif ext == ".su2":
        write_su2(mesh, path)
    elif ext == ".vtk":
        write_vtk(mesh, path, cell_data, binary=binary)
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
    if "MeshFormat" not in sec:
        raise ValueError(f"{Path(path).name} : pas un fichier Gmsh .msh (section $MeshFormat "
                         "absente).")
    fmt = sec["MeshFormat"][0].split()
    version = float(fmt[0])
    if len(fmt) > 1 and fmt[1] != "0":              # M20 : avant KeyError 'Nodes'
        raise ValueError(f"{Path(path).name} : fichier Gmsh binaire, non lu. Le réexporter "
                         "en texte (Gmsh : File → Export, format .msh, option binaire "
                         "décochée).")
    missing = [k for k in ("Nodes", "Elements") if k not in sec]
    if missing:
        raise ValueError(f"{Path(path).name} : section(s) Gmsh {', '.join(missing)} absente(s) "
                         "(fichier incomplet ou tronqué).")
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
def write_vtk(mesh: Mesh2D, path, cell_data: dict | None = None, title: str = "microrans",
              binary: bool = True):
    """VTK legacy (UNSTRUCTURED_GRID), lisible par ParaView ; maillage 2D (z = 0) ou 3D
    (hexaèdres, prismes, pyramides, tétraèdres).

    binary = True (défaut) : données binaires gros-boutistes (format « BINARY » du VTK
    legacy), double précision, relecture exacte ; binary = False : texte, 10 chiffres
    significatifs. Mesuré (10⁶ hexaèdres, U + 4 champs scalaires, machine de test) : texte
    14.0 à 15.1 s et 173 Mo, binaire 0.31 à 0.35 s et 121 Mo (README § 7)."""
    dim = getattr(mesh, "dim", 2)
    pts = mesh.points if dim == 3 else np.column_stack([mesh.points, np.zeros(mesh.n_points)])
    nv_all = np.asarray(mesh.cell_nv)
    if dim == 3:
        from ..mesh3d.mesh import VTK_TYPES as types
    else:
        types = {3: 5, 4: 9}                         # triangle, quadrilatère ; sinon polygone
    if not binary:
        _write_vtk_ascii(mesh, path, cell_data, title, pts, nv_all, types, dim)
        return
    nc = mesh.n_cells
    cn = np.asarray(mesh.cell_nodes)
    keep = np.arange(cn.shape[1])[None, :] < nv_all[:, None]
    conn = np.column_stack([nv_all, cn])[np.column_stack([np.ones(nc, bool), keep])]
    if conn.size >= 2 ** 31:
        raise ValueError("Maillage trop grand pour le format VTK legacy (entiers 32 bits).")
    lut = {int(k): types.get(int(k), 7) for k in np.unique(nv_all)}
    ctypes = np.array([lut[int(k)] for k in nv_all], dtype=">i4") if nc else np.zeros(0, ">i4")

    def block(f, head, arr, dtype):
        f.write(head.encode("utf-8") + b"\n")        # noms accentués : comme le texte
        f.write(np.ascontiguousarray(arr, dtype=dtype).tobytes())
        f.write(b"\n")

    with open(path, "wb") as f:
        f.write(f"# vtk DataFile Version 3.0\n{title}\nBINARY\nDATASET UNSTRUCTURED_GRID\n"
                .encode("utf-8"))
        block(f, f"POINTS {mesh.n_points} double", pts, ">f8")
        block(f, f"CELLS {nc} {conn.size}", conn, ">i4")
        block(f, f"CELL_TYPES {nc}", ctypes, ">i4")
        if cell_data:
            f.write(f"CELL_DATA {nc}\n".encode("utf-8"))
            for name, arr in cell_data.items():
                arr = np.asarray(arr, float)
                safe = re.sub(r"\W", "_", name)
                if arr.ndim == 1:
                    block(f, f"SCALARS {safe} double 1\nLOOKUP_TABLE default", arr, ">f8")
                else:
                    v3 = arr[:, :3] if arr.shape[1] >= 3 else np.column_stack(
                        [arr, np.zeros((len(arr), 3 - arr.shape[1]))])
                    block(f, f"VECTORS {safe} double", v3, ">f8")


def _write_vtk_ascii(mesh, path, cell_data, title, pts, nv_all, types, dim):
    """Format texte (avant le lot E3 : seul format ; inchangé)."""
    lines = ["# vtk DataFile Version 3.0", title, "ASCII", "DATASET UNSTRUCTURED_GRID",
             f"POINTS {mesh.n_points} double"]
    lines += [f"{x:.10g} {y:.10g} {z:.10g}" for x, y, z in pts]
    total = int(np.sum(nv_all + 1))
    lines.append(f"CELLS {mesh.n_cells} {total}")
    if dim == 3:
        rows = mesh.cells_as_lists()
    else:
        rows = [row[:nv] for row, nv in zip(mesh.cell_nodes, nv_all)]
    lines += [f"{len(row)} " + " ".join(str(v) for v in row) for row in rows]
    lines.append(f"CELL_TYPES {mesh.n_cells}")
    lines += [str(types.get(int(nv), 7)) for nv in nv_all]
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
                v3 = arr[:, :3] if arr.shape[1] >= 3 else np.column_stack(
                    [arr, np.zeros((len(arr), 3 - arr.shape[1]))])
                lines += [f"{a:.10g} {b:.10g} {c:.10g}" for a, b, c in v3]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_vtk(path) -> dict:
    """Relit un fichier VTK legacy UNSTRUCTURED_GRID écrit par `write_vtk` (texte ou binaire) :
    {"title", "binary", "points" (n, 3), "cells" (connectivité à plat : n, i0, … par
    cellule), "cell_types", "cell_data" {nom: tableau}}. Lecteur minimal (tests, contrôle des
    sorties) ; ParaView lit les deux formats."""
    raw = Path(path).read_bytes()
    pos = 0

    def line() -> str:
        nonlocal pos
        end = raw.find(b"\n", pos)
        end = len(raw) if end < 0 else end
        out = raw[pos:end].decode("utf-8", errors="replace").strip()
        pos = end + 1
        return out

    line()                                           # « # vtk DataFile Version … »
    res = {"title": line(), "cell_data": {}}
    fmt = line().upper()
    if fmt not in ("ASCII", "BINARY"):
        raise ValueError(f"{path} : format VTK {fmt!r} inconnu (ASCII ou BINARY attendu).")
    binary = res["binary"] = fmt == "BINARY"

    def values(count, big_endian, kind):
        nonlocal pos
        if binary:
            dt = np.dtype(big_endian)
            arr = np.frombuffer(raw, dtype=dt, count=count, offset=pos).astype(kind)
            pos += count * dt.itemsize
            if raw[pos:pos + 1] == b"\n":
                pos += 1
            return arr
        toks: list[str] = []
        while len(toks) < count:
            toks += line().split()
        return np.array(toks[:count], dtype=kind)

    nc = 0
    while pos < len(raw):
        words = line().split()
        if not words or words[0] == "DATASET":
            continue
        kw = words[0]
        if kw == "POINTS":
            res["points"] = values(3 * int(words[1]), ">f8", float).reshape(-1, 3)
        elif kw == "CELLS":
            res["cells"] = values(int(words[2]), ">i4", np.int64)
        elif kw == "CELL_TYPES":
            res["cell_types"] = values(int(words[1]), ">i4", np.int64)
        elif kw == "CELL_DATA":
            nc = int(words[1])
        elif kw == "SCALARS":
            line()                                   # LOOKUP_TABLE default
            res["cell_data"][words[1]] = values(nc, ">f8", float)
        elif kw == "VECTORS":
            res["cell_data"][words[1]] = values(3 * nc, ">f8", float).reshape(-1, 3)
        else:
            raise ValueError(f"{path} : section VTK {kw!r} non prise en charge.")
    return res


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
