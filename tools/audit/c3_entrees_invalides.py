"""Campagne 3 : entrées invalides ou limites (3D, clés récentes, commandes sur des cas 3D).
Chaque variante = un exemple modifié écrit en TOML, lancé comme un utilisateur."""

import copy
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from microrans.tomlio import dumps, loads  # noqa: E402

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
only = sys.argv[2:]
DEL = object()


def ex(name):
    return loads(
        (REPO / "microrans/examples" / f"{name}.toml").read_text(encoding="utf-8")
    )


def patch(cfg, **kv):
    cfg = copy.deepcopy(cfg)
    for path, v in kv.items():
        keys = path.split("__")
        d = cfg
        for k in keys[:-1]:
            if k.isdigit():
                d = d[int(k)]
            else:
                d = d.setdefault(k, {})
        if v is DEL:
            d.pop(keys[-1], None)
        else:
            d[keys[-1]] = v
    return cfg


BOX = ex("conduite_carree_3d")
EXT = ex("canal_turbulent_3d")
CAV3 = ex("cavite_cubique_re100_3d")
NACA = ex("naca0012_sa")
NACA3 = patch(
    NACA,
    mesh__extrude={
        "z0": 0,
        "z1": 0.2,
        "nz": 1,
        "patch_types": {"back": "symmetry", "front": "symmetry"},
    },
    boundary__farfield__U=[1.0, 0.0, 0.0],
    solver__max_iter=3,
)
for k in list(NACA3["boundary"]):
    if "U" in NACA3["boundary"][k] and len(NACA3["boundary"][k]["U"]) == 2:
        NACA3["boundary"][k]["U"] = NACA3["boundary"][k]["U"] + [0.0]
if "initial" in NACA3 and "U" in NACA3["initial"]:
    NACA3["initial"]["U"] = list(NACA3["initial"]["U"]) + [0.0]
SHORT = dict(solver__max_iter=3)

V = [
    # maillage box
    ("box_nz0", patch(BOX, mesh__nz=0, **SHORT), "run2d"),
    ("box_z1_lt_z0", patch(BOX, mesh__z1=-1.0, **SHORT), "run2d"),
    ("box_grading2", patch(BOX, mesh__grading=[1.0, 2.0], **SHORT), "run2d"),
    ("box_sans_z0", patch(BOX, mesh__z0=DEL, **SHORT), "run2d"),
    ("box_nz_texte", patch(BOX, mesh__nz="trente", **SHORT), "run2d"),
    ("box_extrude", patch(BOX, mesh__extrude={"nz": 2}, **SHORT), "run2d"),
    ("box_front_sans_cl", patch(BOX, boundary__front=DEL, **SHORT), "run2d"),
    ("box_names_faute", patch(BOX, mesh__names__fronte="front", **SHORT), "run2d"),
    (
        "box_periodic_faux",
        patch(BOX, mesh__periodic=[["inlet", "outlett"]], **SHORT),
        "run2d",
    ),
    (
        "box_periodic_non_apparie",
        patch(BOX, mesh__periodic=[["inlet", "top"]], **SHORT),
        "run2d",
    ),
    # extrusion
    ("ext_nz0", patch(EXT, mesh__extrude__nz=0, **SHORT), "run2d"),
    ("ext_z1_lt_z0", patch(EXT, mesh__extrude__z1=-1.0, **SHORT), "run2d"),
    ("ext_pas_table", patch(EXT, mesh__extrude=5, **SHORT), "run2d"),
    ("ext_grading0", patch(EXT, mesh__extrude__grading=0.0, **SHORT), "run2d"),
    (
        "ext_periodic_faute",
        patch(EXT, mesh__extrude__periodic=[["z0", "z2"]], **SHORT),
        "run2d",
    ),
    (
        "ext_patch_type_faute",
        patch(
            EXT,
            mesh__extrude__periodic=DEL,
            mesh__extrude__patch_types={"z0": "symetrie", "z1": "symmetry"},
            **SHORT,
        ),
        "run2d",
    ),
    ("ext_nz_dans_mesh", patch(EXT, mesh__nz=4, **SHORT), "run2d"),
    ("ext_z_sans_cl", patch(EXT, mesh__extrude__periodic=DEL, **SHORT), "run2d"),
    # vecteurs à 2 composantes en 3D
    ("U_cl_2comp", patch(CAV3, boundary__lid__U=[1.0, 0.0], **SHORT), "run2d"),
    ("U_init_2comp", patch(BOX, initial__U=[3.0, 0.0], **SHORT), "run2d"),
    ("body_force_2comp", patch(BOX, physics__body_force=[1.0, 0.0], **SHORT), "run2d"),
    ("sonde_2comp", patch(BOX, output__probes=[[0.25, 0.0]], **SHORT), "run2d"),
    (
        "ligne_2comp",
        patch(
            BOX,
            output__lines=[{"name": "l", "start": [0, 0], "end": [0.5, 0]}],
            **SHORT,
        ),
        "run2d",
    ),
    (
        "moment_center_2comp",
        patch(BOX, physics__moment_center=[0.0, 0.0], **SHORT),
        "run2d",
    ),
    (
        "sonde_hors_domaine",
        patch(BOX, output__probes=[[5.0, 5.0, 5.0]], **SHORT),
        "run2d",
    ),
    (
        "U_formule_z",
        patch(CAV3, boundary__lid__U=["sin(3.14159*z)", 0.0, 0.0], **SHORT),
        "run2d",
    ),
    ("U_formule_w", patch(CAV3, boundary__lid__U=["w*2", 0.0, 0.0], **SHORT), "run2d"),
    # options 2D seulement
    (
        "3d_compressible",
        patch(BOX, physics__compressible=True, physics__mach=0.5, **SHORT),
        "run2d",
    ),
    ("3d_couple", patch(BOX, solver__algorithm="coupled", **SHORT), "run2d"),
    ("3d_axisym", patch(BOX, physics__axisymmetric=True, **SHORT), "run2d"),
    ("3d_animation", patch(BOX, output__animate=True, **SHORT), "run2d"),
    (
        "3d_parabolique",
        patch(
            EXT,
            mesh__periodic=DEL,
            mesh__extrude__periodic=DEL,
            mesh__extrude__patch_types={"z0": "symmetry", "z1": "symmetry"},
            physics__body_force=DEL,
            boundary__inlet={"type": "inlet", "flow_rate": 1.0, "profile": "parabolic"},
            boundary__outlet={"type": "outlet"},
            **SHORT,
        ),
        "run2d",
    ),
    (
        "3d_reference_area_neg",
        patch(BOX, physics__reference_area=-1.0, **SHORT),
        "run2d",
    ),
    (
        "3d_gros_maillage",
        patch(BOX, mesh__nx=200, mesh__ny=200, mesh__nz=200, **SHORT),
        "run2d:timeout=25",
    ),
    # sortie
    ("vtk_format_faux", patch(CAV3, output__vtk_format="bin", **SHORT), "run2d"),
    # NACA : bord de fuite
    ("te_faux", patch(NACA, bodies__0__trailing_edge="pointu", **SHORT), "run2d"),
    (
        "te_sur_cercle",
        patch(ex("cylindre_re20"), bodies__0__trailing_edge="sharp", **SHORT),
        "run2d",
    ),
    (
        "venkat_k_neg",
        patch(ex("compressible_rampe_mach2"), solver__venkat_k=-1.0, **SHORT),
        "run2d",
    ),
    (
        "venkat_k_texte",
        patch(ex("compressible_rampe_mach2"), solver__venkat_k="haut", **SHORT),
        "run2d",
    ),
    # commandes sur des cas 3D
    ("naca3d_run", NACA3, "run2d"),
    ("naca3d_polaire", NACA3, "polar --alpha 0 4 4"),
    (
        "cav3d_balayage",
        patch(CAV3, **SHORT),
        "sweep --param physics.reynolds --values 50 100",
    ),
    ("cav3d_mesh_msh", CAV3, "mesh -f msh su2"),
    ("cav3d_mesh", CAV3, "mesh"),
]

env = dict(
    os.environ, PYTHONPATH=str(REPO), QT_QPA_PLATFORM="offscreen", MPLBACKEND="Agg"
)
for vid, cfg, cmd in V:
    if only and vid not in only:
        continue
    w = OUT / vid
    w.mkdir(exist_ok=True)
    cfg = copy.deepcopy(cfg)
    cfg.setdefault("output", {})["directory"] = str(w / "res")
    (w / "cas.toml").write_text(dumps(cfg), encoding="utf-8")
    timeout = 300
    if ":timeout=" in cmd:
        cmd, t = cmd.split(":timeout=")
        timeout = int(t)
    parts = shlex.split(cmd)
    argv = [
        sys.executable,
        "-m",
        "microrans",
        parts[0],
        *(["cas.toml"] if parts[0] != "mesh" else ["cas.toml", "-o", str(w / "mesh")]),
        *parts[1:],
    ]
    if parts[0] != "mesh":
        argv += ["--no-plot"]
    t0 = time.perf_counter()
    try:
        p = subprocess.run(
            argv, cwd=w, env=env, capture_output=True, text=True, timeout=timeout
        )
        rc, so, se = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        rc, so, se = (
            "timeout",
            (e.stdout or b"").decode()
            if isinstance(e.stdout, bytes)
            else (e.stdout or ""),
            (e.stderr or b"").decode()
            if isinstance(e.stderr, bytes)
            else (e.stderr or ""),
        )
    dt = round(time.perf_counter() - t0, 1)
    r = {
        "id": vid,
        "cmd": " ".join(argv[3:]),
        "rc": rc,
        "s": dt,
        "trace": "Traceback" in se or "Traceback" in so,
        "stdout": so[-1500:],
        "stderr": se[-1500:],
    }
    with open(OUT / "c3.jsonl", "a") as fh:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(vid, rc, dt, "TRACE" if r["trace"] else "", flush=True)
print("fin")
