"""Interface en ligne de commande : python -m microrans {rans,urans,verify} ..."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .models import MODELS, TURBULENT_MODELS, canonical_name
from .solver import TIME_SCHEMES_1D


def _model_list(values: list[str]) -> list[str]:
    out = []
    for v in values:
        for name in v.split(","):
            if not name:
                continue
            if name.lower() == "all":
                out.extend(TURBULENT_MODELS)
            else:
                out.append(canonical_name(name))
    return list(dict.fromkeys(out))


def _add_common(p: argparse.ArgumentParser, default_out: str):
    p.add_argument("-m", "--model", nargs="+", default=["sa"],
                   help=f"modèle(s) : {', '.join(MODELS)} ou 'all' (défaut : sa)")
    p.add_argument("--re-tau", type=float, default=395.0,
                   help="Reynolds de frottement nominal Re_τ = u_τ h/ν (défaut : 395)")
    p.add_argument("--n-cells", type=int, default=192, help="nombre de mailles (défaut : 192)")
    p.add_argument("--y1plus", type=float, default=0.2,
                   help="hauteur de la 1re maille en unités de paroi (défaut : 0.2)")
    p.add_argument("--sa-ft2", action="store_true",
                   help="SA : active le terme f_t2 (défaut : SA-noft2)")
    p.add_argument("-o", "--out", default=default_out, help=f"dossier de sortie (défaut : {default_out})")
    p.add_argument("--no-plot", action="store_true", help="ne pas produire de figures")
    p.add_argument("-q", "--quiet", action="store_true", help="moins de messages")


def _model_options(name: str, args) -> dict:
    return {"ft2": True} if (name == "sa" and args.sa_ft2) else {}


def cmd_rans(args) -> int:
    from .cases import run_rans_channel
    from .postprocess import save_rans, save_rans_comparison
    from .reference import load_reference_profile

    models = _model_list(args.model)
    out = Path(args.out)
    reference = load_reference_profile(args.reference) if args.reference else None
    results = []
    status = 0
    for name in models:
        if not args.quiet:
            print(f"RANS {name} : Re_τ = {args.re_tau:g}, {args.n_cells} mailles, y1+ = {args.y1plus:g}")
        res = run_rans_channel(name, re_tau=args.re_tau, n_cells=args.n_cells,
                               y1_plus=args.y1plus, model_options=_model_options(name, args),
                               dt=args.dt, relax=args.relax, max_iter=args.max_iter,
                               tol=args.tol, verbose=args.verbose)
        s = save_rans(res, out / name, plot=not args.no_plot)
        results.append(res)
        if not s["converged"]:
            status = 1
        if not args.quiet:
            print(f"  -> {'convergé' if s['converged'] else 'NON CONVERGÉ'} en {s['iterations']} it. "
                  f"({s['wall_time_s']:.2f} s) | U_b+ = {s['Ub_plus']:.3f} "
                  f"(Dean : {s['Ub_plus_dean']:.3f}, écart {s['Ub_plus_error_vs_dean_pct']:+.1f} %) | "
                  f"τ_w = {s['tau_wall_bottom']:.5f}")
    if len(results) > 1 or reference is not None:
        save_rans_comparison(results, out, plot=not args.no_plot, reference=reference)
    if not args.quiet:
        print(f"Résultats dans {out.resolve()}")
    return status


def cmd_urans(args) -> int:
    from .cases import run_pulsating_channel
    from .postprocess import save_urans, save_urans_comparison

    models = _model_list(args.model)
    out = Path(args.out)
    results = []
    for name in models:
        if not args.quiet:
            print(f"URANS {name} : Re_τ = {args.re_tau:g}, ω⁺ = {args.omega_plus:g}, "
                  f"A = {args.amplitude:g}, schéma {args.scheme}")
        res = run_pulsating_channel(
            name, re_tau=args.re_tau, n_cells=args.n_cells, y1_plus=args.y1plus,
            omega_plus=args.omega_plus, amplitude=args.amplitude, n_periods=args.periods,
            t_transient=args.t_transient, steps_per_period=args.steps_per_period,
            n_average=args.average, scheme=args.scheme, max_inner=args.max_inner,
            inner_tol=args.inner_tol, relax=args.relax,
            model_options=_model_options(name, args), verbose=args.verbose)
        s = save_urans(res, out / name, plot=not args.no_plot)
        results.append(res)
        if not args.quiet:
            print(f"  -> {s['n_periods']} périodes, {s['mean_inner_iterations']:.1f} sous-it./pas "
                  f"({s['wall_time_s']:.1f} s) | ⟨τ_w⟩ = {s['mean_tau_wall']:.5f} (attendu 1) | "
                  f"|τ̂_w| = {s['tau_wall_amplitude']:.4f}, phase {s['tau_wall_phase_deg']:.1f}° | "
                  f"périodicité {s['periodicity_error']:.1e}")
    if len(results) > 1:
        save_urans_comparison(results, out, plot=not args.no_plot)
    if not args.quiet:
        print(f"Résultats dans {out.resolve()}")
    return 0


def cmd_verify(args) -> int:
    from .verification import run_all
    ok = run_all(verbose=True)
    print("Vérification : " + ("tout est OK" if ok else "ÉCHEC"))
    return 0 if ok else 1


def cmd_mesh(args) -> int:
    import json

    from .mesh2d.builder import PRESETS, build_mesh, load_config
    from .mesh2d.io import write_mesh
    from .mesh2d.plot import plot_mesh

    if args.preset:
        if args.preset not in PRESETS:
            raise ValueError(f"Préréglage inconnu '{args.preset}'. Choix : {', '.join(PRESETS)}")
        mesh = PRESETS[args.preset]()
        name = args.preset
    elif args.config:
        args.config = str(resolve_example(args.config))
        cfg = load_config(args.config)
        if args.type:
            cfg.setdefault("mesh", {})["type"] = args.type
        mesh = build_mesh(cfg, base_dir=Path(args.config).parent, verbose=args.verbose)
        name = Path(args.config).stem
    else:
        raise ValueError("Donner un fichier de configuration ou --preset.")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    q = mesh.quality()
    (out / "quality.json").write_text(json.dumps(q, indent=2, ensure_ascii=False), encoding="utf-8")
    for fmt in args.format:
        target = out / name if fmt == "foam" else out / f"{name}.{fmt}"
        write_mesh(mesh, target if fmt != "foam" else target.with_suffix(".foam"))
    if not args.no_plot:
        plot_mesh(mesh, out / f"{name}.png")
        walls = [p.name for p in mesh.patches if p.type == "wall"]
        if walls:
            import numpy as np
            pts = np.vstack([mesh.points[mesh.patch_face_nodes(w)].reshape(-1, 2) for w in walls])
            (x0, y0), (x1, y1) = pts.min(axis=0), pts.max(axis=0)
            L = max(x1 - x0, y1 - y0)
            plot_mesh(mesh, out / f"{name}_zoom.png",
                      zoom=(x0 - 0.3 * L, x1 + 0.3 * L, y0 - 0.3 * L, y1 + 0.3 * L))
    if not args.quiet:
        types = ", ".join(f"{v} {k}" for k, v in q["cell_types"].items())
        print(f"Maillage '{name}' : {q['n_cells']} cellules ({types}), {q['n_points']} sommets")
        for pn, pv in q["patches"].items():
            print(f"  patch {pn:16s} {pv['type']:9s} {pv['faces']} faces")
        print(f"  non-orthogonalité max {q['non_orthogonality_max_deg']:.1f}° "
              f"(moy. {q['non_orthogonality_mean_deg']:.1f}°), asymétrie max {q['skewness_max']:.2f}, "
              f"rapport d'aspect max {q['aspect_ratio_max']:.0f}")
        for w in mesh.check():
            print(f"  ATTENTION : {w}")
        print(f"Fichiers dans {out.resolve()}")
    return 0


def _load_case(args):
    from .mesh2d.builder import load_config

    args.case = str(resolve_example(args.case))
    cfg = load_config(args.case)
    for item in args.set or []:
        key, _, val = item.partition("=")
        sec, _, name = key.partition(".")
        try:
            import json as _json
            val = _json.loads(val)
        except ValueError:
            pass
        cfg.setdefault(sec, {})[name] = val
    return cfg


def cmd_run2d(args) -> int:
    from .fv2d.case import run_case

    cfg = _load_case(args)
    out = args.out or cfg.get("output", {}).get("directory") or f"results/{Path(args.case).stem}"
    if args.continue_run:
        args.restart = str(Path(out) / "checkpoint.npz")
    if args.restart:
        cfg.setdefault("initial", {})["restart"] = str(Path(args.restart).resolve())
    summary = run_case(cfg, base_dir=Path(args.case).parent, out_dir=out,
                       verbose=not args.quiet, plot=not args.no_plot)
    return 0 if summary.get("converged", True) else 1


def cmd_sweep(args) -> int:
    from .fv2d.sweep import ALPHA, run_sweep

    cfg = _load_case(args)
    sw = cfg.get("sweep", {})
    key = ALPHA if getattr(args, "alpha", None) else (args.param or sw.get("parameter"))
    values = (args.alpha or args.range or args.values
              or sw.get("values") or sw.get("range"))
    if not key or values is None:
        print("Indiquer le paramètre et les valeurs : --param physics.reynolds --values 10 20 "
              "(ou --range début fin pas), ou une section [sweep] dans le cas.")
        return 2
    from_range = bool(args.range or args.alpha) or (
        values is sw.get("range") and values is not None)
    if isinstance(values, list) and len(values) == 3 and from_range:
        values = f"{values[0]}:{values[1]}:{values[2]}"
    cont = sw.get("continuation", True) and not args.no_continuation
    out = args.out or cfg.get("output", {}).get("directory") or f"results/{Path(args.case).stem}"
    rows = run_sweep(cfg, key, values, base_dir=Path(args.case).parent, out_dir=out,
                     continuation=cont, verbose=not args.quiet, plot=not args.no_plot)
    cols = [key, "converged", "iterations"] + [k for k in rows[0] if k.split("_")[0] in
                                               ("Cl", "Cd", "Cm") and "_" in k
                                               and not k.startswith(("Cd_p", "Cd_v"))]
    print("  ".join(f"{c:>14s}" for c in cols))
    for r in rows:
        print("  ".join(f"{r[c]:>14.5g}" if isinstance(r[c], float) else f"{str(r[c]):>14s}"
                        for c in cols))
    return 0 if all(r["converged"] for r in rows) else 1


def cmd_schemes(args) -> int:
    from .studies import markdown_table, plot_time_study, time_study_1d, time_study_2d
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print("Étude 1D (canal turbulent pulsé, SA)...")
    r1 = time_study_1d()
    print(markdown_table(r1, ["scheme", "steps", "error", "cpu", "inner"]))
    print("Étude 2D (tourbillon de Taylor-Green advecté)...")
    r2 = time_study_2d()
    print(markdown_table(r2, ["scheme", "dt", "courant", "error", "cpu"]))
    plot_time_study(r1, r2, out / "schemas_temps.png")
    (out / "schemas_temps.md").write_text(
        "## 1D\n\n" + markdown_table(r1, ["scheme", "steps", "error", "cpu", "inner"])
        + "\n## 2D\n\n" + markdown_table(r2, ["scheme", "dt", "courant", "error", "cpu"]),
        encoding="utf-8")
    print(f"Figure et tableaux dans {out.resolve()}")
    return 0


def examples_dir() -> Path:
    """Dossier des cas d'exemple (paquet installé ou exécutable PyInstaller)."""
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    for cand in (Path(__file__).resolve().parent / "examples", base / "microrans" / "examples"):
        if cand.is_dir():
            return cand
    return Path(__file__).resolve().parent / "examples"


def resolve_example(name) -> Path:
    """Chemin d'un fichier de cas ; à défaut, nom d'un exemple fourni (avec ou sans .toml)."""
    p = Path(name)
    if p.exists():
        return p
    for cand in (examples_dir() / p.name, examples_dir() / (p.name + ".toml")):
        if cand.exists():
            return cand
    raise ValueError(f"Fichier de cas introuvable : {name} (exemples : microrans examples)")


def cmd_examples(args) -> int:
    d = examples_dir()
    print(f"Exemples fournis ({d}) :")
    for f in sorted(d.glob("*.toml")):
        first = f.read_text(encoding="utf-8").splitlines()[0].lstrip("# ").strip()
        print(f"  {f.stem:28s} {first}")
    print("Lancer : microrans run2d <nom>   (ou microrans mesh <nom> pour mesh_*)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="microrans",
        description="Micro-solveur RANS/URANS 1D/2D : SA, k-ε, k-ω, k-ω SST ; mailleur 2D.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("rans", help="canal plan établi, stationnaire")
    _add_common(p, "results/rans")
    p.add_argument("--dt", type=float, default=5.0, help="pas de pseudo-temps (défaut : 5)")
    p.add_argument("--relax", type=float, default=0.5,
                   help="sous-relaxation des variables de turbulence (défaut : 0.5)")
    p.add_argument("--max-iter", type=int, default=20000)
    p.add_argument("--tol", type=float, default=1e-10,
                   help="variation relative max par itération pour l'arrêt (défaut : 1e-10)")
    p.add_argument("--reference", help="fichier (y+, U+) à superposer, ex. profil DNS")
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(func=cmd_rans)

    p = sub.add_parser("urans", help="canal à gradient de pression pulsé, instationnaire")
    _add_common(p, "results/urans")
    p.add_argument("--omega-plus", type=float, default=0.01,
                   help="pulsation en unités de paroi ω⁺ = ων/u_τ² (défaut : 0.01)")
    p.add_argument("--amplitude", type=float, default=10.0,
                   help="amplitude A du forçage f = 1 + A sin(ωt) (défaut : 10)")
    p.add_argument("--periods", type=int, default=None,
                   help="nombre total de périodes (défaut : auto via --t-transient)")
    p.add_argument("--t-transient", type=float, default=80.0,
                   help="durée de transitoire avant moyenne, en h/u_τ (défaut : 80)")
    p.add_argument("--steps-per-period", type=int, default=64,
                   help="pas par période (défaut 64 : erreur en temps ~1e-5 avec sdirk2)")
    p.add_argument("--average", type=int, default=5, help="périodes moyennées (défaut : 5)")
    p.add_argument("--scheme", choices=list(TIME_SCHEMES_1D), default="sdirk2",
                   help="schéma en temps (défaut sdirk2, voir README) : implicites euler, bdf2, "
                        "cn, sdirk2, sdirk3 ; explicites rk1..rk4, ab2 (Δt minuscule requis)")
    p.add_argument("--max-inner", type=int, default=30)
    p.add_argument("--inner-tol", type=float, default=1e-6)
    p.add_argument("--relax", type=float, default=1.0,
                   help="sous-relaxation dans les sous-itérations (défaut : 1)")
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(func=cmd_urans)

    p = sub.add_parser("mesh", help="générer / convertir un maillage 2D")
    p.add_argument("config", nargs="?", help="fichier de configuration .toml ou .json")
    p.add_argument("--preset", help="préréglage : cavity, channel, backstep, flatplate, "
                   "cylinder-ogrid, cylinder-tri, cylinder-hybrid, naca0012-ogrid, naca0012-hybrid")
    p.add_argument("--type", help="force le type : blocks, rectangle, ogrid, unstructured, hybrid, file")
    p.add_argument("-f", "--format", nargs="+", default=["msh", "vtk"],
                   choices=["msh", "su2", "vtk", "foam"], help="formats de sortie (défaut : msh vtk)")
    p.add_argument("-o", "--out", default="results/mesh")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(func=cmd_mesh)

    p = sub.add_parser("run2d", help="calcul 2D (RANS/URANS) décrit par un fichier de cas")
    p.add_argument("case", help="fichier de cas .toml / .json, ou nom d'un exemple "
                                "(liste : microrans examples)")
    p.add_argument("-o", "--out", help="dossier de sortie")
    p.add_argument("--set", nargs="+", action="extend", metavar="SECTION.CLE=VALEUR",
                   help="surcharge d'un paramètre, ex. physics.model=sst solver.max_iter=500")
    p.add_argument("--restart", metavar="FICHIER",
                   help="repartir d'un fichier checkpoint.npz (même maillage : reprise exacte ; "
                        "autre maillage : champs interpolés)")
    p.add_argument("--continue", dest="continue_run", action="store_true",
                   help="poursuivre le calcul depuis le checkpoint.npz du dossier de sortie "
                        "(stationnaire : max_iter itérations de plus ; instationnaire : "
                        "jusqu'au nouveau t_end)")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_run2d)

    def sweep_common(p):
        p.add_argument("case", help="fichier de cas .toml / .json, ou nom d'un exemple")
        p.add_argument("-o", "--out", help="dossier de sortie")
        p.add_argument("--set", nargs="+", action="extend", metavar="SECTION.CLE=VALEUR")
        p.add_argument("--no-continuation", action="store_true",
                       help="chaque point part de l'état initial (plus lent, pas d'hystérésis)")
        p.add_argument("--no-plot", action="store_true")
        p.add_argument("-q", "--quiet", action="store_true")
        p.set_defaults(func=cmd_sweep, param=None, values=None, range=None, alpha=None)

    p = sub.add_parser("polar", help="polaire Cl(α), Cd(α), Cm(α) (incidence de l'écoulement)")
    sweep_common(p)
    p.add_argument("--alpha", nargs=3, type=float, metavar=("DEBUT", "FIN", "PAS"),
                   required=True, help="incidences en degrés, ex. --alpha -4 12 2")
    p = sub.add_parser("sweep", help="balayage d'un paramètre quelconque du cas")
    sweep_common(p)
    p.add_argument("--param", help="clé du cas, ex. physics.reynolds, physics.angle_of_attack")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--values", nargs="+", type=float, help="valeurs explicites")
    g.add_argument("--range", nargs=3, type=float, metavar=("DEBUT", "FIN", "PAS"))

    p = sub.add_parser("examples", help="liste des cas d'exemple fournis")
    p.set_defaults(func=cmd_examples)

    p = sub.add_parser("gui", help="interface graphique (nécessite PySide6)")
    p.set_defaults(func=lambda a: __import__("microrans.gui", fromlist=["main"]).main([]))

    p = sub.add_parser("schemes", help="étude précision / coût des schémas en temps")
    p.add_argument("-o", "--out", default="results/schemes")
    p.set_defaults(func=cmd_schemes)

    p = sub.add_parser("verify", help="vérification contre des solutions exactes")
    p.set_defaults(func=cmd_verify)
    return parser


def _utf8_console():
    """Consoles Windows (cp1252...) : sortie UTF-8 tolérante au lieu d'un plantage sur
    les symboles (≈, ν, τ, ⁺...)."""
    for stream in (sys.stdout, sys.stderr):
        enc = (getattr(stream, "encoding", None) or "").lower().replace("-", "")
        if stream is not None and enc != "utf8" and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except (ValueError, OSError):
                pass


def main(argv=None) -> int:
    _utf8_console()
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        print(f"Erreur : {exc}", file=sys.stderr)
        return 2
