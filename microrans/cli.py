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


def cmd_run2d(args) -> int:
    from .fv2d.case import run_case
    from .mesh2d.builder import load_config

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
    out = args.out or cfg.get("output", {}).get("directory") or f"results/{Path(args.case).stem}"
    summary = run_case(cfg, base_dir=Path(args.case).parent, out_dir=out,
                       verbose=not args.quiet, plot=not args.no_plot)
    return 0 if summary.get("converged", True) else 1


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
    p.add_argument("--steps-per-period", type=int, default=128)
    p.add_argument("--average", type=int, default=5, help="périodes moyennées (défaut : 5)")
    p.add_argument("--scheme", choices=list(TIME_SCHEMES_1D), default="bdf2",
                   help="schéma en temps (voir README) : implicites euler, bdf2, cn, sdirk2, "
                        "sdirk3 ; explicites rk1..rk4, ab2 (Δt très petit requis)")
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
    p.add_argument("case", help="fichier de cas .toml ou .json (voir dossier cases/)")
    p.add_argument("-o", "--out", help="dossier de sortie")
    p.add_argument("--set", nargs="+", action="extend", metavar="SECTION.CLE=VALEUR",
                   help="surcharge d'un paramètre, ex. physics.model=sst solver.max_iter=500")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_run2d)

    p = sub.add_parser("verify", help="vérification contre des solutions exactes")
    p.set_defaults(func=cmd_verify)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        print(f"Erreur : {exc}", file=sys.stderr)
        return 2
