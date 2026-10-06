"""Écriture TOML minimale (la bibliothèque standard ne sait que lire le TOML).

Couvre la structure des fichiers de cas : sections, sous-sections, tableaux de tables
([[bodies]]), tables et tableaux « en ligne ». Garantie testée : load(dump(cfg)) == cfg.
"""
from __future__ import annotations

import json
import math
import re

_BARE = re.compile(r"^[A-Za-z0-9_-]+$")


def _key(k: str) -> str:
    return k if _BARE.match(k) else json.dumps(k, ensure_ascii=False)


def _value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, dict):
        return "{ " + ", ".join(f"{_key(k)} = {_value(x)}" for k, x in v.items()
                                if x is not None) + " }"
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_value(x) for x in v) + "]"
    if hasattr(v, "tolist"):                       # tableaux numpy
        return _value(v.tolist())
    raise TypeError(f"Valeur non sérialisable en TOML : {v!r}")


def _is_table_of_tables(v) -> bool:
    return isinstance(v, dict) and len(v) > 0 and all(isinstance(x, dict) for x in v.values())


def dumps(cfg: dict) -> str:
    out: list[str] = []

    def table(path: list[str], d: dict):
        scalars = {k: v for k, v in d.items() if v is not None and not _is_table_of_tables(v)
                   and not (isinstance(v, list) and v and all(isinstance(x, dict) for x in v)
                            and len(path) == 0)}
        if path:
            out.append("")
            out.append("[" + ".".join(_key(p) for p in path) + "]")
        for k, v in scalars.items():
            if isinstance(v, dict) and len(path) == 0:
                continue                           # sections traitées plus bas
            out.append(f"{_key(k)} = {_value(v)}")
        for k, v in d.items():
            if len(path) == 0 and isinstance(v, dict):
                table([k], v)
            elif _is_table_of_tables(v):
                for sub, sv in v.items():
                    table(path + [k, sub], sv)
        if len(path) == 0:
            for k, v in d.items():
                if isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
                    for item in v:
                        out.append("")
                        out.append(f"[[{_key(k)}]]")
                        for kk, vv in item.items():
                            if vv is not None:
                                out.append(f"{_key(kk)} = {_value(vv)}")

    table([], cfg)
    return "\n".join(out).strip() + "\n"


def loads(text: str) -> dict:
    """TOML → dict ; erreur de syntaxe → ValueError en français (ligne, colonne)."""
    try:
        import tomllib
    except ModuleNotFoundError:                    # Python 3.10
        import tomli as tomllib
    try:
        return tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        m = re.search(r"\(at line (\d+), column (\d+)\)", str(exc))
        where = f"ligne {m.group(1)}, colonne {m.group(2)}" if m else "position inconnue"
        raise ValueError(f"syntaxe TOML incorrecte ({where}) : {exc}. Rappels : texte entre "
                         "guillemets droits (\"wall\"), point décimal (0.5), une clé par "
                         "ligne, sections entre crochets ([solver]).") from None


def read_text(path) -> tuple[str, str | None]:
    """Texte d'un fichier de cas et, s'il y a lieu, une remarque à afficher : UTF-8 avec
    BOM (Bloc-notes, PowerShell) lu sans remarque ; fichier qui n'est pas en UTF-8 (Latin-1
    / Windows-1252 d'anciens éditeurs) lu en Windows-1252, avec une remarque."""
    from pathlib import Path
    raw = Path(path).read_bytes()
    try:
        return raw.decode("utf-8-sig"), None
    except UnicodeDecodeError:
        return raw.decode("cp1252", errors="replace"), (
            f"{Path(path).name} n'est pas en UTF-8 : lu comme Latin-1 / Windows-1252 (accents "
            "des commentaires à vérifier) ; l'enregistrer en UTF-8.")
