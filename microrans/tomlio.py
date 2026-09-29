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
    try:
        import tomllib
    except ModuleNotFoundError:                    # Python 3.10
        import tomli as tomllib
    return tomllib.loads(text)
