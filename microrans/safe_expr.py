"""Évaluation SÛRE des formules des fichiers de cas (profils d'entrée, valeurs initiales,
sources… : « 6*y*(1-y) », « where(y < 0.5, 1, 0) »).

Un fichier de cas peut venir de n'importe qui : il ne doit pas pouvoir exécuter de code.
`eval` de Python, même privé de `__builtins__`, se contourne (introspection des objets,
ex. `().__class__.__base__.__subclasses__()`). Ici la formule est analysée en arbre
syntaxique puis évaluée nœud par nœud, avec une liste blanche : nombres, variables
autorisées, opérateurs arithmétiques et de comparaison, fonctions mathématiques de NumPy.
Tout le reste (attributs, indices, appels quelconques, lambda, chaînes…) est refusé.
Les nombres sont convertis en flottants : `9**9**9` déborde immédiatement au lieu de
bloquer le programme sur un calcul d'entier géant.
"""
from __future__ import annotations

import ast
import operator

import numpy as np

MAX_LENGTH = 2000

FUNCTIONS = {name: getattr(np, name) for name in (
    "sqrt", "exp", "log", "log10", "sin", "cos", "tan", "arcsin", "arccos", "arctan",
    "arctan2", "sinh", "cosh", "tanh", "abs", "minimum", "maximum", "where", "hypot", "clip",
    "sign", "floor", "ceil")}
CONSTANTS = {"pi": np.pi, "e": np.e}

_BIN = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.Pow: operator.pow, ast.Mod: operator.mod,
        ast.FloorDiv: operator.floordiv, ast.BitAnd: operator.and_, ast.BitOr: operator.or_}
_UNARY = {ast.UAdd: operator.pos, ast.USub: operator.neg, ast.Invert: operator.invert}
_CMP = {ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt, ast.GtE: operator.ge,
        ast.Eq: operator.eq, ast.NotEq: operator.ne}


class UnsafeExpression(ValueError):
    pass


def evaluate(expr, variables: dict):
    """Valeur de la formule `expr` (texte ou nombre) ; `variables` : {nom: valeur}."""
    if isinstance(expr, (int, float)) and not isinstance(expr, bool):
        return float(expr)
    text = str(expr)
    if len(text) > MAX_LENGTH:
        raise UnsafeExpression(f"Formule trop longue ({len(text)} caractères).")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as exc:
        raise UnsafeExpression(f"Formule invalide « {text} » : {exc.msg}") from None
    except (RecursionError, MemoryError):
        raise UnsafeExpression("Formule trop imbriquée.") from None
    names = {**CONSTANTS, **variables}

    def ev(n):
        if isinstance(n, ast.Constant):
            if isinstance(n.value, bool) or not isinstance(n.value, (int, float)):
                raise UnsafeExpression(f"« {text} » : seuls les nombres sont permis.")
            return float(n.value)
        if isinstance(n, ast.Name):
            if n.id in names:
                return names[n.id]
            raise UnsafeExpression(f"« {text} » : nom inconnu « {n.id} » (permis : "
                                   f"{', '.join(sorted(names))}, fonctions "
                                   f"{', '.join(sorted(FUNCTIONS))}).")
        if isinstance(n, ast.BinOp) and type(n.op) in _BIN:
            return _BIN[type(n.op)](ev(n.left), ev(n.right))
        if isinstance(n, ast.UnaryOp) and type(n.op) in _UNARY:
            return _UNARY[type(n.op)](ev(n.operand))
        if isinstance(n, ast.Compare) and all(type(o) in _CMP for o in n.ops):
            left, out = ev(n.left), None
            for op, comp in zip(n.ops, n.comparators):
                right = ev(comp)
                r = _CMP[type(op)](left, right)
                out = r if out is None else operator.and_(out, r)
                left = right
            return out
        if isinstance(n, ast.Call) and not n.keywords:
            f = n.func
            # « np.exp(x) » accepté (anciens cas) : seulement pour les fonctions permises
            if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) \
                    and f.value.id == "np":
                fname = f.attr
            elif isinstance(f, ast.Name):
                fname = f.id
            else:
                fname = None
            if fname in FUNCTIONS:
                return FUNCTIONS[fname](*[ev(a) for a in n.args])
            raise UnsafeExpression(f"« {text} » : fonction non permise ({ast.unparse(f)}). "
                                   f"Permises : {', '.join(sorted(FUNCTIONS))}.")
        raise UnsafeExpression(f"« {text} » : construction non permise "
                               f"({type(n).__name__}).")

    with np.errstate(all="ignore"):
        try:
            return ev(tree.body)
        except OverflowError:
            raise UnsafeExpression(f"« {text} » : dépassement de capacité.") from None
        except RecursionError:
            raise UnsafeExpression("Formule trop imbriquée.") from None
