"""Sécurité : un fichier de cas venu d'ailleurs ne doit pas pouvoir exécuter de code."""
import numpy as np
import pytest

from microrans.fv2d.case import run_case
from microrans.safe_expr import UnsafeExpression, evaluate

X = {"x": np.linspace(0, 1, 5), "y": np.linspace(0, 2, 5)}


@pytest.mark.parametrize("expr, expected", [
    ("6*y*(1-y)", 6 * X["y"] * (1 - X["y"])),
    ("np.exp(-x) + sqrt(y)", np.exp(-X["x"]) + np.sqrt(X["y"])),
    ("where((x > 0.2) & (y < 1.5), 1, 0)", np.array([0, 1, 1, 0, 0])),
    ("0 < x < 0.6", np.array([0, 1, 1, 0, 0])),
    ("2*pi*e", 2 * np.pi * np.e), (3, 3.0)])
def test_allowed_formulas(expr, expected):
    assert np.allclose(evaluate(expr, X), expected)


@pytest.mark.parametrize("expr", [
    "().__class__.__base__.__subclasses__()",
    "__import__('os').system('echo pwned')",
    "open('/etc/passwd').read()", "x.__class__", "(lambda: 1)()", "[1, 2][0]",
    "'a' * 3", "exp(x, out=x)", "np.load('f.npy')", "True", "9**9**9",
    "-" * 1999 + "1", "x" * 3000, "x if y else 0", "{1: 2}", "f'{x}'"])
def test_forbidden_constructs(expr):
    with pytest.raises(UnsafeExpression):
        evaluate(expr, X)


def test_malicious_case_file_is_refused_without_side_effect(tmp_path):
    marker = tmp_path / "pwned"
    payload = (f"().__class__.__base__.__subclasses__()[0]"
               f" or open('{marker.as_posix()}', 'w')")
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0, "y1": 1, "nx": 4, "ny": 4,
                    "names": {"left": "inlet", "right": "outlet", "bottom": "w", "top": "w"}},
           "physics": {"nu": 0.1},
           "boundary": {"inlet": {"type": "inlet", "U": [payload, 0.0]},
                        "outlet": {"type": "outlet"}, "w": {"type": "wall"}},
           "output": {"vtk": False, "plots": False}}
    with pytest.raises(UnsafeExpression):
        run_case(cfg, out_dir=tmp_path / "out", verbose=False, plot=False)
    assert not marker.exists()


def test_line_names_cannot_escape_output_directory(tmp_path):
    cfg = {"mesh": {"type": "rectangle", "x0": 0, "x1": 1, "y0": 0, "y1": 1, "nx": 4, "ny": 4,
                    "names": {"left": "w", "right": "w", "bottom": "w", "top": "lid"}},
           "physics": {"nu": 0.1},
           "boundary": {"lid": {"type": "wall", "U": [1.0, 0.0]}, "w": {"type": "wall"}},
           "solver": {"max_iter": 5},
           "output": {"vtk": False, "plots": False,
                      "lines": [{"name": "../../evil", "start": [0, 0], "end": [1, 1]}]}}
    with pytest.raises(ValueError, match="nom de fichier"):
        run_case(cfg, out_dir=tmp_path / "out", verbose=False, plot=False)
