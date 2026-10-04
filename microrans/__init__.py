"""microrans — solveur RANS / URANS 1D (canal plan) et 2D (volumes finis, mailleur intégré).

Modèles : Spalart-Allmaras, k-ε Launder-Sharma, k-ω Wilcox 2006, k-ω SST Menter 2003.
"""
import os as _os
import sys as _sys

# Un calcul = un fil BLAS. OpenBLAS (NumPy, SciPy) lance par défaut un fil par cœur sans
# rien gagner ici (cavité cubique 32³ : 20 s pour 73 s de CPU, contre 19 s pour 19 s avec
# un fil ; aucun cas mesuré plus lent avec un fil) et deux calculs simultanés
# s'effondrent (plaque compressible : 137 s au lieu de 9.7 s). À fixer avant le premier
# import de NumPy ; une valeur donnée par l'utilisateur est conservée.
if "numpy" not in _sys.modules:
    for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS"):
        _os.environ.setdefault(_k, "1")

from .grid import Grid, channel_grid, tanh_grid  # noqa: E402
from .models import MODELS, TURBULENT_MODELS, get_model  # noqa: E402
from .solver import Solution, UnsteadyResult, solve_steady, solve_unsteady  # noqa: E402

__version__ = "0.2.0"

__all__ = ["Grid", "channel_grid", "tanh_grid", "MODELS", "TURBULENT_MODELS", "get_model",
           "Solution", "UnsteadyResult", "solve_steady", "solve_unsteady", "__version__"]
