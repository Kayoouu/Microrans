"""microrans — micro-solveur RANS / URANS 1D (canal plan turbulent).

Modèles : Spalart-Allmaras, k-ε Launder-Sharma, k-ω Wilcox 2006, k-ω SST Menter 2003.
"""
from .grid import Grid, channel_grid, tanh_grid
from .models import MODELS, TURBULENT_MODELS, get_model
from .solver import Solution, UnsteadyResult, solve_steady, solve_unsteady

__version__ = "0.1.0"

__all__ = ["Grid", "channel_grid", "tanh_grid", "MODELS", "TURBULENT_MODELS", "get_model",
           "Solution", "UnsteadyResult", "solve_steady", "solve_unsteady", "__version__"]
