"""Matrices creuses CSR écrites uniquement avec des opérations de tableaux (lecture indexée,
`bincount`), pour les bibliothèques « type NumPy » sans module creux (dpnp d'Intel). Même
interface que le sous-ensemble de scipy.sparse / cupyx.scipy.sparse utilisé par
microrans : construction (data, indices, indptr), produit matrice-vecteur, diagonale,
dense. Les indices de colonnes doivent être triés par ligne (c'est le cas des matrices
assemblées par fvm.py).
"""
from __future__ import annotations


class GenericCSR:
    def __init__(self, arg, shape=None, xp=None):
        if isinstance(arg, GenericCSR):
            self.data, self.indices, self.indptr = arg.data, arg.indices, arg.indptr
            self.shape, self.xp = arg.shape, arg.xp
            self._rows = arg._rows
            return
        data, indices, indptr = arg
        self.xp = xp
        self.data, self.indices, self.indptr = data, indices, indptr
        self.shape = tuple(shape)
        self._rows = None

    @property
    def nnz(self):
        return int(self.data.shape[0])

    has_sorted_indices = True

    def sort_indices(self):
        pass

    def tocsr(self):
        return self

    @property
    def rows(self):
        if self._rows is None:
            xp = self.xp
            self._rows = xp.repeat(xp.arange(self.shape[0]), xp.diff(self.indptr))
        return self._rows

    def __matmul__(self, x):
        xp = self.xp
        return xp.bincount(self.rows, weights=self.data * x[self.indices],
                           minlength=self.shape[0])

    def diagonal(self):
        xp = self.xp
        on = self.rows == self.indices
        return xp.bincount(self.rows, weights=xp.where(on, self.data, 0.0),
                           minlength=self.shape[0])

    def to_scipy(self, host):
        import scipy.sparse as sp
        return sp.csr_matrix((host(self.data), host(self.indices), host(self.indptr)),
                             shape=self.shape)

    def toarray(self):
        from .backend import to_host_array
        return self.xp.asarray(self.to_scipy(to_host_array).toarray())


class GenericSparseModule:
    """Module creux minimal lié à un module de tableaux `xp`."""

    def __init__(self, xp):
        self.xp = xp

    def csr_matrix(self, arg, shape=None):
        if isinstance(arg, GenericCSR):
            return GenericCSR(arg)
        return GenericCSR(arg, shape, self.xp)
