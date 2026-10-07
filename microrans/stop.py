"""Arrêt demandé (bouton « Arrêter » de l'interface) pendant une étape longue qui n'a pas de
point d'arrêt à elle (maillage non structuré : avant, mené à son terme, U20). La tâche
installe une fonction de test avec `stop_check` ; les boucles longues appellent
`check_stop()`, qui lève StopRequested si l'arrêt est demandé."""
from __future__ import annotations

import contextvars
from contextlib import contextmanager


class StopRequested(Exception):
    """Arrêt demandé par l'utilisateur."""


_CHECK = contextvars.ContextVar("microrans_stop_check", default=None)


@contextmanager
def stop_check(fn):
    """Dans ce bloc (même fil d'exécution), check_stop() appelle fn() : vrai → arrêt."""
    token = _CHECK.set(fn)
    try:
        yield
    finally:
        _CHECK.reset(token)


def check_stop(what: str = "le maillage"):
    fn = _CHECK.get()
    if fn is not None and fn():
        raise StopRequested(f"Arrêt demandé pendant {what}.")
