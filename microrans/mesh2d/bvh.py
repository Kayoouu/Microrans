"""Arbre de boîtes englobantes (alignées sur les axes) des faces de paroi, pour la distance
exacte à la paroi en 2D (segments) et en 3D (polygones).

La distance d'un point à la boîte d'un groupe de faces est une borne inférieure exacte de
sa distance à chacune de ces faces : une boîte plus loin que la meilleure distance déjà
trouvée est écartée avec toutes ses faces. Contrairement à une recherche par centres de
faces (sphère de rayon d + plus grand rayon de face), la boîte suit l'orientation et la
taille des faces qu'elle contient : une paroi plane vue de loin, ou une paroi mêlant petites
et grandes faces, ne fait examiner que quelques faces par cellule.
"""
from __future__ import annotations

import numpy as np


def build(lo_f, hi_f, leaf: int = 8):
    """Arbre binaire : coupe à la médiane des centres de boîtes selon le plus grand côté,
    feuilles de ≤ `leaf` faces. lo_f, hi_f : (nf, dim) coins des boîtes des faces.
    Renvoie (order, start, end, left, right, lo, hi) : les faces du nœud i sont
    order[start[i]:end[i]] ; left / right = −1 pour une feuille."""
    nf, dim = lo_f.shape
    cen = 0.5 * (lo_f + hi_f)
    order = np.arange(nf)
    start, end, left, right = [0], [nf], [-1], [-1]
    stack = [0]
    while stack:
        i = stack.pop()
        s, e = start[i], end[i]
        if e - s <= leaf:
            continue
        seg = order[s:e]
        cs = cen[seg]
        ax = int(np.argmax(cs.max(axis=0) - cs.min(axis=0)))
        mid = (e - s) // 2
        order[s:e] = seg[np.argpartition(cs[:, ax], mid)]
        left[i], right[i] = len(start), len(start) + 1
        start += [s, s + mid]
        end += [s + mid, e]
        left += [-1, -1]
        right += [-1, -1]
        stack += [left[i], right[i]]
    start, end = np.array(start), np.array(end)
    left, right = np.array(left), np.array(right)
    lo, hi = np.empty((len(start), dim)), np.empty((len(start), dim))
    leaves = np.nonzero(left < 0)[0]
    leaves = leaves[np.argsort(start[leaves])]               # segments contigus, dans l'ordre
    lo[leaves] = np.minimum.reduceat(lo_f[order], start[leaves], axis=0)
    hi[leaves] = np.maximum.reduceat(hi_f[order], start[leaves], axis=0)
    for i in range(len(start) - 1, -1, -1):                  # enfants créés après le parent
        if left[i] >= 0:
            lo[i] = np.minimum(lo[left[i]], lo[right[i]])
            hi[i] = np.maximum(hi[left[i]], hi[right[i]])
    return order, start, end, left, right, lo, hi


def search(P, cells, best, tree, visit) -> None:
    """Descente simultanée de l'arbre pour les points P[cells] : une boîte est ouverte si sa
    distance au point (marge d'arrondi relative 1e-9, du côté prudent) est sous best[cellule] ;
    visit(c, f) reçoit les paires (cellule, face) des feuilles atteintes et met best à jour
    en place (la descente suivante en profite)."""
    order, start, end, left, right, lo, hi = tree
    node = np.zeros(len(cells), dtype=np.int64)
    while len(cells):
        q = P[cells]
        gap = np.maximum(np.maximum(lo[node] - q, q - hi[node]), 0.0)
        keep = np.sqrt(np.sum(gap * gap, axis=1)) * (1.0 - 1e-9) < best[cells]
        cells, node = cells[keep], node[keep]
        leaf = left[node] < 0
        if np.any(leaf):
            lc, ln = cells[leaf], node[leaf]
            cnt = end[ln] - start[ln]
            pos = np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt)
            visit(np.repeat(lc, cnt), order[np.repeat(start[ln], cnt) + pos])
        cells, node = cells[~leaf], node[~leaf]
        cells = np.concatenate([cells, cells])
        node = np.concatenate([left[node], right[node]])
