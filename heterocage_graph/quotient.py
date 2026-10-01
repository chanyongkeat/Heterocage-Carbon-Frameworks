"""
Crystal quotient graph utilities.

A periodic structure is represented by a networkx MultiGraph whose nodes are
the atoms of one cell and whose edges (i, j) carry the integer cell
translation S ('vector') of the bond: atom j sits at x_j + S as seen from i.
Edges are always stored with i <= j and S oriented from i to j.
"""
import numpy as np
import networkx as nx


def get_cycle_sums(G):
    """Cycle-sum vectors of a quotient graph (one row per independent cycle)."""
    SG = nx.Graph(G)                       # simple graph; keeps one edge per pair
    sums = []
    for cycle in nx.cycle_basis(SG):
        s = np.zeros(3)
        for i, j in zip(cycle, cycle[1:] + cycle[:1]):
            s += SG[i][j]['vector'] if i <= j else -SG[i][j]['vector']
        sums.append(s)
    for i, j in SG.edges():                # parallel edges and periodic self-loops
        for e in G[i][j].values():
            sums.append(SG[i][j]['vector'] - e['vector'])
            if i == j:
                sums.append(e['vector'])
    return np.unique(np.array(sums).reshape(-1, 3), axis=0) if sums else np.zeros((0, 3))


def get_dimension(G):
    """Periodic dimensionality D in {0,1,2,3}: rank of the cycle-sum vectors (G2)."""
    sums = get_cycle_sums(G)
    return int(np.linalg.matrix_rank(sums)) if len(sums) else 0


def remove_selfloops(G):
    H = G.copy()
    H.remove_edges_from(list(nx.selfloop_edges(H)))
    return H
