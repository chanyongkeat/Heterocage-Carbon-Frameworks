"""heterocage_graph: graph-theoretic identification, validation and classification
of fullerene-cage frameworks (closed / partially opened / fully opened networks,
endohedral and overlapping cages)."""
from .classifier import (classify, report, identify_cages, bond_graph, is_fullerene,
                         same_connection_mode, COEF, CLASH, SYMPREC)
from .quotient import get_dimension

__all__ = ['classify', 'report', 'identify_cages', 'bond_graph', 'is_fullerene',
           'same_connection_mode', 'get_dimension', 'COEF', 'CLASH', 'SYMPREC']
__version__ = '1.0.0'
