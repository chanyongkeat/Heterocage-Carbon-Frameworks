"""Regression tests:  python -m pytest tests -q"""
import os
import sys

import numpy as np
try:
    import pytest
except ImportError:                  # plain `python tests/test_examples.py` also works
    pytest = None
from ase.io import read, write

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from heterocage_graph import classify, same_connection_mode  # noqa: E402

EX = os.path.join(ROOT, 'examples')
PARETO = {  # folder: (label, closed cages, cages, |B|)
    'REF-48': ('Closed-cage network', 3, 3, 13),
    'REF-37': ('Closed-cage network', 3, 3, 12),
    'R1S-68_032': ('Closed-cage network', 3, 3, 16),
    'R1S-68_024': ('Partially opened', 1, 3, 18),
    'R2B-128_002': ('Partially opened', 1, 6, 29),
    'R2B-88_031': ('Fully opened', 0, 4, 20),
}


def pareto(name, **kw):
    d = os.path.join(EX, 'pareto', name)
    return classify(os.path.join(d, 'CONTCAR'), os.path.join(d, 'POSCAR'), **kw)


def _parametrize(names):
    return pytest.mark.parametrize('name', names) if pytest else (lambda f: f)


@_parametrize(sorted(PARETO))
def test_pareto_labels(name):
    label, n_closed, n_cages, nB = PARETO[name]
    r = pareto(name, symmetry=False)
    assert (r['label'], r['n_closed'], r['n_cages'], r['B']) == (label, n_closed, n_cages, nB)


def test_ref48_symmetry_orbits():
    s = pareto('REF-48')['symmetry']
    assert s['layer_group'] == 'p112' and s['graph_invariant'] and s['partition_equivariant']
    assert s['B_orbits'] == 8          # Burnside: (13 + 3) / 2


def test_endohedral():
    r = classify(os.path.join(EX, 'endohedral', 'POSCAR'), symmetry=False)
    assert r['label'].startswith('Endohedral')
    assert r['endohedral'][0].startswith('C20@C80')


def test_overlap():
    r = classify(os.path.join(EX, 'overlap', 'POSCAR'), symmetry=False)
    assert r['label'].startswith('Overlapping')
    assert sorted(r['cages']) == ['C60', 'C60', 'C80', 'C80'] and len(r['overlap']) == 2


def test_edge_orientation_large_cell():
    """Regression: subgraph views may iterate bonds as (j, i); S must not be flipped."""
    d = os.path.join(ROOT, 'dataset', 'structures', 'REF-245')
    r = classify(os.path.join(d, 'CONTCAR'), os.path.join(d, 'POSCAR'), symmetry=False)
    assert r['label'] == 'Closed-cage network' and r['n_closed'] == 8


def test_connection_mode_invariance(tmp_path):
    """Renumbering atoms and translating one cage by a lattice vector keeps the mode."""
    d = os.path.join(EX, 'pareto', 'REF-48')
    r0 = pareto('REF-48', symmetry=False)
    perm = np.random.default_rng(0).permutation(r0['natoms'])
    moved = r0['cage_atoms'][1]
    for fn in ('POSCAR', 'CONTCAR'):
        a = read(os.path.join(d, fn))
        pos = a.get_positions()
        pos[moved] += a.cell[0]                     # rigid lattice translation of one cage
        a.set_positions(pos)
        write(tmp_path / fn, a[perm], format='vasp')
    r1 = classify(str(tmp_path / 'CONTCAR'), str(tmp_path / 'POSCAR'), symmetry=False)
    assert same_connection_mode(r0, r1) and same_connection_mode(r0, r1, use_labels=False)
    assert not same_connection_mode(r0, pareto('REF-37', symmetry=False), use_labels=False)


if __name__ == '__main__':          # also runnable without pytest
    import tempfile
    import pathlib
    for n in sorted(PARETO):
        test_pareto_labels(n)
        print('ok  pareto', n)
    for t in (test_ref48_symmetry_orbits, test_endohedral, test_overlap, test_edge_orientation_large_cell):
        t()
        print('ok ', t.__name__)
    with tempfile.TemporaryDirectory() as tmp:
        test_connection_mode_invariance(pathlib.Path(tmp))
    print('ok  test_connection_mode_invariance')
