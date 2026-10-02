# GraphTheory — heterocage structure validation and classification

Graph-theory tool that checks whether a relaxed fullerene-cage framework (e.g. C<sub>24</sub>/C<sub>20</sub>
heterocage monolayers) is physically valid and classifies it. For each structure it:

- identifies the cages and checks whether each one stays closed, opens, breaks or distorts after relaxation,
- rejects invalid structures: atomic clashes, overlapping cages, cage-in-cage (endohedral), non-2D networks,
- classifies valid ones as **Closed-cage network**, **Partially opened** or **Fully opened**,
- reports the inter-cage bonds, connection mode and layer-group symmetry.

## Download

```bash
git clone https://github.com/chanyongkeat/Heterocage-Carbon-Frameworks.git
cd Heterocage-Carbon-Frameworks/GraphTheory
pip install -r requirements.txt   # Python >= 3.8, numpy, scipy, networkx, ase, spglib, pandas
```

## Usage

Run from the `GraphTheory/` folder. Any ASE-readable format works (POSCAR/CONTCAR, cif, xyz, ...).

```bash
# relaxed structure + unrelaxed input (recommended; same atom order)
python -m heterocage_graph examples/pareto/REF-48/CONTCAR --ref examples/pareto/REF-48/POSCAR

# single structure
python -m heterocage_graph examples/endohedral/POSCAR

# JSON output
python -m heterocage_graph examples/overlap/POSCAR --json
```

Python API:

```python
from heterocage_graph import classify, report
r = classify('CONTCAR', ref='POSCAR')
print(report(r))
```

Batch runs:

```bash
python scripts/run_examples.py   # classify examples/
python scripts/run_dataset.py    # classify dataset/ -> dataset/graph_classification.csv
```

## Test

```bash
python -m pytest tests -q        # or: python tests/test_examples.py
```

## Contents

| path | function |
|---|---|
| `heterocage_graph/` | the package: bond graph, cage identification, validation, classification, symmetry (`classifier.py`, `quotient.py`, CLI in `__main__.py`) |
| `scripts/` | batch runs over `examples/` and `dataset/` |
| `tests/` | regression tests |
| `examples/` | sample structures (closed, partially/fully opened, endohedral, overlap) |
| `dataset/` | DFT-relaxed C<sub>24</sub>/C<sub>20</sub> structures with energies, band gaps and classification results |
