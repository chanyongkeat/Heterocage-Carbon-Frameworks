# heterocage-graph

Graph-theoretic identification, validation and classification of fullerene-cage frameworks:
2D networks built from covalently linked C<sub>20</sub>/C<sub>24</sub> cages (and, more generally, any
fullerene cages).

Given a relaxed structure (and, optionally, the unrelaxed input it came from), the code

1. builds the **crystal quotient graph** of the periodic bonding network,
2. **identifies the cages** and certifies them with a fullerene-graph test,
3. assigns each cage an **integrity status** (closed / open / broken / distorted),
4. **validates** the structure: atomic clashes, 2D covalent connectivity, **endohedral (cage-in-cage)**
   and **overlapping (interpenetrating)** cages,
5. **classifies** it as a *closed-cage network*, *partially opened* or *fully opened* framework,
6. computes the **inter-cage connection set** and the contracted cage graph, groups structures
   into **connection modes**, and checks that the graph is consistent with the **layer-group symmetry**.

The repository also contains the screening data set of 576 DFT-relaxed C<sub>24</sub>/C<sub>20</sub>
heterocage monolayers used to develop and validate the method (`dataset/`).

---

## Installation

```bash
git clone <repository-url> heterocage-graph
cd heterocage-graph
pip install -r requirements.txt
```

Requirements: Python ≥ 3.8, numpy, scipy, networkx ≥ 3.0, ase, spglib ≥ 2.5 (layer groups);
pandas is only needed for `scripts/run_dataset.py`. Tested with Python 3.10, numpy 2.1, scipy 1.14,
networkx 3.4, ase 3.24, spglib 2.6 and pandas 2.2.

No installation step is needed; run the commands below from the repository root.

## Quick start

```bash
# relaxed structure + unrelaxed input (recommended: cage identity is taken from the input)
python -m heterocage_graph examples/pareto/REF-48/CONTCAR --ref examples/pareto/REF-48/POSCAR

# a single structure (cages are identified on the structure itself)
python -m heterocage_graph examples/endohedral/POSCAR

# machine-readable output
python -m heterocage_graph examples/overlap/POSCAR --json
```

Output for the first command:

```
  label               : Closed-cage network
  cage identification : fullerene-certified
  cages (status)      : C24:closed, C20:closed, C20:closed
  closed / all        : 3 / 3
  |B| by cage pair    : 13  {'C20-C24': 10, 'C20-C20': 3}
  bonds per cage      : [10, 8, 8]   dim(G^) = 2
  dim(G), components  : 2, [(64, 2)]
  clashes (min d, A)  : 0 (None)   over-valent atoms: 0
  layer group         : p112 (No. 3, order 2); sigma(G)=G True; C(sG)=sC(G) True
  B symmetry orbits   : 8  sizes [2, 2, 2, 2, 2, 1, 1, 1]
  flags               : none
```

Python API:

```python
from heterocage_graph import classify, report, same_connection_mode

r = classify('examples/pareto/REF-48/CONTCAR', ref='examples/pareto/REF-48/POSCAR')
print(r['label'], r['cage_status'], r['B'])
print(report(r))
```

Any structure format readable by ASE can be used (POSCAR/CONTCAR, cif, xyz with cell, ...).
With `--ref`, the two files must contain the same atoms **in the same order** (true for a VASP
relaxation: CONTCAR keeps the POSCAR order).

### Why the unrelaxed input?

To decide whether a cage *survived* relaxation, one must know which atoms formed that cage.
In the unrelaxed input every cage is intact, so the partition into cages is unambiguous and can be
certified. A cage that opens during relaxation is no longer a fullerene graph and can no longer be
recognised from the relaxed structure alone (e.g. `R2B-128_002`: 5 of 6 cages open; no valid
partition exists in the relaxed structure). The input partition is therefore imposed on the
relaxed structure, and the status of every cage is judged against its own input graph.

---

## Method

All decisions are explicit predicates on graphs; tags (G1)–(G14) are also used in the source code.

**(G1) Bond graph.** For a cell with lattice matrix $\mathsf L$, atoms $n_i$ at $\vec x_i$ and covalent radii
$r_i$, the bonding network is the labelled quotient multigraph $G=(V,E)$,

$$E=\{(n_i,n_j,\vec S): \lambda(r_i+r_j)\le d_{ij}(\vec S)<k(r_i+r_j)\},\qquad
\mathcal K=\{(n_i,n_j,\vec S): d_{ij}(\vec S)<\lambda(r_i+r_j)\},$$

with $d_{ij}(\vec S)=\lVert\vec x_j+\mathsf L^{\mathsf T}\vec S-\vec x_i\rVert$ and cell translation $\vec S\in\mathbb Z^3$.
Defaults: $k=1.2$ (1.824 Å for C–C) and $\lambda=0.75$ (1.14 Å, below the shortest C–C bond, 1.20 Å).
Contacts in $\mathcal K$ are clashes, not bonds.

**(G2) Dimensionality.** $D(X)=\operatorname{rank}\{\vec s(c)=\sum_{e\in c}\pm\vec S_e\}$ over a cycle basis of $X$
(plus parallel-edge differences). $D=0$: finite molecule, $D=2$: layer.

**(G3) Cage predicate.** A graph $H$ is a closed cage (fullerene graph) iff

$$\Phi(H)=1 \iff H \text{ connected} \wedge \deg_H(v)=3\ \forall v \wedge H\text{ planar} \wedge |f|\in\{5,6\}\ \forall f\in F(H).$$

Euler's relation then gives exactly 12 pentagons and $|V|/2-10$ hexagons. For C<sub>20</sub> ($I_h$) and
C<sub>24</sub> ($D_{6d}$) the fullerene graph is unique, so $\Phi$ coincides with VF2 isomorphism to the
reference molecules in `heterocage_graph/data/`.

**Intra-cage bonds.** Each cage is unwrapped by placing every atom at its image closest to the cage
centre (offsets $\vec o_i$). A bond $(n_i,n_j,\vec S)$ between atoms of the same cage is intra-cage iff
$\vec S=\vec o_j-\vec o_i$; otherwise it binds the cage to **its own periodic image** and is an inter-cage bond.
The internal graph $H_a$ of cage $a$ contains only intra-cage bonds.

**(G4) Decomposition.** Connected components are split recursively by Girvan–Newman edge
betweenness until every subunit satisfies $\Phi(H_a)=1$:
$\mathcal C(X)=\{X\}$ if $\Phi(H_X)=1$, otherwise $\bigcup_{Y\in\mathrm{GN}(X)}\mathcal C(Y)$.
Several bond windows propose partitions; a partition is accepted only if it is **certified**:
(1) all subunits satisfy $\Phi$ (`fullerene-certified`); (2) all subunits are connected 20/24-atom cages
(`size-certified`, slightly distorted inputs); (3) exact cover by induced subgraphs isomorphic to
C<sub>20</sub>/C<sub>24</sub> (`exact-cover`). Otherwise the input is reported as `Unresolved input`.

**(G5) Inter-cage connection set and contracted graph.**

$$\mathcal B=\{(n_i,n_j,\vec S)\in E: n_i\in X_a,\ n_j\in X_b,\ \text{not intra-cage}\},\qquad
\hat G:\ (a,b,\hat{\vec S}=\vec S+\vec o_i-\vec o_j).$$

$\deg_{\hat G}(a)$ is the number of inter-cage bonds of cage $a$; $D(\hat G)$ is the dimensionality of the cage net.

**(G6) Connection mode.** $\hat G_1\simeq\hat G_2$ iff there is a bijection $\varphi$ of cages preserving cage type and
edge multiplicities and integer gauges $\vec c_a$ with
$\{\{\hat{\vec S}\}\}_{\varphi(a)\varphi(b)}=\{\{\hat{\vec S}+\vec c_a-\vec c_b\}\}_{ab}$ for all $a,b$
(VF2 + gauge search; `same_connection_mode`).

**(G7) Cage status** of cage $a$ with relaxed internal graph $H_a$ and input graph $g_a$:
*closed* if $H_a\cong g_a$; *broken* if $H_a$ is disconnected; *open* if its largest face has more than
6 atoms (a window); *distorted* otherwise (rearranged or extra intra-cage bonds).

**(G8)–(G11) Validation.**

| predicate | condition | meaning |
|---|---|---|
| $V_{\rm clash}$ | $\mathcal K=\varnothing$ | no unphysically short contacts |
| $V_{\rm val}$ | $\deg_G(n)\le 4$ for C | valence check (diagnostic flag only) |
| $V_{\rm 2D}$ | $D(G)=2$ and every cage in a 2D component | covalent 2D network, no van der Waals-packed cage |
| $V_{\rm endo}$ | no $X_b+\vec T\subset\operatorname{conv}(X_a)$ | no cage nested inside another (centroid distance < 1 Å also flagged) |
| $V_{\rm ov}$ | $\rho_{ab}(\vec T)\le\varepsilon$ for all non-nested pairs | no interpenetrating cages |

$\rho_{ab}$ is the radius of the largest ball inside $\operatorname{conv}(X_a)\cap\operatorname{conv}(X_b+\vec T)$ (Chebyshev
centre, one linear program over the facet inequalities of both convex hulls); $\varepsilon=0.5$ Å.
In the data set the hulls of bonded cages never intersect ($\rho=0$ for all 557 resolved structures), while the interpenetrating example gives $\rho=1.93$ Å.

**(G12) Classification** ($n_c$ = number of closed cages, $p$ = number of cages):

| label | condition |
|---|---|
| Endohedral cage-in-cage (rejected) | $\neg V_{\rm endo}$ |
| Overlapping cages (rejected) | $\neg V_{\rm ov}$ or $\neg V_{\rm clash}$ |
| Non-2D / van der Waals-packed (rejected) | $\neg V_{\rm 2D}$ |
| **Closed-cage network** | $n_c=p$ |
| **Partially opened** | $0<n_c<p$ |
| **Fully opened** | $n_c=0$ |
| Unresolved input | no certified cage partition of the input |

**(G13)–(G14) Symmetry.** The layer group $\mathcal L$ (spglib, `symprec` = 0.1 Å) acts on the labelled graph by
$\sigma\cdot(n_i,n_j,\vec S)=(n_{\pi_\sigma(i)},n_{\pi_\sigma(j)},\mathsf R_\sigma\vec S+\vec t_\sigma(j)-\vec t_\sigma(i))$.
The code verifies $\sigma(G)=G$ and $\mathcal C(\sigma G)=\sigma\,\mathcal C(G)$ for every operation and reports the
symmetry orbits of $\mathcal B$; by Burnside's lemma
$|\mathcal B/\bar{\mathcal L}|=\frac{1}{|\bar{\mathcal L}|}\sum_{\sigma}|\mathrm{Fix}_{\mathcal B}(\sigma)|$
(REF-48: $(13+3)/2=8$).

### Parameters

| parameter | default | where |
|---|---|---|
| bond factor $k$ | 1.2 | `--coef`, `COEF` |
| clash factor $\lambda$ | 0.75 | `CLASH` |
| endohedral centroid tolerance | 1.0 Å | `CENTER_TOL` |
| overlap depth $\varepsilon$ | 0.5 Å | `HULL_EPS` |
| symmetry tolerance | 0.1 Å | `--symprec`, `SYMPREC` |

The C–C distance distribution of all relaxed structures has its minimum at 1.90–1.95 Å; with
$k=1.1$ (1.672 Å) the elongated (≈1.7 Å) bonds of intact, sp³-bridged C<sub>20</sub> cages are lost and
closed cages are misread as open, whereas results are stable for $1.2\le k\le1.3$.

---

## Examples (`examples/`)

| example | input | expected label | key evidence |
|---|---|---|---|
| `pareto/REF-48` | CONTCAR + POSCAR | Closed-cage network | 3/3 closed, $\lvert\mathcal B\rvert=13$, layer group *p*112, 8 bond orbits |
| `pareto/REF-37` | CONTCAR + POSCAR | Closed-cage network | 3/3 closed, $\lvert\mathcal B\rvert=12$ |
| `pareto/R1S-68_032` | CONTCAR + POSCAR | Closed-cage network | 3/3 closed, $\lvert\mathcal B\rvert=16$ |
| `pareto/R1S-68_024` | CONTCAR + POSCAR | Partially opened | 1/3 closed (two C<sub>24</sub> distorted) |
| `pareto/R2B-128_002` | CONTCAR + POSCAR | Partially opened | 1/6 closed (five cages open, windows of 8–16 atoms) |
| `pareto/R2B-88_031` | CONTCAR + POSCAR | Fully opened | 0/4 closed, windows of 8–14 atoms |
| `endohedral/POSCAR` | single structure | Endohedral cage-in-cage | C<sub>20</sub>@C<sub>80</sub>, centroid distance 0.00 Å |
| `overlap/POSCAR` | single structure | Overlapping cages | two C<sub>80</sub>/C<sub>60</sub> pairs, 60 clashes ($d_{\min}=0.79$ Å), $\rho=1.93$ Å |

The six `pareto/` structures form the rank-0 Pareto front (lowest energy per atom vs. largest PBE band
gap) of the data set; `examples/pareto/pareto_front.txt` lists their energies and gaps.
Run all examples and regenerate `examples/expected_output.txt`:

```bash
python scripts/run_examples.py
```

## Data set (`dataset/`)

576 converged, DFT-relaxed C<sub>24</sub>/C<sub>20</sub> heterocage monolayers from a high-throughput
search, 44–208 atoms per cell, twelve compositions (C<sub>24</sub>)<sub>a</sub>(C<sub>20</sub>)<sub>b</sub>.

```
dataset/
├── dataset_index.csv          id, origin, composition, E0, E0/atom, PBE gap
├── graph_classification.csv   dataset_index.csv + all graph results (scripts/run_dataset.py)
└── structures/<id>/
    ├── POSCAR                 unrelaxed input (fixes cage identity)
    ├── CONTCAR                relaxed structure
    ├── energy.txt             total energy E0 of the relaxed cell (eV)
    └── gap.txt                PBE band gap (eV); a negative value means no gap (metallic)
```

Structure ids encode the search run: `REF` = reference run, `R1S`/`R1B` = run 1 small/big cells,
`R2S`/`R2B` = run 2 small/big cells, followed by the original folder name (`dataset_index.csv`, column
`source_dir`). Folder names in the big-cell runs are not atom counts (e.g. `R1B-148_*` contain 128-atom
cells); use the `natoms` column.

| composition | atoms | structures |
|---|---|---|
| (C<sub>24</sub>)<sub>1</sub>(C<sub>20</sub>)<sub>1</sub> | 44 | 41 |
| (C<sub>24</sub>)<sub>1</sub>(C<sub>20</sub>)<sub>2</sub> | 64 | 75 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>1</sub> | 68 | 69 |
| (C<sub>24</sub>)<sub>1</sub>(C<sub>20</sub>)<sub>3</sub> | 84 | 79 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>2</sub> | 88 | 63 |
| (C<sub>24</sub>)<sub>1</sub>(C<sub>20</sub>)<sub>4</sub> | 104 | 47 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>3</sub> | 108 | 52 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>4</sub> | 128 | 111 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>5</sub> | 148 | 6 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>6</sub> | 168 | 18 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>7</sub> | 188 | 13 |
| (C<sub>24</sub>)<sub>2</sub>(C<sub>20</sub>)<sub>8</sub> | 208 | 2 |

Classification of the data set (`python scripts/run_dataset.py`, about 30 min on 4 cores):

| label | structures |
|---|---|
| Closed-cage network | 283 |
| Partially opened | 231 |
| Fully opened | 6 |
| Non-2D / van der Waals-packed (rejected) | 37 |
| Overlapping cages (rejected) | 0 |
| Endohedral cage-in-cage (rejected) | 0 |
| Unresolved input | 19 |

Connection modes (G6) among the 557 structures with resolved cages: 460 distinct by cage type and
bond multiplicity, 538 including the translation labels (up to gauge).

`graph_classification.csv` columns: `label`, `identification`, `n_cages`, `n_closed`, `n_open`, `n_broken`,
`n_distorted`, `B` and its split `B_C20_C20`/`B_C20_C24`/`B_C24_C24`, `dim_G`, `dim_Ghat`, `clashes`,
`max_hull_depth_A`, `over_valent_atoms`, `endohedral`, `overlap`, `layer_group`, `layer_group_no`,
`B_symmetry_orbits`, `connection_mode_topology`, `connection_mode_labelled`, `flags`.

---

## Repository layout

```
heterocage_graph/
├── classifier.py     bond graph, cage identification, status, validation, classification, symmetry
├── quotient.py       quotient-graph dimensionality (cycle-sum rank)
├── __main__.py       command line interface
└── data/             C20.xyz, C24.xyz reference molecules
scripts/
├── run_examples.py   classify examples/ and write examples/expected_output.txt
└── run_dataset.py    classify dataset/ and write dataset/graph_classification.csv
tests/
└── test_examples.py  regression tests (pytest, or python tests/test_examples.py)
examples/             Pareto-front structures, endohedral and overlap test cases
dataset/              the full data set
```

## Tests

```bash
python -m pytest tests -q        # or: python tests/test_examples.py
```

## Limitations

- Cage identity relies on a certified partition of the input. Inputs whose cages are already fused
  are reported as `Unresolved input` rather than guessed.
- Unwrapping requires each cage radius to be below half the shortest lattice vector; otherwise a
  warning flag is raised and the status of that cage should be checked by hand (4 structures of the
  data set: `R2S-64_001`, `R2S-68_001`, `R2B-88_001`, `R2B-88_026`, all with a shortest lattice vector of 5.3–6.2 Å).
- The size-certified and exact-cover fallbacks use the C<sub>20</sub>/C<sub>24</sub> reference library; other cage
  sizes are identified through the fullerene predicate (G3) only.
- Hull-based endohedral/overlap tests assume convex, cage-like subunits.

## Citation

If you use this code or data, please cite: *<citation to be added on publication>*.

## License

*<license to be added>*
