"""
Graph-theoretic identification, validation and classification of cage frameworks.

Every decision is an explicit predicate on the crystal quotient graph
G = (V, E, S); the tags (Gn) are the equation labels used in README.md.

  (G1)  bond graph       lam (r_i + r_j) <= d_ij(S) < k (r_i + r_j)  -> bond (n_i, n_j, S)
                         d_ij(S) < lam (r_i + r_j)                   -> clash set K
  (G2)  dimensionality   D(X) = rank of the cycle-sum vectors        (quotient.get_dimension)
  (G3)  cage predicate   Phi(H): connected, 3-regular, planar, faces in {5, 6} (=> 12 pentagons)
  (G4)  decomposition    Girvan-Newman recursion until every subunit satisfies Phi
  (G5)  B and G^         inter-cage connection set and contracted quotient graph
  (G6)  connection mode  isomorphism of G^ preserving cage type, multiplicity and S up to gauge
  (G7)  cage status      closed / broken / open / distorted
  (G8)  V_clash, V_val   no clash; deg(n) <= 4 for carbon (V_val is a diagnostic only)
  (G9)  V_2D             D(G) = 2 and every cage lies in a 2D component
  (G10) V_endo           no cage nested inside the convex hull of another (or its image)
  (G11) V_ov             no hull interpenetration (Chebyshev radius of the intersection <= eps)
  (G12) classification   endohedral / overlapping / non-2D / CC / PO / FO; accepted iff CC
  (G13) group action     sigma (n_i, n_j, S) = (n_pi(i), n_pi(j), R S + t(j) - t(i)); sigma(G) = G
  (G14) equivariance     C(sigma G) = sigma C(G); orbits of B (Burnside)
"""
import os
import itertools
from collections import Counter

import numpy as np
import networkx as nx
from ase.io import read
from ase.data import covalent_radii
from ase.neighborlist import neighbor_list
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, Delaunay
from networkx.algorithms.isomorphism import GraphMatcher

from .quotient import get_dimension, remove_selfloops

COEF = 1.2          # k:   bond if d < k (r_i + r_j)       -> 1.824 A for C-C
CLASH = 0.75        # lam: clash if d < lam (r_i + r_j)   -> 1.140 A (< shortest C-C bond, 1.20 A)
MAX_VALENCE = {6: 4}
CENTER_TOL = 1.0    # A, coincident cage centroids
HULL_EPS = 0.50     # A, inscribed radius of a hull intersection counted as overlap (bonded cages <= 0.35 A)
SYMPREC = 0.1       # A, spglib tolerance
# bond windows (lower, upper) x (r_i + r_j) that PROPOSE partitions; (G3) certifies them
LADDER = [(CLASH, c) for c in (1.05, 1.10, 1.15, 1.20, 1.00)] + [(0.855, c) for c in (1.00, 1.05)]
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
LIBRARY = {n: os.path.join(DATA, f'C{n}.xyz') for n in (20, 24)}    # reference molecules


def load(path):
    atoms = read(path)
    atoms.wrap()
    return atoms


# ======================================================================= (G1)
def bond_graph(atoms, coef=COEF, lam=CLASH):
    """Quotient multigraph; edges stored with i <= j and S oriented i -> j.

    Returns (G, clashes); clashes are (i, j, S, d) contacts shorter than lam (r_i + r_j).
    """
    r = np.array([covalent_radii[z] for z in atoms.get_atomic_numbers()])
    i, j, S, d = neighbor_list('ijSd', atoms, list(r * coef), max_nbins=10)
    G = nx.MultiGraph()
    G.add_nodes_from(range(len(atoms)))
    clashes = []
    for a, b, s, dd in zip(i, j, S, d):
        if a > b or (a == b and tuple(s) < (0, 0, 0)):     # keep each bond once
            continue
        if dd < lam * (r[a] + r[b]):
            clashes.append((int(a), int(b), tuple(int(x) for x in s), float(dd)))
            continue
        G.add_edge(int(a), int(b), vector=np.asarray(s, int))
    return G, clashes


def oriented_edges(G, nodes=None):
    """Yield (i, j, S) with i <= j and S oriented i -> j.

    S is stored oriented from the smaller to the larger index (bond_graph), but a
    (sub)graph view may iterate an edge as (j, i): only the endpoints are reordered.
    """
    for u, v, dat in (G.subgraph(nodes) if nodes is not None else G).edges(data=True):
        yield min(u, v), max(u, v), np.asarray(dat['vector'], int)


# ============================================== cage unwrapping and internal graph
def _lattice_shifts(atoms, reach=2):
    """Integer translations over the periodic directions, |t_k| <= reach."""
    rng = [range(-reach, reach + 1) if p else (0,) for p in atoms.pbc]
    return np.array(list(itertools.product(*rng)), int)


def cage_offsets(atoms, cage, n_iter=20):
    """Integer offsets o_i placing every atom of a cage at its image closest to the cage centre.

    Cartesian nearest image (valid for oblique cells), iterated with the centre.
    Unambiguous while the cage radius is below half the shortest lattice vector
    (checked by `cage_is_compact`).
    """
    L = np.asarray(atoms.cell)
    f = atoms.get_scaled_positions(wrap=False)[cage]
    T = _lattice_shifts(atoms)
    c = f[0] @ L
    off = np.zeros((len(cage), 3), int)
    for _ in range(n_iter):
        cand = (f[:, None, :] + T[None, :, :]) @ L                # (n, shifts, 3)
        off = T[np.argmin(np.linalg.norm(cand - c, axis=2), axis=1)]
        c_new = ((f + off) @ L).mean(0)
        if np.allclose(c_new, c, atol=1e-6):
            break
        c = c_new
    return {int(n): off[k] for k, n in enumerate(cage)}


def internal_graph(G, cage, off):
    """Intra-cage bond graph: bonds joining atoms within the same unwrapped copy.

    A bond (i, j, S) between two atoms of the cage is intra-cage iff
    S = o_j - o_i; otherwise it links the cage to its own periodic image and
    belongs to the inter-cage connection set B.
    """
    H = nx.Graph()
    H.add_nodes_from(cage)
    for i, j, S in oriented_edges(G, cage):
        if i != j and np.array_equal(S, off[j] - off[i]):
            H.add_edge(i, j)
    return H


def cage_is_compact(atoms, cage, off):
    """Cage radius below half the shortest lattice vector: nearest images are unambiguous."""
    L = np.asarray(atoms.cell)
    P = (atoms.get_scaled_positions(wrap=False)[cage] + np.array([off[n] for n in cage])) @ L
    radius = np.linalg.norm(P - P.mean(0), axis=1).max()
    T = _lattice_shifts(atoms)
    lengths = np.linalg.norm(T[np.any(T != 0, axis=1)] @ L, axis=1)
    return bool(radius < 0.5 * lengths.min()) if len(lengths) else True


# ======================================================================= (G3)
def faces(H):
    """Face-size histogram of a connected planar graph, else None."""
    if H.number_of_nodes() < 4 or not nx.is_connected(H):
        return None
    ok, emb = nx.check_planarity(H)
    if not ok:
        return None
    seen, out = set(), []
    for u, v in emb.edges():
        if (u, v) not in seen:
            out.append(len(emb.traverse_face(u, v, mark_half_edges=seen)))
    return Counter(out)


def is_fullerene(H):
    """Phi(H): connected, cubic, planar, faces in {5, 6}; Euler then gives f5 = 12."""
    if H.number_of_nodes() < 20 or any(d != 3 for _, d in H.degree()):
        return False
    f = faces(H)
    return f is not None and set(f) <= {5, 6} and f[5] == 12


def library_graph(n, coef=COEF):
    """Bond graph of the isolated reference molecule C_n (data/C{n}.xyz), if present."""
    if n not in LIBRARY or not os.path.exists(LIBRARY[n]):
        return None
    mol = read(LIBRARY[n])
    mol.set_pbc(False)
    return nx.Graph(remove_selfloops(bond_graph(mol, coef)[0]))


# ======================================================================= (G4)
def gn_recursion(G, atoms):
    """C(X) = {X} if Phi(internal graph of X), else split by Girvan-Newman edge betweenness."""
    queue = [set(c) for c in nx.connected_components(remove_selfloops(G))]
    cages = []
    while queue:
        X = sorted(queue.pop())
        if len(X) >= 20 and is_fullerene(internal_graph(G, X, cage_offsets(atoms, X))):
            cages.append(X)
        elif len(X) < 40:                       # cannot contain two cages: keep as fragment
            cages.append(X)
        else:
            simple = nx.Graph(remove_selfloops(G.subgraph(X)))
            queue.extend(set(p) for p in next(nx.algorithms.community.girvan_newman(simple)))
    return sorted(cages, key=lambda c: c[0])


def _exact_cover(H, refs, natoms):
    """Disjoint induced subgraphs isomorphic to library cages covering all atoms."""
    cands = []
    for ref in refs.values():
        cands += list({frozenset(m) for m in GraphMatcher(H, ref).subgraph_isomorphisms_iter()})
    by_atom = {}
    for k, c in enumerate(cands):
        for a in c:
            by_atom.setdefault(a, []).append(k)
    sol = []

    def search(free):
        if not free:
            return True
        a = min(free, key=lambda x: len(by_atom.get(x, [])))
        for k in by_atom.get(a, []):
            if cands[k] <= free:
                sol.append(k)
                if search(free - cands[k]):
                    return True
                sol.pop()
        return False

    if search(frozenset(range(natoms))):
        return sorted((sorted(int(x) for x in cands[k]) for k in sol), key=lambda c: c[0])
    return None


def identify_cages(atoms, coef=COEF):
    """Certified partition of the atoms into cages -> (cages, how).

    pass 1 'fullerene-certified': every subunit satisfies Phi at the analysis cutoff
    pass 2 'size-certified'     : every subunit is a connected 20/24-atom cage (distorted inputs)
    pass 3 'exact-cover'        : cover by induced subgraphs isomorphic to C20/C24
    otherwise (None, 'unresolved')
    """
    Gval, _ = bond_graph(atoms, coef)
    trials = []
    for lo, hi in LADDER:
        t = gn_recursion(bond_graph(atoms, hi, lo)[0], atoms)
        if all(is_fullerene(internal_graph(Gval, c, cage_offsets(atoms, c))) for c in t):
            return t, 'fullerene-certified'
        trials.append(t)
    for t in trials:
        if all(len(c) in LIBRARY and nx.is_connected(internal_graph(Gval, c, cage_offsets(atoms, c)))
               for c in t):
            return t, 'size-certified'
    refs = {n: g for n in LIBRARY if (g := library_graph(n, coef)) is not None}
    cov = _exact_cover(nx.Graph(remove_selfloops(Gval)), refs, len(atoms))
    if cov:
        return cov, 'exact-cover'
    return None, 'unresolved'


# ======================================================================= (G7)
def cage_status(H, H_input):
    """closed: H isomorphic to the input cage graph (or the library C_n / any fullerene)."""
    n = H.number_of_nodes()
    ref = H_input if (H_input is not None and is_fullerene(H_input)) else library_graph(n)
    closed = nx.is_isomorphic(H, ref) if ref is not None else is_fullerene(H)
    f = faces(H)
    if closed:
        status = 'closed'
    elif not nx.is_connected(H):
        status = 'broken'
    elif f is not None and max(f) > 6:
        status = 'open'
    else:
        status = 'distorted'
    return status, (dict(sorted(f.items())) if f else None)


# ======================================================================= (G5)
def contracted_graph(G, cages, off):
    """Inter-cage connection set B and contracted quotient graph G^ (self-loops allowed).

    Edge (a, b) of G^ carries T = S + o_i - o_j, the translation of cage b's
    copy relative to cage a's, stored with a <= b.
    """
    cage_of = {n: a for a, c in enumerate(cages) for n in c}
    Ghat = nx.MultiGraph()
    Ghat.add_nodes_from(range(len(cages)))
    B = []
    for i, j, S in oriented_edges(G):
        a, b = cage_of[i], cage_of[j]
        T = S + off[i] - off[j]
        if a == b and not T.any():
            continue                                    # intra-cage bond
        if a > b:
            a, b, T = b, a, -T
        if a == b and tuple(T) < (0, 0, 0):
            T = -T
        B.append((i, j, tuple(int(x) for x in S)))
        Ghat.add_edge(a, b, vector=T.astype(int))
    return Ghat, B


# ======================================================================= (G6)
def _labels(Ghat):
    lab = {}
    for a, b, d in Ghat.edges(data=True):
        T = np.asarray(d['vector'], int)
        a, b, T = (a, b, T) if a <= b else (b, a, -T)
        lab.setdefault((a, b), []).append(tuple(T))
    return lab


def _gauge_match(G1, G2, phi):
    """Exists c: multiset{T + c_a - c_b} on (a, b) == labels of (phi a, phi b)?"""
    L1, L2 = _labels(G1), _labels(G2)

    def target(a, b):
        x, y = phi[a], phi[b]
        return sorted(L2.get((x, y), [])) if x <= y else sorted(tuple(-v for v in t) for t in L2.get((y, x), []))

    def ok(a, b, c):
        shifted = sorted(tuple(np.array(t) + c[a] - c[b]) for t in L1[(a, b)])
        if a == b:   # self-loop: S and -S are the same edge
            shifted = sorted(max(t, tuple(-v for v in t)) for t in shifted)
            return shifted == sorted(max(t, tuple(-v for v in t)) for t in target(a, b))
        return shifted == target(a, b)

    simple = nx.Graph()
    simple.add_nodes_from(G1.nodes())
    simple.add_edges_from(k for k in L1 if k[0] != k[1])
    tree = []
    roots = []
    for comp in nx.connected_components(simple):
        r = min(comp)
        roots.append(r)
        tree += list(nx.bfs_edges(simple, r))

    def dfs(k, c):
        if k == len(tree):
            return all(ok(a, b, c) for (a, b) in L1)
        a, b = tree[k]                       # c[a] known, choose c[b]
        key = (a, b) if a <= b else (b, a)
        t0 = np.array(L1[key][0]) if a <= b else -np.array(L1[key][0])
        for m in set(map(tuple, (target(a, b) if a <= b else target(b, a)))):
            m = np.array(m) if a <= b else -np.array(m)
            cb = c[a] + t0 - m               # from t0 + c_a - c_b = m
            c2 = dict(c)
            c2[b] = cb
            if ok(*key, c2) and dfs(k + 1, c2):
                return True
        return False

    return dfs(0, {r: np.zeros(3, int) for r in roots})


def same_connection_mode(r1, r2, use_labels=True, max_iso=20000):
    """(G6) for two `classify` results: VF2 on G^ with cage-type and multiplicity labels,
    plus (use_labels) translation labels up to the per-cage gauge c_a."""
    def simple(r):
        H = nx.Graph()
        for a, k in enumerate(r['cages']):
            H.add_node(a, kind=k)
        for (a, b), ts in _labels(r['contracted_graph']).items():
            H.add_edge(a, b, mult=len(ts))
        return H
    H1, H2 = simple(r1), simple(r2)
    gm = GraphMatcher(H1, H2, node_match=lambda x, y: x['kind'] == y['kind'],
                      edge_match=lambda x, y: x['mult'] == y['mult'])
    if not use_labels:
        return gm.is_isomorphic()
    for k, phi in enumerate(gm.isomorphisms_iter()):
        if _gauge_match(r1['contracted_graph'], r2['contracted_graph'], phi):
            return True
        if k >= max_iso:
            break
    return False


# ================================================================ (G10)-(G11)
def _hull_intersection_depth(A, B):
    """Radius of the largest ball inside conv(A) & conv(B) (LP; < 0 if disjoint)."""
    eqs = np.vstack([ConvexHull(A).equations, ConvexHull(B).equations])   # n.x + h <= 0, |n| = 1
    res = linprog(c=[0, 0, 0, -1], A_ub=np.hstack([eqs[:, :3], np.ones((len(eqs), 1))]),
                  b_ub=-eqs[:, 3], bounds=[(None, None)] * 3 + [(None, 10.0)], method='highs')
    return -res.fun if res.success else -np.inf


def _in_plane_directions(atoms):
    """Periodic directions; c is dropped when it is a vacuum direction (gap > 8 A)."""
    per = [k for k in range(3) if atoms.pbc[k]]
    z = np.sort(atoms.get_scaled_positions()[:, 2])
    gap = np.max(np.diff(np.concatenate([z, z[:1] + 1]))) * atoms.cell.lengths()[2]
    return [k for k in per if not (k == 2 and gap > 8.0)]


def geometric_checks(atoms, cages, off):
    """Endohedral nesting (G10) and hull interpenetration (G11) over periodic images."""
    cell = np.asarray(atoms.cell)
    f = atoms.get_scaled_positions(wrap=False)
    X = [(f[c] + np.array([off[n] for n in c])) @ cell for c in cages]
    cen = [x.mean(0) for x in X]
    rad = [np.linalg.norm(x - x.mean(0), axis=1).max() for x in X]
    per = _in_plane_directions(atoms)
    shifts = [np.array(t) for t in itertools.product(*[(-1, 0, 1) if k in per else (0,) for k in range(3)])]
    endo, overlap, max_depth = [], [], 0.0
    for a, b in itertools.combinations(range(len(cages)), 2):
        for t in shifts:
            Xb = X[b] + t @ cell
            dc = float(np.linalg.norm(cen[b] + t @ cell - cen[a]))
            if dc > rad[a] + rad[b]:
                continue
            in_ab = Delaunay(X[a]).find_simplex(Xb) >= 0
            in_ba = Delaunay(Xb).find_simplex(X[a]) >= 0
            if in_ab.all() or in_ba.all() or dc < CENTER_TOL:
                outer, inner = (a, b) if rad[a] >= rad[b] else (b, a)
                endo.append(dict(outer=outer, inner=inner, shift=tuple(int(x) for x in t),
                                 centroid_distance=round(dc, 3)))
                continue
            depth = _hull_intersection_depth(X[a], Xb)
            max_depth = max(max_depth, depth)
            if depth > HULL_EPS:
                overlap.append(dict(a=a, b=b, shift=tuple(int(x) for x in t),
                                    atoms_inside=int(in_ab.sum() + in_ba.sum()),
                                    depth=round(float(depth), 3), centroid_distance=round(dc, 3)))
    return endo, overlap, max_depth


# ================================================================ (G13)-(G14)
def symmetry_checks(atoms, G, cages, B, symprec=SYMPREC):
    """Layer-group action on G, partition equivariance and orbits of B."""
    import spglib
    cell = (np.asarray(atoms.cell), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    try:
        ds = spglib.get_layergroup(cell, aperiodic_dir=2, symprec=symprec)
    except Exception:
        return None
    if ds is None:
        return None
    frac = atoms.get_scaled_positions()
    lat = np.asarray(atoms.cell)
    edges = Counter()
    for i, j, S in oriented_edges(G):
        edges[(i, j, tuple(S))] += 1
        if i != j:
            edges[(j, i, tuple(-S))] += 1
    part = {frozenset(c) for c in cages}
    graph_ok = part_ok = True
    ops = []
    for R, tau in zip(ds.rotations, ds.translations):
        perm, shift = np.empty(len(atoms), int), np.zeros((len(atoms), 3), int)
        for i, x in enumerate(frac @ R.T + tau):
            d = x - frac
            dist = np.linalg.norm((d - np.round(d)) @ lat, axis=1)
            k = int(np.argmin(dist))
            if dist[k] > 3 * symprec:
                return dict(layer_group=ds.international, layer_no=int(ds.number),
                            order=len(ds.rotations), mapped=False)
            perm[i], shift[i] = k, np.round(d[k]).astype(int)
        ops.append((R, perm, shift))
        for (i, j, S), m in edges.items():
            S2 = tuple(int(x) for x in R @ np.array(S) + shift[j] - shift[i])
            if edges.get((perm[i], perm[j], S2), 0) != m:
                graph_ok = False
        if {frozenset(int(perm[n]) for n in c) for c in cages} != part:
            part_ok = False

    def key(i, j, S):
        return min((i, j, S), (j, i, tuple(-x for x in S)))
    parent = {key(i, j, S): key(i, j, S) for i, j, S in B}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for R, perm, shift in ops:
        for (i, j, S) in list(parent):
            b2 = key(int(perm[i]), int(perm[j]), tuple(int(x) for x in R @ np.array(S) + shift[j] - shift[i]))
            if b2 in parent:
                parent[find((i, j, S))] = find(b2)
    orbits = Counter(find(b) for b in parent)
    return dict(layer_group=ds.international, layer_no=int(ds.number), order=len(ds.rotations),
                mapped=True, graph_invariant=graph_ok, partition_equivariant=part_ok,
                B_orbits=len(orbits), B_orbit_sizes=sorted(orbits.values(), reverse=True))


# ======================================================================= (G12)
def classify(structure, ref=None, coef=COEF, symprec=SYMPREC, symmetry=True):
    """Identify, validate and classify one structure.

    structure : path of the structure to classify (e.g. a relaxed CONTCAR)
    ref       : optional unrelaxed input with the same atom order; cage identity is
                taken from it, so cages that open during relaxation are still tracked
    """
    atoms = load(structure)
    src = load(ref) if ref else atoms
    if len(src) != len(atoms):
        raise ValueError('ref and structure must contain the same atoms in the same order')
    G, clashes = bond_graph(atoms, coef)
    out = dict(structure=str(structure), ref=str(ref) if ref else None, natoms=len(atoms), coef=coef)

    cages, how = identify_cages(src, coef)
    out['identification'] = how

    Z = atoms.get_atomic_numbers()
    deg = Counter()
    for i, j, _ in oriented_edges(G):
        deg[i] += 1
        deg[j] += 1
    over_valent = [n for n in range(len(atoms)) if deg[n] > MAX_VALENCE.get(int(Z[n]), 8)]
    comps = [sorted(c) for c in nx.connected_components(remove_selfloops(G))]
    comp_dim = {n: get_dimension(G.subgraph(c)) for c in comps for n in c}
    out.update(clashes=len(clashes),
               min_clash_distance=round(min(c[3] for c in clashes), 3) if clashes else None,
               over_valent_atoms=len(over_valent), dim_G=get_dimension(G),
               components=[(len(c), comp_dim[c[0]]) for c in comps])
    flags = (['atomic clash'] if clashes else []) + (['valence violation'] if over_valent else [])

    if cages is None:
        out.update(flags=flags + ['unresolved cage identity'], label='Unresolved input')
        return out

    G_src = bond_graph(src, coef)[0] if ref else G
    kinds = [f'C{len(c)}' for c in cages]
    status, face_sets, offs = [], [], {}
    for c in cages:
        off = cage_offsets(atoms, c)
        offs.update(off)
        if not cage_is_compact(atoms, c, off):
            flags.append(f'C{len(c)} cage radius exceeds half the shortest lattice vector: unwrapping unreliable')
        H_in = internal_graph(G_src, c, cage_offsets(src, c))
        s, fc = cage_status(internal_graph(G, c, off), H_in)
        status.append(s)
        face_sets.append(fc)
    Ghat, B = contracted_graph(G, cages, offs)
    out.update(cages=kinds, cage_atoms=cages, cage_status=status, cage_faces=face_sets,
               B=len(B), contracted_graph=Ghat,
               B_pairs=dict(Counter('-'.join(sorted((kinds[a], kinds[b]))) for a, b in Ghat.edges())),
               cage_bonds=[Ghat.degree(a) for a in Ghat.nodes()],
               dim_Ghat=get_dimension(Ghat) if Ghat.number_of_edges() else 0)

    not_2d = [a for a, c in enumerate(cages) if comp_dim[c[0]] < 2]
    if out['dim_G'] != 2 or not_2d:
        flags.append(f'not a 2D covalent network ({len(not_2d)} cage(s) in D<2 components)')
    endo, overlap, depth = geometric_checks(atoms, cages, offs)
    out['max_hull_depth'] = round(float(depth), 3)
    out['endohedral'] = [f"{kinds[e['inner']]}@{kinds[e['outer']]} (centroid distance "
                         f"{e['centroid_distance']} A)" for e in endo]
    out['overlap'] = [f"{kinds[o['a']]}/{kinds[o['b']]}{'' if not any(o['shift']) else ' image ' + str(o['shift'])}: "
                      f"{o['atoms_inside']} atoms inside the other hull, depth {o['depth']} A" for o in overlap]
    if endo:
        flags.append('endohedral cage-in-cage')
    if overlap:
        flags.append('overlapping cages')
    out['symmetry'] = symmetry_checks(atoms, G, cages, B, symprec) if symmetry else None

    n_closed = status.count('closed')
    if endo:
        label = 'Endohedral cage-in-cage (rejected)'
    elif overlap or clashes:
        label = 'Overlapping cages (rejected)'
    elif out['dim_G'] != 2 or not_2d:
        label = 'Non-2D / van der Waals-packed (rejected)'
    elif n_closed == len(cages):
        label = 'Closed-cage network'
    elif n_closed == 0:
        label = 'Fully opened'
    else:
        label = 'Partially opened'
    out.update(flags=flags, n_closed=n_closed, n_cages=len(cages), label=label)
    return out


def report(r):
    """Human-readable summary of a `classify` result."""
    w = 76
    lines = ['=' * w, f"{r['structure']}  ({r['natoms']} atoms, k = {r['coef']}"
             + (f", cage identity from {r['ref']}" if r['ref'] else '') + ')', '-' * w,
             f"  label               : {r['label']}",
             f"  cage identification : {r['identification']}"]
    if 'cages' in r:
        lines += [f"  cages (status)      : {', '.join(f'{k}:{s}' for k, s in zip(r['cages'], r['cage_status']))}",
                  f"  closed / all        : {r['n_closed']} / {r['n_cages']}",
                  f"  |B| by cage pair    : {r['B']}  {r['B_pairs']}",
                  f"  bonds per cage      : {r['cage_bonds']}   dim(G^) = {r['dim_Ghat']}"]
    lines += [f"  dim(G), components  : {r['dim_G']}, {r['components']}",
              f"  clashes (min d, A)  : {r['clashes']} ({r['min_clash_distance']})   "
              f"over-valent atoms: {r['over_valent_atoms']}"]
    for k in ('endohedral', 'overlap'):
        for n, s in enumerate(r.get(k) or []):
            lines.append(f"  {k if n == 0 else '':<20}: {s}")
    s = r.get('symmetry')
    if s and s.get('mapped'):
        lines += [f"  layer group         : {s['layer_group']} (No. {s['layer_no']}, order {s['order']}); "
                  f"sigma(G)=G {s['graph_invariant']}; C(sG)=sC(G) {s['partition_equivariant']}",
                  f"  B symmetry orbits   : {s['B_orbits']}  sizes {s['B_orbit_sizes']}"]
    elif s:
        lines.append(f"  layer group         : {s['layer_group']} (operations not resolvable at symprec)")
    lines.append(f"  flags               : {', '.join(r['flags']) if r['flags'] else 'none'}")
    return '\n'.join(lines)
