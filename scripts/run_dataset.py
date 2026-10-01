"""Classify all structures in dataset/structures and group them into connection modes.

    python scripts/run_dataset.py [--nproc 4]

Each structure folder holds POSCAR (unrelaxed input, fixes cage identity),
CONTCAR (relaxed), energy.txt (E0, eV) and gap.txt (PBE gap, eV).  Writes
dataset/graph_classification.csv (one row per structure, joined with
dataset/dataset_index.csv).
"""
import os
import sys
import glob
import argparse
from collections import Counter
from multiprocessing import Pool

import pandas as pd
import networkx as nx

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from heterocage_graph import classify, same_connection_mode  # noqa: E402

INT_COLUMNS = ('n_cages', 'n_closed', 'n_open', 'n_broken', 'n_distorted', 'B', 'B_C20_C20',
               'B_C20_C24', 'B_C24_C24', 'dim_G', 'dim_Ghat', 'clashes', 'over_valent_atoms',
               'layer_group_no', 'B_symmetry_orbits')


def run_one(folder):
    try:
        return os.path.basename(folder), classify(os.path.join(folder, 'CONTCAR'),
                                                  os.path.join(folder, 'POSCAR'))
    except Exception as exc:                                   # keep the batch alive
        return os.path.basename(folder), dict(label='ERROR', flags=[repr(exc)])


def wl_key(r):
    """Weisfeiler-Lehman hash of G^ (cage type, multiplicity): buckets before exact VF2."""
    H = nx.Graph()
    for a, k in enumerate(r['cages']):
        H.add_node(a, kind=k)
    for a, b in r['contracted_graph'].edges():
        H.add_edge(a, b, mult=H.edges[a, b]['mult'] + 1 if H.has_edge(a, b) else 1)
    for a, b in H.edges():
        H.edges[a, b]['mult'] = str(H.edges[a, b]['mult'])
    return nx.weisfeiler_lehman_graph_hash(H, node_attr='kind', edge_attr='mult')


def group_modes(results, use_labels):
    """Connection-mode index per structure (G6); -1 if cage identity is unresolved."""
    reps, mode = [], {}
    for sid, r in results.items():
        if 'contracted_graph' not in r:
            mode[sid] = -1
            continue
        h = wl_key(r)
        for k, (hk, rep) in enumerate(reps):
            if hk == h and same_connection_mode(r, results[rep], use_labels):
                mode[sid] = k
                break
        else:
            reps.append((h, sid))
            mode[sid] = len(reps) - 1
    return mode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    folders = sorted(f.rstrip('/') for f in glob.glob(os.path.join(ROOT, 'dataset', 'structures', '*', '')))
    with Pool(a.nproc) as p:
        results = dict(p.map(run_one, folders, chunksize=2))
    topo = group_modes(results, use_labels=False)
    full = group_modes(results, use_labels=True)

    rows = []
    for sid, r in sorted(results.items()):
        st = Counter(r.get('cage_status', []))
        pairs = r.get('B_pairs', {})
        sym = r.get('symmetry') or {}
        has_b = 'B' in r
        rows.append(dict(
            id=sid, label=r['label'], identification=r.get('identification'),
            n_cages=r.get('n_cages'), n_closed=st['closed'] if has_b else None,
            n_open=st['open'] if has_b else None, n_broken=st['broken'] if has_b else None,
            n_distorted=st['distorted'] if has_b else None, B=r.get('B'),
            B_C20_C20=pairs.get('C20-C20', 0) if has_b else None,
            B_C20_C24=pairs.get('C20-C24', 0) if has_b else None,
            B_C24_C24=pairs.get('C24-C24', 0) if has_b else None,
            dim_G=r.get('dim_G'), dim_Ghat=r.get('dim_Ghat'), clashes=r.get('clashes'),
            max_hull_depth_A=r.get('max_hull_depth'),
            over_valent_atoms=r.get('over_valent_atoms'),
            endohedral='; '.join(r.get('endohedral') or []), overlap='; '.join(r.get('overlap') or []),
            layer_group=sym.get('layer_group'), layer_group_no=sym.get('layer_no'),
            B_symmetry_orbits=sym.get('B_orbits'),
            connection_mode_topology=topo[sid], connection_mode_labelled=full[sid],
            flags='; '.join(r.get('flags', []))))
    df = pd.read_csv(os.path.join(ROOT, 'dataset', 'dataset_index.csv')).merge(
        pd.DataFrame(rows), on='id', how='left')
    for c in INT_COLUMNS:
        df[c] = df[c].astype('Int64')
    out = os.path.join(ROOT, 'dataset', 'graph_classification.csv')
    df.to_csv(out, index=False)
    print(df['label'].value_counts().to_string())
    ok = df['connection_mode_topology'] >= 0
    print(f"connection modes among {ok.sum()} resolved structures: "
          f"{df.loc[ok, 'connection_mode_topology'].nunique()} (cage type + multiplicity), "
          f"{df.loc[ok, 'connection_mode_labelled'].nunique()} (+ translation labels up to gauge)")
    print('written', os.path.relpath(out, ROOT))


if __name__ == '__main__':
    main()
