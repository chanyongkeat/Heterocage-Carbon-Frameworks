"""Classify every structure in examples/ and write examples/expected_output.txt.

    python scripts/run_examples.py
"""
import os
import sys
import glob

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from heterocage_graph import classify, report  # noqa: E402


def jobs():
    ex = os.path.join(ROOT, 'examples')
    for d in sorted(glob.glob(os.path.join(ex, 'pareto', '*', ''))):
        yield os.path.join(d, 'CONTCAR'), os.path.join(d, 'POSCAR')       # relaxed + input
    for case in ('endohedral', 'overlap'):
        yield os.path.join(ex, case, 'POSCAR'), None                      # single structure


def main():
    os.chdir(ROOT)
    blocks, summary = [], []
    for s, ref in jobs():
        r = classify(os.path.relpath(s, ROOT), os.path.relpath(ref, ROOT) if ref else None)
        blocks.append(report(r))
        summary.append(f"{os.path.relpath(s, ROOT):<42}{r['label']}")
        print(blocks[-1], flush=True)
    table = '\n'.join(['=' * 76, f"{'structure':<42}label", *summary])
    print(table)
    with open(os.path.join(ROOT, 'examples', 'expected_output.txt'), 'w') as f:
        f.write('\n'.join(blocks) + '\n' + table + '\n')


if __name__ == '__main__':
    main()
