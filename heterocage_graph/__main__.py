"""Command line:  python -m heterocage_graph STRUCTURE [--ref UNRELAXED] [--coef 1.2] [--json]"""
import argparse
import json

from .classifier import classify, report, COEF, SYMPREC


def main():
    ap = argparse.ArgumentParser(prog='python -m heterocage_graph', description=__doc__)
    ap.add_argument('structure', help='structure file readable by ASE (POSCAR, CONTCAR, cif, ...)')
    ap.add_argument('--ref', help='unrelaxed input with the same atom order; fixes cage identity')
    ap.add_argument('--coef', type=float, default=COEF, help='bond factor k in d < k (r_i + r_j)')
    ap.add_argument('--symprec', type=float, default=SYMPREC, help='spglib tolerance (A)')
    ap.add_argument('--json', action='store_true', help='print the result as JSON')
    a = ap.parse_args()
    r = classify(a.structure, a.ref, a.coef, a.symprec)
    if a.json:
        r.pop('contracted_graph', None)
        print(json.dumps(r, indent=1, default=str))
    else:
        print(report(r))


if __name__ == '__main__':
    main()
