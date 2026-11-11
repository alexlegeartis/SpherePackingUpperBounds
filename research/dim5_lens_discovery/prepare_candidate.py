"""Convert a numerical kernel to exact rational symmetric PD source blocks.

Only symmetrizes exact finite decimals and sets the isolated origin to the
exact Schur value plus 1e-12. No off-origin shift is permitted here. Every
complete continuum condition must subsequently be checked. Never treat this
operation as certification.
"""
if not __debug__:
    raise RuntimeError('Verification must not run under python -O.')
import argparse
import json
from fractions import Fraction
from pathlib import Path
from flint import fmpq as Q, fmpq_mat as Mat


def fq(value):
    return Q(str(Fraction(str(value))))


def positive_definite(A):
    N = A.nrows()
    L = Mat(N, N)
    ds = []
    for j in range(N):
        L[j, j] = 1
        d = A[j, j] - sum((L[j, k] ** 2 * ds[k] for k in range(j)), Q(0))
        if d <= 0:
            return False
        ds.append(d)
        for i in range(j + 1, N):
            L[i, j] = (A[i, j] - sum((L[i, k] * L[j, k] * ds[k]
                                     for k in range(j)), Q(0))) / d
    return True


ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('source', type=Path)
ap.add_argument('output', type=Path)
args = ap.parse_args()
if args.output.exists():
    ap.error('Output already exists; choose a new file')
x = json.loads(args.source.read_text())[0]
matrices = []
repairs = []
for ell, raw in enumerate(x['matrices']):
    N = len(raw)
    A = Mat([[(fq(raw[i][j]) + fq(raw[j][i])) / 2
              for j in range(N)] for i in range(N)])
    offset = int(ell == 0)
    B = Mat([[A[i + offset, j + offset] for j in range(N - offset)]
             for i in range(N - offset)])
    shift = Q(0)
    if not positive_definite(B):
        raise RuntimeError(f'Block {ell} requires an off-origin shift; this package forbids it')
    if ell == 0:
        b = Mat([[A[i + 1, 0]] for i in range(N - 1)])
        A[0, 0] = (b.transpose() * B.solve(b))[0, 0] + Q(1, 10**12)
    assert positive_definite(A)
    matrices.append([[str(A[i, j]) for j in range(N)] for i in range(N)])
    repairs.append({'block': ell, 'off_origin_diagonal_shift': str(shift)})
x['matrices'] = matrices
x['exact_source_preparation'] = {
    'source': args.source.name,
    'symmetrization': 'Exact average of the two finite-decimal matrix entries',
    'origin_replacement': 'Exact Schur value plus 1/10^12',
    'repairs': repairs,
    'warning': 'No complete geometric condition was certified by this operation.'}
args.output.write_text(json.dumps([x], indent=2) + '\n')
print(json.dumps(x['exact_source_preparation'], indent=2))
