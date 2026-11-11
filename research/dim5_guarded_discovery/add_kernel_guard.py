"""Add an exact PSD annulus-constant guard. Full verification is mandatory."""
if not __debug__:
    raise RuntimeError('Verification must not run under python -O.')
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as P

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('candidate', type=Path)
ap.add_argument('guard', type=Path)
ap.add_argument('epsilon', type=F)
ap.add_argument('output', type=Path)
args = ap.parse_args()
x = json.loads(args.candidate.read_text())[0]
g = json.loads(args.guard.read_text())
assert x['dimension'] == g['dimension'] and F(x['truncation_squared']) == F(9, 4)
assert F(x['exclusion_squared']) == 1 and args.epsilon >= 0
cs = list(map(Q, g['gegenbauer_coefficients']))
assert cs[0] == Q(1, 2) and min(cs) >= 0
assert len(cs) <= len(x['matrices'])
n = x['dimension']
X = P([0, 1])
G = [P([1]), X]
for k in range(1, len(cs) - 1):
    G.append(((2*k+n-2)*X*G[-1]-k*G[-2])/Q(k+n-2))
pol = sum((a*b for a, b in zip(cs, G)), P([])) + Q(1, 2)
seq = [pol, pol.derivative()]
while seq[-1].degree() > 0:
    rem = -(seq[-2] % seq[-1])
    if not rem:
        break
    z = rem.numer()
    seq.append(P(z / abs(z.content())))
def variation(t):
    signs = [1 if q(t) > 0 else -1 for q in seq if q(t)]
    return sum(a != b for a, b in zip(signs, signs[1:]))
assert pol(-1) < 0 and pol(Q(7, 9)) < 0
assert variation(Q(-1)) == variation(Q(7, 9))
D = sum(cs, Q(0))
a = D*D/(4*cs[0])
b = -D/2
assert a*cs[0] == b*b and 2*b+D == 0
eps = Q(str(args.epsilon))
M = [[[Q(str(F(str(v)))) for v in row] for row in A] for A in x['matrices']]
M[0][0][0] += eps*a
M[0][0][1] += eps*b
M[0][1][0] += eps*b
M[0][1][1] += eps*cs[0]
for k in range(1, len(cs)):
    M[k][0][0] += eps*cs[k]
x['matrices'] = [[[str(v) for v in row] for row in A] for A in M]
x['kernel_guard'] = {'source': args.guard.name, 'epsilon': str(args.epsilon),
                     'exact_pair_correction': '0', 'uniform_sign_upper': '-1/2'}
x['certified'] = False
x['warning'] = 'Exact kernel guard added; complete fresh verification is mandatory.'
for key in ['objective', 'density_candidate', 'max_constraint_violation',
            'status', 'history', 'solver_iterations', 'solver_seconds']:
    x.pop(key, None)
args.output.write_text(json.dumps([x], indent=2) + '\n')
