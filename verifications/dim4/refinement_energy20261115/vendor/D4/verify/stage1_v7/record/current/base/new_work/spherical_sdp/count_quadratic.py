"""Exact nonnegative quadratic cone for nested, integer neighbor counts.

Coefficient order is [a, d_0,...,d_(k-1), B_00,B_01,...,B_(k-1,k-1)],
with the symmetric upper triangle in row-major order. The polynomial is
q(n)=a+2*d.n+n^T*B*n. No PSD condition is imposed on B.

Use cone.feasibility_rows() @ w >= 0 in a numerical optimizer. Pair and cyclic
rows use the same w. After rationalizing w, cone.verify(w) checks every feasible
integer count vector exactly. The other global packing hypotheses still need
their own exact verification.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


class NestedCountCone:
    def __init__(self, radii, bounds):
        self.radii = tuple(map(F, radii))
        self.bounds = tuple(bounds)
        self.k = len(self.radii)
        assert self.k == len(self.bounds) and self.k > 0
        assert all(isinstance(m, int) and m >= 0 for m in self.bounds)
        assert all(a < b for a, b in zip(self.radii, self.radii[1:]))
        assert all(a <= b for a, b in zip(self.bounds, self.bounds[1:]))
        self.upper = tuple((i, j) for i in range(self.k) for j in range(i, self.k))
        self.nv = 1+self.k+len(self.upper)
        self.counts = tuple(n for n in itertools.product(*(range(m+1) for m in self.bounds))
                            if all(a <= b for a, b in zip(n, n[1:])))

    @classmethod
    def accepted(cls, path=ROOT/'accepted_count_envelope.json'):
        data = json.loads(Path(path).read_text())
        assert data['status'] == 'EXACT_COUNT_ENVELOPE'
        assert data['dimension'] == 4 and data['minimum_center_separation'] == '1'
        rows = sorted(data['widest_per_count'], key=lambda r: F(r['R']))
        return cls([r['R'] for r in rows], [r['M'] for r in rows])

    def row_for_count(self, n):
        return [1]+[2*x for x in n]+[(1 if i == j else 2)*n[i]*n[j]
                                         for i, j in self.upper]

    def feasibility_rows(self):
        return np.array([self.row_for_count(n) for n in self.counts], dtype=np.int64)

    def pair_row_from_flags(self, flags):
        assert len(flags) == self.k and all(x in (0, 1) for x in flags)
        return [0]+[2*x for x in flags]+[(1 if i == j else 2)*flags[i]*flags[j]
                                             for i, j in self.upper]

    def triangle_row_from_flags(self, flags):
        assert len(flags) == 3 and all(len(f) == self.k for f in flags)
        out = [0]*self.nv
        for a, b in [(0, 1), (0, 2), (1, 2)]:
            for p, (i, j) in enumerate(self.upper, 1+self.k):
                out[p] += (flags[a][i]*flags[b][j] if i == j else
                           flags[a][i]*flags[b][j]+flags[a][j]*flags[b][i])
        return out

    def pair_rows(self, distances):
        # Numerical proposal helper. Exact branch proofs use *_from_flags.
        flags = np.asarray(distances)[:, None] <= np.array(list(map(float, self.radii)))
        return np.array([self.pair_row_from_flags(list(map(int, f))) for f in flags])

    def triangle_rows(self, side_triples):
        flags = np.asarray(side_triples)[:, :, None] <= np.array(list(map(float, self.radii)))
        return np.array([self.triangle_row_from_flags(f.astype(int).tolist()) for f in flags])

    def origin_row(self):
        return [1]+[0]*(self.nv-1)

    def verify(self, coefficients):
        w = list(map(F, coefficients))
        assert len(w) == self.nv
        values = [(sum((F(x)*y for x, y in zip(self.row_for_count(n), w)), F(0)), n)
                  for n in self.counts]
        minimum, witness = min(values)
        assert minimum >= 0, ('Negative count energy', minimum, witness)
        return {'status': 'EXACT_NESTED_COUNT_POLYNOMIAL_NONNEGATIVE',
                'count_vectors_checked': len(values), 'minimum': str(minimum),
                'minimizer': list(witness)}

    def cross_product(self, i, j):
        """(M_i-N_i)*N_j: nonnegative by the cap bound and N_j>=0."""
        assert 0 <= i < self.k and 0 <= j < self.k
        w = [F(0)]*self.nv
        w[1+j] = F(self.bounds[i], 2)
        u = (min(i, j), max(i, j))
        w[1+self.k+self.upper.index(u)] = F(-1 if i == j else -1, 1 if i == j else 2)
        self.verify(w)
        return w

    def cap_times_shell(self, i, j, k):
        """(M_i-N_i)*(N_k-N_j) for j<k."""
        assert 0 <= j < k < self.k
        return [a-b for a, b in zip(self.cross_product(i, k), self.cross_product(i, j))]

    def integer_gap(self, m, k):
        """(m.n-k)*(m.n-k-1), nonnegative for integral m,n,k."""
        assert len(m) == self.k and all(isinstance(x, int) for x in m) and isinstance(k, int)
        return [F(k*(k+1))]+[F(-(2*k+1)*x, 2) for x in m]+[F(m[i]*m[j]) for i, j in self.upper]


def main():
    cone = NestedCountCone.accepted()
    generators = []
    for i, j in itertools.product(range(cone.k), repeat=2):
        w = cone.cross_product(i, j)
        generators.append({'polynomial': '(M_'+str(i)+'-N_'+str(i)+')*N_'+str(j),
                           'coefficients': list(map(str, w)), 'verification': cone.verify(w)})
    for i in range(cone.k):
        for j in range(cone.k):
            for k in range(j+1, cone.k):
                w = cone.cap_times_shell(i, j, k)
                generators.append({'polynomial': '(M_'+str(i)+'-N_'+str(i)+')*(N_'+str(k)+'-N_'+str(j)+')',
                                   'coefficients': list(map(str, w)), 'verification': cone.verify(w)})
    data = {'status': 'EXACT_NESTED_COUNT_CONE_SPECIFICATION',
            'radii': list(map(str, cone.radii)), 'bounds': list(cone.bounds),
            'coefficient_count': cone.nv, 'integer_count_vectors': len(cone.counts),
            'coefficient_order': ['a']+['d_'+str(i) for i in range(cone.k)]+
                                  ['B_'+str(i)+str(j) for i, j in cone.upper],
            'generators': generators}
    (ROOT/'nested_count_cone.json').write_text(json.dumps(data, indent=2)+'\n')
    print('Exact nested count cone:', cone.nv, 'coefficients,', len(cone.counts),
          'feasible integer vectors,', len(generators), 'verified product rays.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run without -O.')
    main()
