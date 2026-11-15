"""Exact quadratic energies on nested integer counts with proved implications.

Coefficient order is [a,d_0,...,d_(k-1),B_00,B_01,...] with upper-triangular
B entries; q(N)=a+2*d.N+N^T*B*N.  This helper exhausts the finite integer
domain.  The calling global verifier must independently verify the ordinary
count bounds and every affine geometric implication before using q>=0.
"""
from fractions import Fraction as Q
import itertools
import math
import numpy as np


class ConditionalCountPolynomial:
    def __init__(self, record):
        self.radii = tuple(map(Q, record['radii']))
        self.bounds = tuple(record['bounds'])
        self.k = len(self.radii)
        assert 1 <= self.k <= 6 and len(self.bounds) == self.k
        assert all(isinstance(m, int) and 0 <= m <= 30 for m in self.bounds)
        assert all(a < b for a, b in zip(self.radii, self.radii[1:]))
        assert all(a <= b for a, b in zip(self.bounds, self.bounds[1:]))
        self.upper = tuple(itertools.combinations_with_replacement(range(self.k), 2))
        self.nv = 1+self.k+len(self.upper)
        self.coefficients = tuple(map(Q, record['coefficients']))
        assert len(self.coefficients) == self.nv
        self.implications = []
        for implication in record['implications']:
            i, j = implication['inner_index'], implication['outer_index']
            assert isinstance(i, int) and isinstance(j, int) and 0 <= i < j < self.k
            assert implication['maximum'] == 24 and self.bounds[i] <= 24
            proof = implication['proof']
            assert Q(proof['inner_radius']) == self.radii[i]
            assert Q(proof['outer_radius']) == self.radii[j]
            assert implication['proof_kind'] == 'affine_latitude'
            self.implications.append((i, j))

    def counts(self):
        for count in itertools.combinations_with_replacement(range(max(self.bounds)+1), self.k):
            if any(n > bound for n, bound in zip(count, self.bounds)):
                continue
            if any(count[i] == 24 and count[j] != 24 for i, j in self.implications):
                continue
            yield count

    def count_row(self, count):
        return [1]+[2*x for x in count]+[(1 if i == j else 2)*count[i]*count[j]
                                        for i, j in self.upper]

    def pair_row(self, flags):
        assert len(flags) == self.k and all(x in [0, 1] for x in flags)
        return [0]+[2*x for x in flags]+[(1 if i == j else 2)*flags[i]*flags[j]
                                       for i, j in self.upper]

    def triangle_row(self, flags):
        assert len(flags) == 3 and all(len(row) == self.k for row in flags)
        assert all(x in [0, 1] for row in flags for x in row)
        result = [0]*self.nv
        for a, b in [(0, 1), (0, 2), (1, 2)]:
            for index, (i, j) in enumerate(self.upper, 1+self.k):
                result[index] += (flags[a][i]*flags[b][j] if i == j else
                                  flags[a][i]*flags[b][j]+flags[a][j]*flags[b][i])
        return result

    def flags(self, breaks, region):
        lo, hi = breaks[region:region+2]
        assert all(hi <= r or lo >= r for r in self.radii)
        return [int(hi <= r) for r in self.radii]

    def verify(self):
        if not __debug__:
            raise RuntimeError('Run without -O.')
        denominator = math.lcm(*(v.denominator for v in self.coefficients))
        integers = [int(v*denominator) for v in self.coefficients]
        counts = list(self.counts())
        rows = np.array([self.count_row(n) for n in counts], dtype=np.int64)
        # Prove absence of integer overflow before using NumPy's fast path.
        bound = sum(int(np.max(abs(rows[:, i])))*abs(value) for i, value in enumerate(integers))
        dtype = np.int64 if bound < 2**63 else object
        values = rows.astype(dtype)@np.array(integers, dtype=dtype)
        index = int(np.argmin(values)); minimum = Q(int(values[index]), denominator)
        assert minimum >= 0, (minimum, counts[index])
        assert self.coefficients[0] >= 0  # Also follows from the zero count.
        return dict(status='EXACT_CONDITIONAL_INTEGER_POLYNOMIAL_NONNEGATIVE',
                    count_vectors_checked=len(counts), minimum=str(minimum),
                    minimizer=list(counts[index]), common_denominator=str(denominator),
                    integer_accumulation_bound=str(bound),
                    dependency='Every count bound and affine implication requires independent geometric verification.')
