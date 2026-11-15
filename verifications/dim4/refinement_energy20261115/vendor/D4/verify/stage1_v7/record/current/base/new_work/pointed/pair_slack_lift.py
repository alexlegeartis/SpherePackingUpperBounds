"""Exact enclosures for recycling a verified outer pair slack.

This is a local construction, not a global density certificate.  Suppose
the unlifted kernel has zero origin coupling on (7/4,2], positive-semidefinite
diagonal D(r), and f(r)+D(r)<=0 there.  Put w(r)=-f(r)-D(r) on that annulus,
zero elsewhere.  For a unit-separated packing, any ball of squared radius
a2<1/4 contains at most one center.  Averaging such balls with centers at
distance one from the origin, pointed relative to each pivot y, gives

  sum_y w(|y|) [1 - occupancy of the pointed ball] >= 0.

In that example the origin is outside every probe.  An outer pivot cannot select itself,
since its distance to a probe center is at least |y|-1>3/4.  The new energy
has zero constant term, pair increment w(r), and nonpositive cyclic terms.
Thus it cancels the outer pair slack and leaves the density objective
unchanged.  Fixed convex combinations of probe directions are allowed.

The implementation also permits arbitrary fixed probe distances and
capacities M, with the usual variance or R4 obtuse-vector occupancy proof.
It checks self exclusion by the exact quadratic inequality
  r^2-2*d*sigma*r+d^2-a2>0 on [7/4,2].
Each normalized convex weight theta contributes theta*w/M times the
M-point occupancy energy, hence pair increment theta*w.  Distant probes
can select outer neighbors as well as inner ones.

The radial function is Gaussian-polynomial plus polynomial.  The code
below encloses it using outward integer Bernstein coefficients and exact
Taylor bounds for the exponential.  It never replaces a negative sampled
pair value with a proof of the whole annulus.  The optional nonnegativity
clamp is enabled only after the calling global checker has itself proved
every unlifted outer pair interval.
"""
from fractions import Fraction as Q
from functools import lru_cache
import itertools
import math
import numpy as np

from interval_certificate_general_closed import probability_bounds
from interval_certificate_full_closed import angle_bounds
from pair_interval_certificate_closed import exponential_bounds


def product_interval(a, b):
    values = [x*y for x, y in itertools.product(a, b)]
    return min(values), max(values)


def validate_probes(records, core):
    """Validate fixed convex weights and exact small-ball capacities."""
    assert core.BREAKS[-1] == 2 and Q(7, 4) in core.BREAKS
    probes = []
    for row in records:
        assert row['budget'] == 'NEGATIVE_UNLIFTED_OUTER_PAIR'
        assert Q(row['inner_radius']) == Q(7, 4)
        assert Q(row['outer_radius']) == 2
        d, a2, sigma, fraction = [Q(row[k]) for k in ['d', 'a2', 'sigma', 'fraction']]
        M = int(row['M']); assert Q(row['M']) == M and M >= 1
        assert d > 0 and a2 > 0 and -1 <= sigma <= 1 and fraction >= 0
        anchor = int(d*d <= a2)
        total = M+anchor
        variance = a2 < Q(total, 2*(total+1))
        rankin = total == 5 and a2 < Q(1, 2)
        assert variance or rankin, 'No verified occupancy bound for this probe.'
        vertex = min(Q(2), max(Q(7, 4), d*sigma))
        self_gap = vertex*vertex-2*d*sigma*vertex+d*d-a2
        assert self_gap > 0, 'An outer pivot could select itself.'
        probes.append(dict(M=M, d=d, a2=a2, sigma=sigma, fraction=fraction,
                           self_exclusion_squared_gap=self_gap,
                           anchor_in_ball=bool(anchor),
                           capacity_proof='variance' if variance else 'R4_obtuse_vectors'))
    assert sum((p['fraction'] for p in probes), Q(0)) <= 1
    return probes


class ScalarPairBudget:
    """Exact bounds for -f-D on existing radial branches.

    The actual function is zero at r=7/4 and positive-support branches
    start immediately to its right.  Checking both closed branch extensions
    is sufficient for this one-sided convention, even when w jumps there.
    """

    def __init__(self, core, scalar):
        self.core = core
        self.alpha = Q(scalar['alpha'])
        assert self.alpha == 2 and core.BREAKS[-1] == 2
        assert Q(7, 4) in core.BREAKS
        self.outer_regions = [i for i in range(core.NS) if core.BREAKS[i] >= Q(7, 4)]
        self.nonnegative_proven = False
        self.roots = {}
        self.polynomials = {}
        f = list(map(Q, scalar['F']))
        fp = [Q(0)]*(2*len(f)-1)
        for k, value in enumerate(f):
            fp[2*k] = value*self.alpha**k
        matrices, denominator, _ = core.gram()
        for region in self.outer_regions:
            b = core.radial(region)
            start = region*core.NC
            m = np.zeros((core.NC, core.NC), dtype=object)
            for ell, matrix in enumerate(matrices):
                matrix = matrix[1:, 1:] if ell == 0 else matrix
                m += matrix[start:start+core.NC, start:start+core.NC]
            polynomial = b.T@m@b
            diagonal = [Q(0)]*(2*core.DEG+1)
            for i, j in itertools.product(range(core.NC), repeat=2):
                diagonal[i+j] += Q(int(polynomial[i, j]), denominator)
            interval = tuple(core.BREAKS[region:region+2])
            self.polynomials[region] = (fp, diagonal)
            self.roots[region] = tuple(self._integer_bernstein(p, interval)
                                       for p in [fp, diagonal])

    def establish_nonnegativity(self, freshly_computed_pair_reports):
        """Caller must supply results of its own completed exact pair run."""
        for region in self.outer_regions:
            rows = [r for r in freshly_computed_pair_reports if r['region'] == region]
            intervals = sorted((Q(r['left']), Q(r['right'])) for r in rows)
            assert rows and all(r['status'] == 'CERTIFIED' for r in rows)
            assert intervals[0][0] == self.core.BREAKS[region]
            assert intervals[-1][1] == self.core.BREAKS[region+1]
            assert all(a[1] == b[0] for a, b in zip(intervals, intervals[1:]))
        self.nonnegative_proven = True
        self.radial_bounds.cache_clear()

    def _integer_bernstein(self, polynomial, interval):
        denominator = math.lcm(*(v.denominator for v in polynomial))
        coefficients = np.array([int(v*denominator) for v in polynomial], dtype=object)
        coefficients, denominator = self.core.to_bernstein(coefficients, denominator, [interval])
        biggest = max(abs(int(v)) for v in coefficients.flat)
        precision = min(50, 58-(biggest.bit_length()-denominator.bit_length()+1))
        assert precision >= 0, 'Insufficient outward-integer dynamic range.'
        scale = 1 << precision
        low = np.array([int(v)*scale//denominator for v in coefficients.flat], dtype=np.int64)
        high = np.array([-((-int(v)*scale)//denominator) for v in coefficients.flat], dtype=np.int64)
        assert max(abs(low).max(), abs(high).max()) < 2**59
        return low, high, scale

    @lru_cache(maxsize=32768)
    def _dyadic(self, region, level, index):
        if level == 0:
            assert index == 0
            return self.roots[region]
        assert level > 0 and 0 <= index < 2**level
        parent = self._dyadic(region, level-1, index//2)
        side = index%2
        return tuple((self.core.split(low, 0, False)[side],
                      self.core.split(high, 0, True)[side], scale)
                     for low, high, scale in parent)

    @staticmethod
    @lru_cache(maxsize=32768)
    def _exp(t):
        return exponential_bounds(t)

    @lru_cache(maxsize=65536)
    def radial_bounds(self, radius, region):
        radius = tuple(map(Q, radius))
        lo, hi = radius
        left, right = self.core.BREAKS[region:region+2]
        assert left <= lo <= hi <= right
        if region not in self.outer_regions:
            return Q(0), Q(0)
        if lo == hi:
            values = []
            for polynomial in self.polynomials[region]:
                result = Q(0)
                for value in reversed(polynomial):
                    result = result*lo+value
                values.append((result, result))
        else:
            normalized_width = (hi-lo)/(right-left)
            index = (lo-left)/(hi-lo)
            denominator = normalized_width.denominator
            if (normalized_width.numerator == 1 and denominator&(denominator-1) == 0
                    and index.denominator == 1):
                arrays = self._dyadic(region, denominator.bit_length()-1, int(index))
            else:
                arrays = tuple(self._integer_bernstein(p, radius)
                               for p in self.polynomials[region])
            values = [(Q(int(low.min()), scale), Q(int(high.max()), scale))
                      for low, high, scale in arrays]
        exphi = self._exp(self.alpha*hi*hi)[1]
        explo = self._exp(self.alpha*lo*lo)[0]
        f = product_interval(values[0], (1/exphi, 1/explo))
        # D is a sum of squared radial functions with nonnegative weights.
        d = max(Q(0), values[1][0]), max(Q(0), values[1][1])
        result = -f[1]-d[1], -f[0]-d[0]
        if self.nonnegative_proven:
            assert result[1] >= 0, 'The supplied pair proof contradicts the radial enclosure.'
            result = max(Q(0), result[0]), result[1]
        return result


def correction_bounds(kind, bounds, regions, probes, budget):
    """Enclose the magnitude subtracted from the unlifted cyclic kernel."""
    core = budget.core
    assert kind in ['full', 'outside']
    if kind == 'full':
        assert len(bounds) == len(regions) == 3
        pairs = [(i, j, angle_bounds(bounds[i], bounds[j], bounds[k]))
                 for i, j, k in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]]
        radial = [budget.radial_bounds(tuple(bounds[i]), regions[i]) for i in range(3)]
    else:
        assert len(bounds) == 3 and len(regions) == 2
        a, b, z = bounds
        assert 1 <= a[0] <= a[1] <= 2 and 1 <= b[0] <= b[1] <= 2
        assert 0 <= z[0] <= z[1] <= 1
        # For r,s in[1,2], h=((r+s)^2-4)/(2rs) is coordinatewise
        # nondecreasing and nonnegative.  Thus t=-1+z*h has these extrema.
        h = lambda r, s: ((r+s)**2-core.BREAKS[-1]**2)/(2*r*s)
        cosine = Q(-1)+z[0]*h(a[0], b[0]), Q(-1)+z[1]*h(a[1], b[1])
        assert -1 <= cosine[0] <= cosine[1] <= 1
        pairs = [(0, 1, cosine)]
        radial = [budget.radial_bounds(tuple(bounds[i]), regions[i]) for i in range(2)]
    lower = upper = Q(0)
    for probe in probes:
        if probe['fraction'] == 0:
            continue
        for i, j, cosine in pairs:
            pi = probability_bounds(tuple(bounds[i]), cosine, probe)
            pj = probability_bounds(tuple(bounds[j]), cosine, probe)
            qi = product_interval(radial[j], pi)
            qj = product_interval(radial[i], pj)
            lower += probe['fraction']*(qi[0]+qj[0])/(2*probe['M'])
            upper += probe['fraction']*(qi[1]+qj[1])/(2*probe['M'])
    return lower, upper
