"""Exact bounds for converting a full pair residual into a count energy.

Let q(r) be the COMPLETE unlifted pair expression on the kernel support.
For a packing neighborhood C suppose 0 <= psi(r) < 1 and sum_C psi <= 1.
Once q <= 0 has been independently proved, put w=-q/(1-psi). Then

    (sum_C w) (1-sum_C psi) >= 0.

The added origin constant is zero, its pair contribution is w(1-psi)=-q,
and its off-diagonal kernel is -(w(r)psi(s)+w(s)psi(r))/2. Thus the actual
lifted pair is identically zero and the cyclic correction is nonpositive.

This module supports the initial gauge psi=1_[1,R]/M. It encloses the whole
q, including the restored origin row, all count-polynomial terms, and the
negative pointed self diagonal. No sign is assumed for its polynomial
part. The nonnegativity clamp is enabled only by fresh complete pair
proofs, not by any status or coefficient supplied in a candidate JSON.
"""
from fractions import Fraction as Q
from functools import lru_cache
import itertools

from pair_slack_lift import ScalarPairBudget, product_interval
from pair_interval_certificate_closed import diagonal_bounds


def validate_gauge(record, core, proved_count_radii, proved_count_bounds):
    """Tie an indicator gauge to a count theorem replayed by the caller."""
    assert record['gauge'] == 'indicator_over_count_bound'
    assert record['weight_rule'] == 'negative_of_unlifted_pair_divided_by_one_minus_gauge'
    radius, support = Q(record['radius']), Q(record['support_radius'])
    maximum = int(record['count_bound'])
    assert Q(record['count_bound']) == maximum and maximum >= 2
    assert 1 < radius < support == core.BREAKS[-1]
    assert radius in core.BREAKS
    assert Q(record.get('self_factor_inside', Q(maximum-1, maximum))) == Q(maximum-1, maximum)
    assert Q(record.get('self_factor_outside', 1)) == 1
    assert len(proved_count_radii) == len(proved_count_bounds)
    assert any(radius <= Q(r) and int(m) <= maximum
               for r, m in zip(proved_count_radii, proved_count_bounds)), \
        'The gauge needs a matching independently verified local count bound.'
    flags = [Q(int(core.BREAKS[i+1] <= radius), maximum) for i in range(core.NS)]
    assert all(core.BREAKS[i+1] <= radius or core.BREAKS[i] >= radius
               for i in range(core.NS))
    return dict(radius=radius, count_bound=maximum, support_radius=support,
                flags=flags, origin_constant=Q(0), actual_pair_multiplier=Q(0))


class FullPairCountBudget(ScalarPairBudget):
    """Enclose -q/(1-psi) using exact univariate polynomial and cap bounds.

    pair_polynomials[region] is the restored polynomial part of q. It
    includes 2h, the PSD diagonal, nested and conditional count terms, and
    the linear parts of all pointed and tent energies. The pointed self
    reduction and optional count cuts are subtracted separately from q.
    The calling checker constructs these polynomials from the same exact
    identities used by its mandatory base-pair proofs.
    """
    def __init__(self, core, scalar, pair_polynomials, caps, basis_intervals,
                 gauge, cuts=()):
        self.core = core
        self.alpha = Q(scalar['alpha'])
        assert self.alpha > 0
        self.outer_regions = list(range(core.NS))  # inherited proof gate: ALL regions
        self.nonnegative_proven = False
        self.gauge = gauge
        self.caps = tuple(caps)
        self.basis_intervals = tuple(basis_intervals)
        self.cuts = tuple((Q(row['R']), Q(row['weight'])) for row in cuts)
        assert all(weight >= 0 for _, weight in self.cuts)
        assert set(pair_polynomials) == set(self.outer_regions)
        f = list(map(Q, scalar['F']))
        fp = [Q(0)]*(2*len(f)-1)
        for k, value in enumerate(f):
            fp[2*k] = value*self.alpha**k
        self.polynomials, self.roots = {}, {}
        for region in self.outer_regions:
            dp = list(map(Q, pair_polynomials[region]))
            interval = tuple(core.BREAKS[region:region+2])
            self.polynomials[region] = (fp, dp)
            self.roots[region] = tuple(self._integer_bernstein(p, interval) for p in [fp, dp])

    @lru_cache(maxsize=65536)
    def radial_bounds(self, radius, region):
        radius = tuple(map(Q, radius))
        lo, hi = radius
        left, right = self.core.BREAKS[region:region+2]
        assert left <= lo <= hi <= right
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
        # Unlike a PSD diagonal, the complete restored polynomial has no
        # predetermined sign. Clamping it at zero would be incorrect.
        d = values[1]
        reduction = diagonal_bounds(radius, region, self.caps, self.core,
                                    self.basis_intervals)
        cut_lower = sum((weight for cut, weight in self.cuts if hi <= cut), Q(0))
        cut_upper = sum((weight for cut, weight in self.cuts if lo <= cut), Q(0))
        residual = (-f[1]-d[1]+reduction[0]+cut_lower,
                    -f[0]-d[0]+reduction[1]+cut_upper)
        self_factor = 1-self.gauge['flags'][region]
        assert self_factor > 0
        result = residual[0]/self_factor, residual[1]/self_factor
        if self.nonnegative_proven:
            assert result[1] >= 0, 'Fresh base-pair proof contradicts the enclosure.'
            result = max(Q(0), result[0]), result[1]
        return result


def correction_bounds(kind, bounds, regions, unused_probes, budget):
    """Enclose the positive amount subtracted from a cyclic expression."""
    assert not unused_probes
    assert kind in ['full', 'outside'] and len(bounds) == 3
    count = 3 if kind == 'full' else 2
    assert len(regions) == count
    flags = [budget.gauge['flags'][region] for region in regions]
    if not any(flags):
        return Q(0), Q(0)
    weights = [budget.radial_bounds(tuple(bounds[i]), region)
               for i, region in enumerate(regions)]
    lower = upper = Q(0)
    for i, j in itertools.combinations(range(count), 2):
        lower += (weights[i][0]*flags[j]+weights[j][0]*flags[i])/2
        upper += (weights[i][1]*flags[j]+weights[j][1]*flags[i])/2
    return lower, upper
