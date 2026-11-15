"""Exact lower bounds reusing freshly constructed complete-pair majorants.

This is an optional enclosure wrapper, not a verifier or a certificate.
The calling verifier must construct each supplied P from the SAME actual
pair polynomial that it freshly proves, with the contract

    P(r) >= exp(alpha*r*r) * (q(r) + J(r)).

Here q already includes every signed energy and active count cut, while
J is the nonnegative pointed self reduction.  In particular the wrapper
MUST NOT add any count cut to its new formula.  Original upper bounds,
fresh pair gates, and all cyclic coverage and rejection gates are retained.

Only direct coefficient arrays and exact intervals are accepted.  There
is no path or saved-status loader.  Freshness is a calling-verifier
obligation, just as it is for FullPairCountBudget.establish_nonnegativity.
"""
from fractions import Fraction as Q
from functools import lru_cache
from numbers import Integral

from pair_slack_lift import ScalarPairBudget
from pair_interval_certificate_closed import diagonal_bounds


def _rational(value):
    assert not isinstance(value, bool)
    assert isinstance(value, (str, Integral, Q)), 'Use exact rational inputs.'
    return Q(str(value)) if isinstance(value, Integral) else Q(value)


def _integer(value):
    assert isinstance(value, Integral) and not isinstance(value, bool)
    return int(value)


def _key(record):
    region = _integer(record['region'])
    left, right = _rational(record['left']), _rational(record['right'])
    assert left < right
    return region, left, right


class CancellationPreservingPairBudget:
    """Refine only the lower enclosure of a freshly proved full pair budget.

    majorants is a list of dictionaries with exactly these keys:
      region, left, right, coefficients, denominator.
    coefficients / denominator is the physical-r power polynomial P.
    Each list entry is made in the calling verifier's current pair loop,
    after its corresponding fresh report returned CERTIFIED.

    fresh_pair_reports must be the reports from that same invocation.  All
    exact regions, all active count-cut partitions, and no duplicate or
    missing subinterval are required.  Both closed extensions at a cut
    are included in a crossing or singleton query.  Their minimum is safe
    even when the actual count indicator selects only the closed left side.
    """

    def __init__(self, base, majorants, fresh_pair_reports):
        assert base.nonnegative_proven is True, 'Fresh complete pair proofs are required.'
        self.base = base
        self.core = base.core
        self.alpha = Q(base.alpha)
        self.gauge = base.gauge
        assert self.alpha > 0
        assert isinstance(majorants, (list, tuple)) and majorants
        assert isinstance(fresh_pair_reports, (list, tuple)) and fresh_pair_reports
        assert set(base.polynomials) == set(range(self.core.NS))
        assert len(self.gauge['flags']) == self.core.NS
        assert all(Q(0) <= Q(flag) < 1 for flag in self.gauge['flags'])

        expected = set()
        for region in range(self.core.NS):
            left, right = map(Q, self.core.BREAKS[region:region+2])
            assert left < right
            boundaries = sorted({left, right} | {
                Q(cut) for cut, _ in base.cuts if left < Q(cut) < right})
            expected.update((region, a, b) for a, b in zip(boundaries, boundaries[1:]))
        reports = {}
        for report in fresh_pair_reports:
            key = _key(report)
            assert key not in reports, 'Duplicate fresh pair report.'
            assert report['status'] == 'CERTIFIED'
            reports[key] = report
        assert set(reports) == expected, 'Fresh pair reports must cover exactly every piece.'

        self.powers, self.roots, self.pieces = {}, {}, {}
        allowed_keys = {'region', 'left', 'right', 'coefficients', 'denominator'}
        for source in majorants:
            assert isinstance(source, dict) and set(source) == allowed_keys
            key = _key(source)
            assert key in expected and key not in self.powers
            denominator = _integer(source['denominator'])
            assert denominator > 0
            raw = source['coefficients']
            # Iteration copies every scalar, so subsequent caller mutation
            # cannot change the certified polynomial or its cached bounds.
            if hasattr(raw, 'ndim'):
                assert raw.ndim == 1
            coefficients = tuple(_integer(value) for value in raw)
            assert coefficients
            polynomial = tuple(Q(value, denominator) for value in coefficients)
            self.powers[key] = polynomial
            self.roots[key] = ScalarPairBudget._integer_bernstein(self, polynomial, key[1:])
        assert set(self.powers) == expected, 'Majorants must cover exactly the fresh pair partition.'
        for region in range(self.core.NS):
            self.pieces[region] = tuple(sorted(key for key in expected if key[0] == region))
        self.contract = 'FRESH_COMPLETE_PAIR_MAJORANT_WITH_COUNT_CUT_ALREADY_INCLUDED'

    @property
    def nonnegative_proven(self):
        return self.base.nonnegative_proven

    @lru_cache(maxsize=32768)
    def _dyadic_majorant(self, key, level, index):
        assert type(level) is int and level >= 0
        assert type(index) is int and 0 <= index < 2**level
        if level == 0:
            return self.roots[key]
        low, high, scale = self._dyadic_majorant(key, level-1, index//2)
        side = index % 2
        return (self.core.split(low, 0, False)[side],
                self.core.split(high, 0, True)[side], scale)

    def _majorant_upper(self, key, radius):
        left, right = key[1:]
        lo, hi = radius
        assert left <= lo <= hi <= right
        if lo == hi:
            value = Q(0)
            for coefficient in reversed(self.powers[key]):
                value = value*lo+coefficient
            return value
        normalized_width = (hi-lo)/(right-left)
        index = (lo-left)/(hi-lo)
        denominator = normalized_width.denominator
        if (normalized_width.numerator == 1 and denominator & (denominator-1) == 0
                and index.denominator == 1):
            _, high, scale = self._dyadic_majorant(key, denominator.bit_length()-1, int(index))
        else:
            _, high, scale = ScalarPairBudget._integer_bernstein(self, self.powers[key], radius)
        return Q(int(high.max()), scale)

    @lru_cache(maxsize=65536)
    def _refined_bounds(self, radius, region):
        lo, hi = radius
        left, right = map(Q, self.core.BREAKS[region:region+2])
        assert left <= lo <= hi <= right
        old_lower, old_upper = self.base.radial_bounds(radius, region)
        assert Q(0) <= old_lower <= old_upper
        self_factor = 1-Q(self.gauge['flags'][region])
        assert self_factor > 0
        new_lowers = []
        for key in self.pieces[region]:
            a, b = max(lo, key[1]), min(hi, key[2])
            if a > b:
                continue
            interval = (a, b)
            pmax = self._majorant_upper(key, interval)
            e_lower = self.base._exp(self.alpha*a*a)[0]
            e_upper = self.base._exp(self.alpha*b*b)[1]
            assert 0 < e_lower <= e_upper
            # For Pmax < 0 the upper bound divides by the LARGER
            # positive exponential endpoint; using e_lower is unsound.
            quotient_upper = pmax/(e_lower if pmax >= 0 else e_upper)
            j_lower, _ = diagonal_bounds(interval, region, self.base.caps,
                                         self.core, self.base.basis_intervals)
            assert j_lower >= 0
            new_lowers.append((j_lower-quotient_upper)/self_factor)
        assert new_lowers, 'The exact closed pair partition must cover the query.'
        lower = max(old_lower, min(new_lowers))
        assert lower <= old_upper, 'The new lower bound contradicts the original enclosure.'
        return lower, old_upper

    def radial_bounds(self, radius, region):
        # This check precedes cache lookup, so revoking the fresh proof
        # gate cannot return a previously cached nonnegative enclosure.
        assert self.base.nonnegative_proven is True
        assert self.base.core is self.core and Q(self.base.alpha) == self.alpha
        assert _integer(region) == region and 0 <= region < self.core.NS
        radius = tuple(_rational(value) for value in radius)
        assert len(radius) == 2
        return self._refined_bounds(radius, int(region))

    def clear_cache(self):
        self._dyadic_majorant.cache_clear()
        self._refined_bounds.cache_clear()
