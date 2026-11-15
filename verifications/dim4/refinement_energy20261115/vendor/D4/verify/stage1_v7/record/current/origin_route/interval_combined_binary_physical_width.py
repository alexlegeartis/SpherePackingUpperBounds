"""Exact cyclic verification with a previously proved full pair/count gauge.

This separate checker extends the audited closed-boundary cap/tent/count
helper by subtracting the independently validated antipodal correction.
The allocation is also present in the complete pair budget, via the
mandatory linear-origin identity in its calling global verifier.
The positive Gaussian-polynomial budget must first be established by the
calling verifier's fresh exact unlifted pair checks.  Accepted historical
helper sources are not modified.
"""
from fractions import Fraction as Q
import itertools
from functools import lru_cache
import time
import hashlib
from pathlib import Path
import sys
import numpy as np

_GEOMETRY = Path(__file__).resolve().parent.parent/'geometry_route'
_BINARY_PINS = {
    'fixed_binary_bernstein_v3.py': 'dee3a3fafa6ba7cb811bae7d2568b8633ed68eed9702155132cac79f70075f07',
    'fixed_binary_bernstein_two_product_v3.py': '39d585defddfbee9935daa25743dafb27c01fd8dda9f642b17da9ef657a27176',
    'binary_bernstein_signed_exponent_v3.py': '61da8cac898370a5a9daf1a599cba79109a8d08e53a7ddf665e4e17bf0fa37c6',
}
for _name, _digest in _BINARY_PINS.items():
    assert hashlib.sha256((_GEOMETRY/_name).read_bytes()).hexdigest() == _digest
sys.path.insert(0, str(_GEOMETRY))
import fixed_binary_bernstein_two_product_v3 as _binary
import binary_bernstein_signed_exponent_v3 as _handoff
assert Path(_binary.__file__).resolve() == _GEOMETRY/'fixed_binary_bernstein_two_product_v3.py'
assert Path(_handoff.__file__).resolve() == _GEOMETRY/'binary_bernstein_signed_exponent_v3.py'
_COEFFICIENT_PRECISION = 192
from signed_spherical_exact import build_signed
from interval_certificate_tent_closed import correction_bounds as tent_bounds
from interval_certificate_full_closed import correction_bounds as cap_bounds
from full_pair_count_gauge import correction_bounds as lift_bounds
from two_capacity_antipodal_allocation_exact import correction_bounds as allocation_bounds


def check(core, kind, regions, eta, tents, basis_intervals, budget, probes,
          node_limit=180000, depth_limit=66, caps=(), allocations=(), spherical_rays=()):
    """Prove one closed cyclic branch, with no inference about other branches."""
    started = time.time()
    assert kind in ['full', 'outside'] and budget.core is core
    assert budget.nonnegative_proven, 'Fresh exact base-pair proofs are required first.'
    C, denominator, domain, signed_identity = build_signed(core, kind, regions, eta, spherical_rays)
    # Every discarded nonnegative correction has an EXACT zero upper
    # bound on the entire closed branch, including every descendant box.
    retained_caps, zero_cap_indices = [], []
    for index, cap in enumerate(caps):
        lower, upper = cap_bounds(kind, domain, regions, [cap], core, basis_intervals)
        assert 0 <= lower <= upper
        if upper == 0:
            zero_cap_indices.append(index)
        else:
            retained_caps.append(cap)
    cap_screen = dict(status='EXACT_WHOLE_BRANCH_ZERO_CAP_SCREEN',
                      original_count=len(caps), retained_count=len(retained_caps),
                      zero_cap_indices=zero_cap_indices)
    caps = tuple(retained_caps)
    # Every Bernstein coefficient is enclosed by outward integer
    # arithmetic. No sampled coefficient or sign estimate enters here.
    coefficient_low, coefficient_high, coefficient_scale = _binary.to_bernstein_interval(
        C, denominator, domain, precision=_COEFFICIENT_PRECISION)
    del C
    try:
        low, high, scale, precision = _handoff.to_int64_signed_exponent(
            coefficient_low, coefficient_high, coefficient_scale)
    except OverflowError:
        return dict(status='UNKNOWN', reason='coefficient dynamic range')
    assert type(scale) is Q and scale > 0
    assert low.dtype == high.dtype == np.dtype(np.int64)
    assert max(abs(low).max(), abs(high).max()) < 2**59
    enclosure_record = dict(
        method='EXACT_OUTWARD_TWO_PRODUCT_BERNSTEIN_V3',
        all_coefficients_enclosed=True, coefficient_precision=_COEFFICIENT_PRECISION,
        output_exponent=precision, exact_positive_output_scale=str(scale),
        source_sha256=dict(_BINARY_PINS))
    del coefficient_low, coefficient_high
    angular_degree = core.gram()[2]

    def factor(point):
        value = 2*point[0]*point[1]
        if kind == 'full':
            value *= point[2]
        return value**angular_degree

    @lru_cache(maxsize=8192)
    def all_corrections(bounds):
        tl, tu = tent_bounds(kind, bounds, regions, tents, core, basis_intervals)
        cl, cu = cap_bounds(kind, bounds, regions, caps, core, basis_intervals) if caps else (Q(0), Q(0))
        ll, lu = lift_bounds(kind, bounds, regions, probes, budget)
        al, au = allocation_bounds(kind, bounds, regions, allocations, core,
                                   basis_intervals)
        assert min(tl, cl, ll, al) >= 0
        return tl+cl+ll+al, tu+cu+lu+au, tl+cl+al, tl+cl+ll

    stack = [(low, high, ((0, 1),)*3, (0,)*3)]
    nodes = leaves = correction_leaves = lift_dependent_leaves = maximum_depth = 0
    antipodal_dependent_leaves = 0
    failures = []
    unknown = False
    while stack:
        lo, hi, box, levels = stack.pop()
        nodes += 1
        maximum_depth = max(maximum_depth, sum(levels))
        maximum = int(hi.max())
        if maximum <= 0:
            leaves += 1
            continue
        bounds = [(left+(right-left)*Q(x, 2**level), left+(right-left)*Q(y, 2**level))
                  for (left, right), (x, y), level in zip(domain, box, levels)]
        correction_lower, _, unlifted_lower, without_antipodal = all_corrections(tuple(bounds))
        factor_lower = factor([b[0] for b in bounds])
        if Q(maximum, scale) <= factor_lower*correction_lower:
            leaves += 1
            correction_leaves += 1
            lift_dependent_leaves += int(Q(maximum, scale) > factor_lower*unlifted_lower)
            antipodal_dependent_leaves += int(Q(maximum, scale) > factor_lower*without_antipodal)
            continue
        failure = None
        for corner in itertools.product((0, 1), repeat=3):
            index = tuple(-1 if side else 0 for side in corner)
            base_lower = int(lo[index])
            if base_lower <= 0:
                continue
            point = [bounds[i][side] for i, side in enumerate(corner)]
            _, correction_upper, _, _ = all_corrections(tuple((p, p) for p in point))
            if Q(base_lower, scale) > factor(point)*correction_upper:
                failure = dict(point=list(map(str, point)), base_integer_lower=base_lower,
                               correction_upper=str(correction_upper))
                break
        if failure:
            failures.append(failure)
            break
        if nodes >= node_limit or sum(levels) >= depth_limit:
            unknown = True
            break
        axis = max(range(3), key=lambda j: (Q(domain[j][1]) - Q(domain[j][0])) / (1 << levels[j]))
        ll, lr = core.split(lo, axis, False)
        ul, ur = core.split(hi, axis, True)
        nextlevels = list(levels)
        nextlevels[axis] += 1
        leftbox, rightbox = list(box), list(box)
        a, b = box[axis]
        leftbox[axis], rightbox[axis] = (2*a, a+b), (a+b, 2*b)
        stack.extend([(lr, ur, tuple(rightbox), tuple(nextlevels)),
                      (ll, ul, tuple(leftbox), tuple(nextlevels))])
    status = 'VIOLATION' if failures else ('UNKNOWN' if unknown or stack else 'CERTIFIED')
    record = dict(status=status, kind=kind, regions=list(regions), nodes=nodes, leaves=leaves,
                  signed_spherical_identity=signed_identity,
                  exact_cap_screen=cap_screen,
                  bernstein_enclosure=enclosure_record,
                  correction_cache=all_corrections.cache_info()._asdict(),
                  correction_leaves=correction_leaves, lift_dependent_leaves=lift_dependent_leaves,
                  antipodal_dependent_leaves=antipodal_dependent_leaves,
                  maxdepth=maximum_depth, binary_scale=precision, failures=failures,
                  seconds=time.time()-started,
                  pointed_helper='EXACT_CAP_TENT_TWO_CAPACITY_ANTIPODAL_AND_VERIFIED_FULL_PAIR_COUNT_GAUGE')
    print('EXACT BINARY COMBINED COUNT GAUGE CYCLIC', record, flush=True)
    return record
