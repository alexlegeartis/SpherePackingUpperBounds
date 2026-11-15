"""Stage 1's subtracted corrections on an exterior branch, as exact symbolic values.

On an exterior branch ('outside', regions) the stage 1 lift helper
(record/current/origin_route/interval_combined_binary_six_graphs_correlated_small_cap_bounds.py,
function check) proves

    C(r, s, z) / (dd (2 r s)^18)  <=  tents + caps + gauge + allocations + graphs,

so stage 1's kernel there is  K = C/(dd (2rs)^18) - (sum of these corrections), and
-K = -C/(dd (2rs)^18) + corrections. This module writes each correction as a
symalg.Sym in (r, s, t), with z = 2 r s (1 + t) / ((r + s)^2 - (7/4)^2).

Sources (paths under stage1_v7/record/current/):
  caps        base/new_work/pointed/interval_certificate_general_closed.py: correction_bounds
              (weight/2 (w(s) P(r,t) + w(r) P(s,t))) and probability_bounds (P). The inner
              lambda c(r) = (r^2 + b)/(2 d r) is called directly. The branch structure of P
              is transcribed: sigma = -1: P = [c + t <= 0]; sigma = 1: P = [c - t <= 0];
              otherwise P = clip((1 - (c - sigma t)/sqrt((1 - sigma^2)(1 - t^2)))/2).
  radial w    base/new_work/pointed/interval_certificate.py: radial_bounds; its inner lambda
              bernstein(z) is called directly, z = (r - a)/(b - a) on the basis interval.
  tents       base/new_work/pointed/interval_certificate_tent_closed.py: correction_bounds
              (weight/2 (w(s) relu(value(r,s,z)) + w(r) relu(value(s,r,z)))); the inner
              lambda value of outside_tent_bounds is called directly.
  allocations base/new_work/pointed/independent_antipodal_allocation_exact.py and
              two_capacity_antipodal_allocation_exact.py: correction_bounds
              (weight (u(r) + u(s)) phi/2); the inner function value of
              squared_distance_bounds is called directly; u = (2r - 2)^4 on the basis
              interval [1, 3/2] (profile_bounds); phi = capacity_point_bounds profile.
  graphs      six_exterior_graph_route/*_bounds.py: weight times the closed radial and
              angular indicators.
  gauge       base/new_work/pointed/full_pair_count_gauge.py: correction_bounds
              ((w(r) psi(s) + w(s) psi(r))/2), FullPairCountBudget.radial_bounds
              (w = -q/(1 - psi), q = exp(-alpha r^2) fp(r) + dp(r) - reduction(r)),
              pair_interval_certificate_closed.py: diagonal_bounds (reduction).
singleton_checks evaluates each term at exact rational points and requires it to lie
in the enclosure that the stage 1 routine itself returns for the singleton box.
"""
from fractions import Fraction as Q
import sys
import types

import symalg as sa
from symalg import Sym

AA = Q(133, 250)


def _module(name, suffix):
    m = sys.modules[name]
    assert m.__file__.replace('\\', '/').endswith(suffix), (name, m.__file__)
    return m


def inner_function(outer, name, cells):
    codes = [c for c in outer.__code__.co_consts if isinstance(c, types.CodeType) and c.co_name == name]
    assert len(codes) == 1, (outer.__name__, name)
    code = codes[0]
    assert set(code.co_freevars) == set(cells), (code.co_freevars, cells)
    closure = tuple(types.CellType(cells[v]) for v in code.co_freevars)
    return types.FunctionType(code, outer.__globals__, name, None, closure)


class Stage1:
    def __init__(self, core, call, L):
        self.core = core
        self.call = call
        self.L = L
        self.R = core.BREAKS
        self.cutoff = core.BREAKS[-1]
        self.deg = core.DEG
        self.basis = call['basis_intervals']
        self.lift = L['lift_helper']
        assert self.lift.__file__.endswith('origin_route/interval_combined_binary_six_graphs_correlated_small_cap_bounds.py')
        # the modules whose functions the real check() calls
        mod = lambda f, suffix: _module(f.__module__, suffix)
        self.tc = mod(self.lift.tent_bounds, 'pointed/interval_certificate_tent_closed.py')
        self.alloc = mod(self.lift.allocation_bounds, 'pointed/two_capacity_antipodal_allocation_exact.py')
        self.gauge_mod = mod(self.lift.lift_bounds, 'pointed/full_pair_count_gauge.py')
        self.full_closed = mod(self.lift.cap_bounds, 'pointed/interval_certificate_full_closed.py')
        self.gc = mod(self.full_closed.exterior_correction_bounds, 'pointed/interval_certificate_general_closed.py')
        self.ic = mod(self.gc.radial_bounds, 'pointed/interval_certificate.py')
        self.graph_bounds = mod(self.lift.cyclic_decrement_lower, 'six_exterior_graph_route/graph_repair_bounds.py')
        self.outer_graph = mod(self.lift.outer_graph_lower, 'six_exterior_graph_route/outer_exterior_graph_bounds.py')
        self.narrow_graph = mod(self.lift.narrow_graph_lower, 'six_exterior_graph_route/narrow_exterior_graph_bounds.py')
        self.moderate_graph = mod(self.lift.moderate_graph_lower, 'six_exterior_graph_route/moderate_exterior_graph_bounds.py')
        base = call['budget'].base
        diag = type(base).radial_bounds.__wrapped__.__globals__['diagonal_bounds']
        assert diag.__globals__['probability_bounds'] is self.gc.probability_bounds
        assert diag.__globals__['radial_bounds'] is self.ic.radial_bounds
        assert self.gc.correction_bounds.__globals__['radial_bounds'] is self.ic.radial_bounds
        assert self.tc.correction_bounds.__globals__['radial_bounds'] is self.ic.radial_bounds
        assert self.alloc.profile_bounds is sys.modules['independent_antipodal_allocation_exact'].profile_bounds
        self.bern = {j: inner_function(self.ic.radial_bounds, '<lambda>', dict(degree=self.deg, j=j))
                     for j in range(self.deg + 1)}
        self.sqd = inner_function(self.alloc.squared_distance_bounds, 'value', dict(cutoff=self.cutoff))

    # ---------------------------------------------------------- pieces
    def zvar(self, r, s, t):
        return (2 * r * s * (1 + t)) / ((r + s) ** 2 - self.cutoff ** 2)

    def radial(self, rho, region, w):
        a, b = self.basis[region]
        zz = (rho - a) / (b - a)
        total = Sym.const(0)
        for j, weight in enumerate(w[region * (self.deg + 1):(region + 1) * (self.deg + 1)]):
            if not weight:
                continue
            total = total + self.bern[j](zz) * weight
        return total

    def prob(self, rho, t, cap):
        d, a2, sigma = cap['d'], cap['a2'], cap['sigma']
        b = d * d - a2
        c = inner_function(self.gc.probability_bounds, '<lambda>', dict(b=b, d=d))(rho)
        if sigma == -1:
            return sa.ind_le0(c + t)
        if sigma == 1:
            return sa.ind_le0(c - t)
        if t == -1:
            return sa.ind_le0(c + sigma)
        D = ((1 - sigma * sigma) * (1 - t * t)).sqrt()
        return sa.clip01((1 - (c - sigma * t) / D) / 2)

    def caps(self, r, s, t, regions):
        total = Sym.const(0)
        for cap in self.call['caps']:
            wa = self.radial(r, regions[0], cap['w'])
            wb = self.radial(s, regions[1], cap['w'])
            if wa.is_zero() and wb.is_zero():
                continue
            term = wb * self.prob(r, t, cap) + wa * self.prob(s, t, cap)
            total = total + term * (cap['weight'] / 2)
        return total

    def tents(self, r, s, t, z, regions):
        total = Sym.const(0)
        for tent in self.call['tents']:
            value = inner_function(self.tc.outside_tent_bounds, '<lambda>', dict(cutoff=self.cutoff, d=tent['d']))
            wa = self.radial(r, regions[0], tent['w'])
            wb = self.radial(s, regions[1], tent['w'])
            term = wb * sa.relu(value(r, s, z)) + wa * sa.relu(value(s, r, z))
            total = total + term * (tent['weight'] / 2)
        return total

    def profile(self, rho, region):
        interval = tuple(map(Q, self.basis[region]))
        if interval == self.alloc.OUTER_INTERVAL:
            return Sym.const(0)
        assert interval == self.alloc.INNER_INTERVAL
        return (2 * rho - 2) ** 4

    def allocations(self, r, s, t, z, regions):
        total = Sym.const(0)
        u = self.profile(r, regions[0]) + self.profile(s, regions[1])
        for a in self.call['allocations']:
            assert a['capacity_profile'] in ('classical_improvement', 'two_level'), a['capacity_profile']
            phi = sa.capacity(a['capacity_profile'], self.sqd(r, s, z))
            total = total + (u * phi) * (Q(a['weight']) / 2)
        return total

    def _contained(self, region, interval):
        lo, hi = self.R[region], self.R[region + 1]
        a, b = interval
        if a <= lo and hi <= b:
            return True
        if hi < a or lo >= b:   # (lo, hi] misses [a, b]
            return False
        raise AssertionError('region straddles a graph radius')

    def graphs(self, r, s, t, regions):
        gw = self.call['graph_weights']
        assert gw[0] == 0 and gw[1] == 0  # graphs A and B carry no weight
        gb = self.graph_bounds
        specs = [(gw[2], gb.ANTIPODAL_RADIAL, gb.ANTIPODAL_ANGULAR),
                 (self.call['outer_graph_weight'], self.outer_graph.RADIAL, self.outer_graph.ANGULAR),
                 (self.call['narrow_graph_weight'], (Q(33, 28), Q(7, 4)), (Q(-1), Q(-183, 200))),
                 (self.call['moderate_graph_weight'], self.moderate_graph.RADIAL, self.moderate_graph.ANGULAR)]
        total = Sym.const(0)
        for w, radial, angular in specs:
            if w == 0:
                continue
            if self._contained(regions[0], radial) and self._contained(regions[1], radial):
                total = total + sa.band(t, *angular) * w
        return total

    def q(self, rho, region):
        base = self.call['budget'].base
        assert not base.cuts
        assert base.caps == tuple(self.call['caps']) or list(base.caps) == list(self.call['caps'])
        fp, dp = base.polynomials[region]
        fpoly = sum((rho ** k * Q(c) for k, c in enumerate(fp) if c), Sym.const(0))
        dpoly = sum((rho ** k * Q(c) for k, c in enumerate(dp) if c), Sym.const(0))
        out = (-(rho * rho) * base.alpha).exp() * fpoly + dpoly
        lo, hi = self.R[region], self.R[region + 1]
        for cap in base.caps:
            p = self.gc.probability_bounds((lo, hi), (Q(1), Q(1)), cap)
            if 'pivot_lower_cutoff' in cap:
                cut = cap['pivot_lower_cutoff']
                if hi <= cut:
                    p = (Q(cap['M']), Q(cap['M']))
                else:
                    assert lo > cut, 'region straddles a pivot cutoff'
            assert p[0] == p[1], 'self probability undecided on a region'
            if p[1] == 0:
                continue
            out = out - self.radial(rho, region, cap['w']) * (cap['weight'] * p[0])
        return out

    def gauge(self, r, s, regions):
        flags = self.call['budget'].base.gauge['flags']
        fr, fs = Q(flags[regions[0]]), Q(flags[regions[1]])
        if not fr and not fs:
            return Sym.const(0)
        wr = self.q(r, regions[0]) * (-1 / (1 - fr))
        ws = self.q(s, regions[1]) * (-1 / (1 - fs))
        return (wr * fs + ws * fr) / 2

    def cell(self, i, k, endpoint=False):
        Sym.BOX = ((self.R[i], self.R[i + 1]), (self.R[k], self.R[k + 1]), (Q(-1), -AA))
        r, s, t = sa.var('r'), sa.var('s'), sa.var('t')
        if endpoint:
            t = Sym.const(-1)
        z = self.zvar(r, s, t)
        comp = dict(caps=self.caps(r, s, t, (i, k)), tents=self.tents(r, s, t, z, (i, k)),
                    allocations=self.allocations(r, s, t, z, (i, k)), graphs=self.graphs(r, s, t, (i, k)),
                    gauge=self.gauge(r, s, (i, k)))
        return comp


def build_all(core, call, L, nreg):
    st = Stage1(core, call, L)
    cells, components, endpoints = {}, {}, {}
    for i in range(nreg):
        for k in range(nreg):
            comp = st.cell(i, k)
            components[(i, k)] = comp
            total = Sym.const(0)
            for v in comp.values():
                total = total + v
            cells[(i, k)] = total
            endpoints[(i, k)] = sum(st.cell(i, k, endpoint=True).values(), Sym.const(0))
    return dict(cells=cells, components=components, endpoints=endpoints, helper=st)


# ------------------------------------------------------------------ singleton checks
def _points(R, i, k):
    out = []
    for u, v, tt in [(Q(1, 3), Q(2, 3), Q(-9, 10)), (Q(1), Q(1, 5), Q(-7, 10)),
                     (Q(1, 7), Q(1), Q(-6, 10)), (Q(1, 2), Q(1, 2), Q(-77, 100)),
                     (Q(3, 4), Q(1, 9), Q(-99, 100))]:
        r = R[i] + (R[i + 1] - R[i]) * u
        s = R[k] + (R[k + 1] - R[k]) * v
        z = (1 + tt) * 2 * r * s / ((r + s) ** 2 - Q(49, 16))
        assert 0 <= z <= 1
        out.append((r, s, z, tt))
    return out


def singleton_checks(core, call, L, nreg, pos):
    from flint import arb
    st = pos['helper']
    R = core.BREAKS
    budget = call['budget']
    count, worst = 0, []
    for (i, k), comp in pos['components'].items():
        regions = (i, k)
        if i > k:
            continue  # stage 1 routines take the branch's own region order
        for r, s, z, t in _points(R, i, k):
            b = [(r, r), (s, s), (z, z)]
            ref = dict(
                caps=st.lift.cap_bounds('outside', b, regions, call['caps'], core, st.basis),
                tents=st.tc.correction_bounds('outside', b, regions, call['tents'], core, st.basis),
                allocations=st.alloc.correction_bounds('outside', b, regions, call['allocations'], core, st.basis),
                gauge=st.gauge_mod.correction_bounds('outside', b, regions, [], budget))
            g = (st.graph_bounds.cyclic_decrement_lower('outside', b, call['graph_weights'])
                 + st.outer_graph.cyclic_decrement_lower('outside', b, call['outer_graph_weight'])
                 + st.narrow_graph.cyclic_decrement_lower('outside', b, call['narrow_graph_weight'])
                 + st.moderate_graph.cyclic_decrement_lower('outside', b, call['moderate_graph_weight']))
            ref['graphs'] = (g, g)
            box = ((r, r), (s, s), (t, t))
            for name, (lo, hi) in ref.items():
                v = sa.eval_arb(comp[name], box)
                enc = sa._ball(lo, hi)
                assert v.is_finite() and v.overlaps(enc), (name, (i, k), str(r), str(s), str(t), v, float(lo), float(hi))
                worst.append(float(hi - lo))
                count += 1
    return dict(points=count, max_stage1_enclosure_width=max(worst))
