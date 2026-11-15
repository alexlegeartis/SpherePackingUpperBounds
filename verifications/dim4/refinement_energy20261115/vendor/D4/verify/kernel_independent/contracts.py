"""Exact component-data contracts and closed-boundary regression checks.

The contracts compare the objects actually passed to stage 1's checker.
They are not a second interval sign proof. See COMPONENT_IDENTITIES.md for
the mathematical meaning of the formulas and their boundary conventions.
"""
from fractions import Fraction as Q
import itertools
import sys


def compare(K, state):
    equations = 0

    def eq(left, right, name):
        nonlocal equations
        if left != right:
            raise ValueError('component contract differs: '+name)
        equations += 1

    core = state['core']
    lift = state['lift_helper']
    eq(list(map(Q, core.BREAKS)), K.breaks, 'radial breaks')
    eq(list(map(tuple, state['basis_intervals'])), K.basis, 'basis intervals')
    eq(core.DEG, K.deg, 'radial degree')
    eq(state['angular_degree'], K.L, 'angular degree')
    for p in range(K.np):
        ref = core.radial(p)
        for j in range(K.nc):
            eq(list(map(Q, ref[j])), K.phi_power[p][j], 'Bernstein polynomial %s/%s' % (p, j))

    for name, mine, other, fields in [
        ('cap', K.caps, state['caps'], ('d', 'a2', 'sigma', 'M', 'weight', 'w')),
        ('tent', K.tents, state['tents'], ('d', 'capacity', 'weight', 'w')),
        ('allocation', K.allocs, state['allocations'], ('weight', 'w')),
    ]:
        eq(len(mine), len(other), name+' count')
        for i, (a, b) in enumerate(zip(mine, other)):
            for field in fields:
                eq(a[field], b[field], '%s[%s].%s' % (name, i, field))
            if name == 'cap':
                eq(a['cutoff'], b.get('pivot_lower_cutoff'), 'cap cutoff')
                if a['cutoff'] is not None:
                    # With pivot radius <= C and another packing point at
                    # distance >= 1, the closed probe ball is not reached.
                    assert a['sigma'] == 1 and a['d']**2 <= a['a2'] < Q(1,2)
                    assert 1 <= a['cutoff']
                    assert a['cutoff'] < (1+a['d']**2-a['a2'])/a['d']
            elif name == 'allocation':
                eq(a['profile'], b['capacity_profile'], 'allocation profile')
                if a['profile'] == 'two_level':
                    eq(a['inner'], Q(b['exact_inner_radius']), 'strict inner radius')
                    eq(a['outer'], Q(b['exact_outer_radius']), 'closed outer radius')

    base = state['budget'].base
    eq(list(base.caps), state['caps'], 'gauge cap reductions')
    eq(tuple(base.basis_intervals), tuple(K.basis), 'gauge radial basis')
    eq(base.alpha, K.alpha, 'scalar alpha')
    eq(base.cuts, (), 'no additional count cuts')
    eq(base.gauge['radius'], K.gauge['R'], 'gauge radius')
    eq(base.gauge['count_bound'], K.gauge['bound'], 'gauge count bound')
    eq(list(base.gauge['flags']), [K.gauge_psi(p) for p in range(K.np)], 'gauge flags')
    eq(state['lift_probes'], [], 'no additional gauge probes')
    weights = list(state['graph_weights'])+[state['outer_graph_weight'],
               state['narrow_graph_weight'], state['moderate_graph_weight']]
    eq(weights, [g['delta'] for g in K.graphs], 'graph weights')

    # Both allocation branch boundaries, including the strict inner and
    # closed outer convention, plus the classical polynomial/root join.
    allocation_module = sys.modules[lift.allocation_bounds.__module__]
    allocation_boundary_checks = 0
    for al in K.allocs:
        cuts = [Q(0), Q(1, 4), Q(1, 2)]
        if al['profile'] == 'two_level':
            cuts += [al['inner']**2, al['outer']**2]
        cuts = sorted(set(cuts))
        points = cuts+[(a+b)/2 for a,b in zip(cuts,cuts[1:])]+[Q(3, 4)]
        for x in points:
            lo, hi = allocation_module.capacity_point_bounds(x, al['profile'], bits=200)
            v = K.alloc_capacity(al, x)
            if isinstance(v, Q):
                assert lo <= v <= hi, ('allocation boundary', al['profile'], x)
            else:
                import kernel
                assert kernel.to_arb(lo) <= v and v <= kernel.to_arb(hi)
            allocation_boundary_checks += 1

    # The active graph kernels are rectangular closed indicators in (r,s,t).
    # On full branches they vanish: even their minimum selected distance is
    # strictly larger than the support. The other two graph weights are zero.
    active = [g for g in K.graphs if g['delta']]
    assert all(not g['mixed'] for g in active)
    for g in active:
        r, s = g['win'][0][0], g['win'][1][0]
        assert r*r+s*s-2*r*s*g['t'][1] > K.support**2
    radii = sorted({Q(1), K.support} | {e for g in active for win in g['win'] for e in win[:2]})
    angles = sorted({Q(-1), Q(1)} | {e for g in active for e in g['t']})
    rpoints = radii+[(a+b)/2 for a,b in zip(radii,radii[1:])]
    tpoints = angles+[(a+b)/2 for a,b in zip(angles,angles[1:])]
    graph_boundary_checks = 0
    for r,s,t in itertools.product(rpoints,rpoints,tpoints):
        if r*r+s*s-2*r*s*t < K.support**2:
            assert K.graphs_k(r,s,t) == 0
            continue
        z = (1+t)*2*r*s/((r+s)**2-K.support**2)
        bounds = [(r,r),(s,s),(z,z)]
        reference = (lift.cyclic_decrement_lower('outside',bounds,state['graph_weights'])
                     +lift.outer_graph_lower('outside',bounds,state['outer_graph_weight'])
                     +lift.narrow_graph_lower('outside',bounds,state['narrow_graph_weight'])
                     +lift.moderate_graph_lower('outside',bounds,state['moderate_graph_weight']))
        assert reference == K.graphs_k(r,s,t), ('graph boundary',r,s,t)
        graph_boundary_checks += 1

    # Domain ownership at every jump is checked independently of samples.
    boundary_checks = 0
    for p, r in enumerate(K.breaks):
        assert K.piece_of(r) == max(0,p-1)
        if p < K.np:
            assert K.piece_of((r+K.breaks[p+1])/2) == p
        boundary_checks += 1
    assert K.piece_of(Q(999,1000)) is None and K.piece_of(K.support+Q(1,1000)) is None
    return dict(status='PASS', exact_parameter_equations=equations,
                radial_basis_polynomials=K.np*K.nc,
                radial_boundary_checks=boundary_checks+2,
                allocation_boundary_checks=allocation_boundary_checks,
                graph_boundary_checks=graph_boundary_checks,
                scope='Exact component parameters and radial polynomials; boundary regression checks. '
                      'Continuum formula interpretation is the written COMPONENT_IDENTITIES.md proof; '
                      'nonpolynomial enclosure samples alone are not an identity proof.')
