"""Independent exact implementation of the V7 kernel K of the D4 paper.

This module reads only the four files in data/v7 and the written
definitions (paper Section 3, verify/stage1_v7/COMPONENTS.md and the proof
notes in verify/stage1_v7/proofs). It imports no stage 1 or stage 5 code.
Conventions taken from reading stage 1 code are listed in CONVENTIONS.md.

K(x, y) depends on r = |x|, s = |y| and t = cos(angle(x, y)). It is the sum
of these components (COMPONENTS.md rows):

  gram       harmonic Gram kernel, origin row v and origin entry kappa (row 1)
  caps       98 pointed caps, pivot caps truncated below their cutoff (row 3)
  tents      4 negative-pole tents (row 4)
  alloc      2 antipodal allocations (row 5)
  condray    2 conditional radial rays (23-N_r)(N_R-N_r) (row 6)
  joint      2 conditional count polynomials (row 7)
  nested     nested count quadratic (row 8)
  sph        2 signed spherical energies (row 9)
  graphs     six graph energies delta_j (C_j - e_j) (rows 10-14)
  gauge      pair count gauge (row 15)

The broad19 rays of origin.json (row 2) are not added: by COMPONENTS.md row 2
they are already Gram rays of kernel.json. The origin credit (row 16) is not
part of U_raw and is not implemented.

Each component c has an origin constant A_c, a linear term h_c(r) = K_c(x,0)
and a pair kernel k_c(r, s, t) = K_c(x, y) for nonzero x, y in the support
1 <= r, s <= 7/4. Outside the support every nonorigin value is zero.

Radial convention: a radius r in (B_p, B_{p+1}] belongs to piece p, and
r = 1 belongs to piece 0. Every count indicator is 1[|x| <= R] (closed
balls), and every R is a break, so indicators are constant on each piece.

All values are exact rationals, except the parts that are not rational at
rational points: the cap probabilities and the classical allocation curve
(square roots) and f and the gauge weight (exp). For these, the point
functions return a python-flint arb ball, and an exact representation is
returned where the value is rational.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import math

import flint

PREC = 320


def exact_q(value):
    """Accept exact JSON integers/strings; never silently convert a float."""
    if isinstance(value, bool) or not isinstance(value, (int, str, Q)):
        raise ValueError('expected an exact integer or rational string')
    return Q(value)


def exact_int(value, name, minimum=0):
    q = exact_q(value)
    if q.denominator != 1 or q < minimum:
        raise ValueError(name + ': expected an integer >= ' + str(minimum))
    return q.numerator


def fq(x):
    """Fraction -> fmpq."""
    x = Q(x)
    return flint.fmpq(x.numerator, x.denominator)


def qf(x):
    """fmpq -> Fraction."""
    return Q(int(x.p), int(x.q))


def to_arb(x):
    if isinstance(x, flint.arb):
        return x
    x = Q(x)
    return flint.arb(flint.fmpq(x.numerator, x.denominator))


def chebyshev_u(L):
    """Power coefficients of U_0..U_L by the three-term recurrence."""
    rows = [[Q(1)], [Q(0), Q(2)]]
    while len(rows) <= L:
        a, b = rows[-2], rows[-1]
        nxt = [Q(0)]*(len(b)+1)
        for i, c in enumerate(b):
            nxt[i+1] += 2*c
        for i, c in enumerate(a):
            nxt[i] -= c
        rows.append(nxt)
    return rows[:L+1]


def poly_eval(coefficients, x):
    value = Q(0)
    for c in reversed(coefficients):
        value = value*x+c
    return value


class Kernel:
    def __init__(self, datadir, perturb=None):
        datadir = Path(datadir)
        kd = json.loads((datadir/'kernel.json').read_text())
        od = json.loads((datadir/'origin.json').read_text())
        sd = json.loads((datadir/'scalar.json').read_text())
        gd = json.loads((datadir/'graph_energies.json').read_text())
        assert exact_int(kd['dimension'], 'dimension') == 4
        self.deg = exact_int(kd['local_degree'], 'local_degree', 1)
        self.L = exact_int(kd['angular_degree'], 'angular_degree')
        assert (self.deg, self.L) == (6, 18), 'this checker is for the degree 6/18 V7 package'
        self.breaks = [exact_q(b) for b in kd['breaks']]
        assert len(self.breaks) == 8 and self.breaks[0] == 1 and self.breaks[-1] == Q(7, 4)
        assert all(a < b for a, b in zip(self.breaks, self.breaks[1:]))
        self.np = len(self.breaks)-1
        self.nc = self.deg+1
        self.n = self.nc*self.np
        self.support = self.breaks[-1]
        self.basis = [tuple(map(exact_q, row)) for row in kd['radial_basis_intervals']]
        assert len(self.basis) == self.np
        for p, (a, b) in enumerate(self.basis):
            assert a <= self.breaks[p] < self.breaks[p+1] <= b

        # Gram blocks: M_l = sum w v v^T / D, computed as V^T diag(w) V.
        self.D = exact_int(kd['common_denominator'], 'common_denominator', 1)
        assert self.D > 0 and len(kd['blocks']) == self.L+1
        self.M = []
        for ell, block in enumerate(kd['blocks']):
            n = exact_int(block['dimension'], 'block.dimension', 1)
            assert n == self.n+(1 if ell == 0 else 0)
            rows, wrows = [], []
            for ray in block['rays']:
                w, v = ray['weight'], ray['vector']
                assert type(w) is int and w > 0 and len(v) == n
                assert all(type(z) is int for z in v)
                if ell == 0:
                    # The origin coordinate is 0 in every ray; the origin row and
                    # entry come from origin.json (v, kappa).
                    assert v[0] == 0
                    v = v[1:]
                rows.append(v)
                wrows.append([w*z for z in v])
            V = flint.fmpz_mat(rows)
            W = flint.fmpz_mat(wrows)
            self.M.append(V.transpose()*W)  # integer matrix; divide by D
        self.rays = sum(len(b['rays']) for b in kd['blocks'])

        # Origin data of the Gram kernel.
        self.v = [exact_q(x) for x in od['v']]
        self.u = [exact_q(x) for x in od['u']]
        self.kappa = exact_q(od['kappa'])
        assert len(self.v) == len(self.u) == self.n and self.kappa >= 0

        # Scalar function f(x) = exp(-alpha |x|^2) F(alpha |x|^2).
        assert exact_int(sd['dimension'], 'scalar.dimension') == 4
        self.alpha = exact_q(sd['alpha'])
        self.F = [exact_q(x) for x in sd['F']]
        self.G = [exact_q(x) for x in sd['G']]
        assert self.alpha > 0 and len(self.F) == len(self.G) == 50
        assert exact_q(sd['cutoff']) == self.support

        R = exact_q
        self.caps = []
        for c in od['pointed_caps']:
            cap = dict(d=R(c['d']), a2=R(c['a2']), sigma=R(c['sigma']), M=exact_int(c['M'], 'cap.M'),
                       weight=R(c['weight']), w=[R(x) for x in c['w']],
                       cutoff=R(c['pivot_lower_cutoff']) if 'pivot_lower_cutoff' in c else None)
            assert len(cap['w']) == self.n and min(cap['w']) >= 0 and cap['weight'] >= 0
            assert cap['d'] > 0 and cap['a2'] > 0 and -1 <= cap['sigma'] <= 1
            self.caps.append(cap)
        self.tents = []
        for c in od['pointed_tents']:
            assert R(c.get('sigma', -1)) == -1
            self.tents.append(dict(d=R(c['d']), capacity=R(c['capacity']),
                                   weight=R(c['weight']), w=[R(x) for x in c['w']]))
            assert len(self.tents[-1]['w']) == self.n and min(self.tents[-1]['w']) >= 0
            assert self.tents[-1]['weight'] >= 0 and self.tents[-1]['d'] > 0
        self.allocs = []
        for c in od['independent_antipodal_allocations']:
            assert R(c['origin_constant']) == 0
            prof = c['capacity_profile']
            assert prof in ('classical_improvement', 'two_level')
            a = dict(weight=R(c['weight']), w=[R(x) for x in c['w']], profile=prof)
            assert len(a['w']) == self.n and min(a['w']) >= 0 and a['weight'] >= 0
            if prof == 'two_level':
                a['inner'], a['outer'] = R(c['inner_radius']), R(c['outer_radius'])
                assert 0 < a['inner'] < a['outer']
                assert c['inner_boundary'] == 'strict' and c['outer_boundary'] == 'closed'
            self.allocs.append(a)
        self.condrays = []
        for c in od['conditional_radial_rays']:
            assert R(c['pair_coefficient']) == 23 and R(c['mixed_kernel_coefficient']) == Q(-1, 2)
            self.condrays.append(dict(r=R(c['inner_radius']), R=R(c['outer_radius']),
                                      weight=R(c['weight'])))
            assert self.condrays[-1]['weight'] >= 0
        nq = od['nested_count_quadratic']
        self.counts = [dict(name='nested', radii=[R(x) for x in nq['radii']], weight=Q(1),
                            coefficients=[R(x) for x in nq['coefficients']])]
        for i, c in enumerate(od['conditional_count_polynomials']):
            self.counts.append(dict(name='joint%d' % i, radii=[R(x) for x in c['radii']],
                                    weight=R(c['weight']),
                                    coefficients=[R(x) for x in c['coefficients']]))
        for c in self.counts:
            k = len(c['radii'])
            assert len(c['coefficients']) == 1+k+k*(k+1)//2
            assert len(set(c['radii'])) == k and c['weight'] >= 0
            assert all(r in self.breaks for r in c['radii'])
        self.sph = []
        for c in od['appended_spherical_local_energy_rays']:
            assert c['scalar_basis'] == 'U_j(t)/(j+1)'
            g = [R(x) for x in c['scalar_coefficients']]
            ray = dict(R=R(c['R']), weight=R(c['weight']), g=g, B=R(c['B']),
                       kappa=R(c['kappa_s']))
            assert ray['B']-sum(g) == ray['kappa']  # B - g(1) = kappa_s
            assert ray['weight'] >= 0 and ray['R'] in self.breaks
            self.sph.append(ray)
        # Graph energies: weights from graph_energies.json, windows and
        # capacities from paper Table 1 (tab:graphs) and COMPONENTS.md 10-14.
        w3 = [R(x) for x in gd['weights']]
        assert len(w3) == 3 and min(w3) >= 0
        I = lambda a, b, lo_open=False: (R(a), R(b), lo_open)
        self.graphs = [
            dict(name='A', delta=w3[0], C=862, win=(I(1, '103/100'), I('103/100', '3/2', True)),
                 mixed=True, t=(R(-4)/5, R(-1)/2)),
            dict(name='B', delta=w3[1], C=1102, win=(I(1, '103/100'), I('3/2', '7/4')),
                 mixed=True, t=(R('-19/20'), R(-3)/4)),
            dict(name='E9', delta=w3[2], C=1358, win=(I('33/28', '7/4'), I('33/28', '7/4')),
                 mixed=False, t=(R(-1), R(-7)/8)),
            dict(name='near', delta=R(gd['outer_graph_weight']), C=120,
                 win=(I(1, '103/100'), I(1, '103/100')), mixed=False, t=(R(-1), R(-3)/4)),
            dict(name='E4', delta=R(gd['narrow_graph_weight']), C=604,
                 win=(I('33/28', '7/4'), I('33/28', '7/4')), mixed=False, t=(R(-1), R('-183/200'))),
            dict(name='moderate', delta=R(gd['moderate_graph_weight']), C=1880,
                 win=(I('7/5', '3/2'), I('7/5', '3/2')), mixed=False, t=(R(-4)/5, R('-11/20'))),
        ]
        gg = od['full_pair_count_gauge']
        assert gg['gauge'] == 'indicator_over_count_bound'
        assert all(g['delta'] >= 0 for g in self.graphs)
        self.gauge = dict(R=R(gg['radius']), bound=exact_int(gg['count_bound'], 'gauge.count_bound', 2))
        assert self.gauge['R'] in self.breaks

        if perturb is not None:
            perturb(self)

        # Chebyshev and Bernstein bases.
        self.U = chebyshev_u(self.L)
        self.Z = [[c/(ell+1) for c in row] for ell, row in enumerate(self.U)]
        self.phi_power = [self._bernstein_power(p) for p in range(self.np)]

    # ------------------------------------------------------------------
    # radial pieces and indicators
    def piece_of(self, r):
        r = Q(r)
        if r < 1 or r > self.support:
            return None
        if r == 1:
            return 0
        for p in range(self.np):
            if self.breaks[p] < r <= self.breaks[p+1]:
                return p
        raise AssertionError

    def _ind(self, p, R):
        """1[|x| <= R] on piece p (R is a break)."""
        lo, hi = self.breaks[p:p+2]
        assert hi <= R or lo >= R, 'indicator radius is not a break'
        return int(hi <= R)

    def count_flags(self, c, p):
        return [self._ind(p, R) for R in c['radii']]

    def condray_flags(self, ray, p):
        inner = self._ind(p, ray['r'])
        shell = self._ind(p, ray['R'])-inner
        return inner, shell

    def gauge_psi(self, p):
        return Q(self._ind(p, self.gauge['R']), self.gauge['bound'])

    # ------------------------------------------------------------------
    # radial basis: degree-deg Bernstein polynomials on the basis interval
    def _bernstein_power(self, p):
        a, b = self.basis[p]
        n = self.deg
        x = flint.fmpq_poly([0, 1])
        z = (x-fq(a))*fq(1/(b-a))
        rows = []
        for i in range(n+1):
            poly = math.comb(n, i)*z**i*(1-z)**(n-i)
            coeffs = [qf(c) for c in poly.coeffs()]
            rows.append(coeffs+[Q(0)]*(n+1-len(coeffs)))
        return rows  # rows[i][a] = coefficient of r^a in basis function i

    def phi_point(self, p, r):
        """Basis values on piece p at r, from the Bernstein formula directly."""
        a, b = self.basis[p]
        z = (Q(r)-a)/(b-a)
        n = self.deg
        return [math.comb(n, i)*z**i*(1-z)**(n-i) for i in range(n+1)]

    def radial_poly(self, vec, p):
        """Power coefficients (in r) of vec . phi on piece p."""
        block = vec[self.nc*p:self.nc*(p+1)]
        out = [Q(0)]*self.nc
        for i, c in enumerate(block):
            if c:
                for a, e in enumerate(self.phi_power[p][i]):
                    out[a] += c*e
        return out

    def radial_point(self, vec, p, r):
        block = vec[self.nc*p:self.nc*(p+1)]
        return sum((c*b for c, b in zip(block, self.phi_point(p, r))), Q(0))

    # ------------------------------------------------------------------
    # origin constant A = K(0,0)
    def origin_parts(self):
        parts = {'gram_kappa': self.kappa}
        for c in self.counts:
            parts['count_'+c['name']] = c['weight']*c['coefficients'][0]
        for g in self.graphs:
            parts['graph_'+g['name']] = g['delta']*g['C']
        # caps, tents, allocations, conditional rays, spherical rays and the
        # gauge have origin constant 0.
        return parts

    def A(self):
        return sum(self.origin_parts().values(), Q(0))

    def U_raw(self):
        # f(0) = F(0); paper eq. (normal): Delta_4 <= (f(0) + A)/8.
        return (self.F[0]+self.A())/8

    # ------------------------------------------------------------------
    # linear term h(r) = K(x, 0)
    def h_poly_parts(self, p):
        """Power coefficients in r of each polynomial part of h on piece p.

        Caps enter with their untruncated linear weight M w/2 here; the
        truncation below a pivot cutoff is in h_cap_truncation.
        """
        parts = {}
        parts['gram'] = self.radial_poly(self.v, p)
        acc = [Q(0)]*self.nc
        for cap in self.caps:
            for a, c in enumerate(self.radial_poly(cap['w'], p)):
                acc[a] += cap['weight']*cap['M']*c/2
        parts['caps_linear_untruncated'] = acc
        acc = [Q(0)]*self.nc
        for tent in self.tents:
            for a, c in enumerate(self.radial_poly(tent['w'], p)):
                acc[a] += tent['weight']*tent['capacity']*c/2
        parts['tents'] = acc
        acc = [Q(0)]*self.nc
        for al in self.allocs:
            for a, c in enumerate(self.radial_poly(al['w'], p)):
                acc[a] += al['weight']*c/2
        parts['alloc'] = acc
        parts['sph'] = [sum((ray['weight']*ray['B']*self._ind(p, ray['R'])/2
                             for ray in self.sph), Q(0))]
        parts['condray'] = [sum((Q(23)*ray['weight']*self.condray_flags(ray, p)[1]/2
                                 for ray in self.condrays), Q(0))]
        for c in self.counts:
            k = len(c['radii'])
            n = self.count_flags(c, p)
            parts['count_'+c['name']] = [c['weight']*sum((c['coefficients'][1+i]*n[i]
                                                         for i in range(k)), Q(0))]
        return parts

    def h_cap_truncation_point(self, r):
        """Correction to the untruncated cap linear term: -beta M w(r)/2 for r <= C."""
        p = self.piece_of(r)
        out = Q(0)
        for cap in self.caps:
            if cap['cutoff'] is not None and Q(r) <= cap['cutoff']:
                out -= cap['weight']*cap['M']*self.radial_point(cap['w'], p, r)/2
        return out

    def h_point(self, r):
        """Exact K(x,0) for |x| = r, without the gauge (gauge: gauge_w/2)."""
        p = self.piece_of(r)
        if p is None:
            return Q(0)
        total = Q(0)
        for coeffs in self.h_poly_parts(p).values():
            total += poly_eval(coeffs, Q(r))
        return total+self.h_cap_truncation_point(r)

    # ------------------------------------------------------------------
    # polynomial part of the pair kernel on a piece pair, as an mpoly in (r,s,t)
    def pair_kernel_poly(self, p, q, ctx):
        """Sum of the components of K(x,y) that are polynomial on pieces (p,q):
        gram, sph, condray, nested, joint. Returns a fmpq_mpoly in ctx(r,s,t)."""
        r, s, t = ctx.gens()
        total = ctx.from_dict({})
        # Gram: sum_l Z_l(t) * Phi_p^T M_l[p,q] Phi_q / D
        Pp = flint.fmpq_mat([[fq(c) for c in row] for row in self.phi_power[p]])
        Pq = flint.fmpq_mat([[fq(c) for c in row] for row in self.phi_power[q]])
        invD = flint.fmpq(1, self.D)
        for ell in range(self.L+1):
            Mb = flint.fmpq_mat(self.nc, self.nc,
                                [self.M[ell][self.nc*p+i, self.nc*q+j]
                                 for i in range(self.nc) for j in range(self.nc)])
            C = Pp.transpose()*Mb*Pq*invD  # C[a,b] = coefficient of r^a s^b
            radial = ctx.from_dict({(a, b, 0): C[a, b] for a in range(self.nc)
                                    for b in range(self.nc) if C[a, b] != 0})
            ang = ctx.from_dict({(0, 0, m): fq(c) for m, c in enumerate(self.Z[ell]) if c})
            total += radial*ang
        # Signed spherical: -beta g(t) chi(r) chi(s)
        for ray in self.sph:
            if self._ind(p, ray['R']) and self._ind(q, ray['R']):
                g = [Q(0)]*(len(ray['g']))
                for j, gj in enumerate(ray['g']):
                    for m, c in enumerate(self.Z[j]):
                        g[m] += gj*c
                total -= ctx.from_dict({(0, 0, m): fq(ray['weight']*c) for m, c in enumerate(g) if c})
        const = Q(0)
        # Conditional rays: -beta (a(x) b(y) + b(x) a(y))/2
        for ray in self.condrays:
            ap, bp = self.condray_flags(ray, p)
            aq, bq = self.condray_flags(ray, q)
            const -= ray['weight']*(ap*bq+bp*aq)/2
        # Count quadratics: weight * n(x)^T B n(y), B symmetric from upper triangle.
        for c in self.counts:
            const += c['weight']*self.count_bilinear(c, self.count_flags(c, p), self.count_flags(c, q))
        if const:
            total += ctx.from_dict({(0, 0, 0): fq(const)})
        return total

    def count_bilinear(self, c, n1, n2):
        k = len(c['radii'])
        idx = 1+k
        value = Q(0)
        for i in range(k):
            for j in range(i, k):
                b = c['coefficients'][idx]
                idx += 1
                value += b*(n1[i]*n2[j] if i == j else n1[i]*n2[j]+n1[j]*n2[i])
        return value

    def pair_kernel_point(self, r, s, t):
        """Exact value of the same polynomial components at one point."""
        p, q = self.piece_of(r), self.piece_of(s)
        if p is None or q is None:
            return Q(0)
        fr, fs = self.phi_point(p, r), self.phi_point(q, s)
        total = Q(0)
        for ell in range(self.L+1):
            z = poly_eval(self.Z[ell], Q(t))
            acc = Q(0)
            for i in range(self.nc):
                if fr[i]:
                    row = sum((Q(int(self.M[ell][self.nc*p+i, self.nc*q+j]))*fs[j]
                               for j in range(self.nc)), Q(0))
                    acc += fr[i]*row
            total += z*acc/self.D
        for ray in self.sph:
            if Q(r) <= ray['R'] and Q(s) <= ray['R']:
                total -= ray['weight']*sum((gj*poly_eval(self.Z[j], Q(t))
                                            for j, gj in enumerate(ray['g'])), Q(0))
        for ray in self.condrays:
            a = lambda x: int(Q(x) <= ray['r'])
            b = lambda x: int(ray['r'] < Q(x) <= ray['R'])
            total -= ray['weight']*(a(r)*b(s)+b(r)*a(s))/2
        for c in self.counts:
            n1 = [int(Q(r) <= R) for R in c['radii']]
            n2 = [int(Q(s) <= R) for R in c['radii']]
            total += c['weight']*self.count_bilinear(c, n1, n2)
        return total

    # ------------------------------------------------------------------
    # pair term q0(r) = f(r) + 2 K(x,0) + K(x,x) without the gauge
    def pair_poly(self, p):
        """D_p(r): polynomial part of 2h + k(x,x) on piece p, with caps entering
        through their untruncated linear weight only (as a comparison object).
        Returns (power coefficients, parts)."""
        ctx = flint.fmpq_mpoly_ctx.get(('r', 's', 't'), 'lex')
        k = self.pair_kernel_poly(p, p, ctx)
        diag = [Q(0)]*(2*self.nc-1)
        for (a, b, m), c in k.to_dict().items():
            diag[a+b] += qf(c)  # s = r, t = 1
        total = list(diag)
        for name, coeffs in self.h_poly_parts(p).items():
            for a, c in enumerate(coeffs):
                total[a] += 2*c
        return total

    def F_of_alpha_r2(self):
        """F(alpha r^2) as power coefficients in r."""
        out = [Q(0)]*(2*len(self.F)-1)
        for k, c in enumerate(self.F):
            out[2*k] = c*self.alpha**k
        return out

    def f_point(self, r):
        x = self.alpha*Q(r)**2
        poly = to_arb(poly_eval(self.F, x))
        return poly*(-to_arb(x)).exp()

    def cap_self_probability(self, cap, r):
        """P(x | x): v.xhat = sigma is fixed, so the event is sigma >= c(r)."""
        d, a2, sigma = cap['d'], cap['a2'], cap['sigma']
        return Q(int(sigma*2*d*Q(r) >= Q(r)**2+d*d-a2))

    def cap_pair_point(self, r):
        """Exact sum over caps of 2 h_cap + k_cap(x,x) = beta w~(r) (M - P(x|x))."""
        p = self.piece_of(r)
        out = Q(0)
        for cap in self.caps:
            if cap['cutoff'] is not None and Q(r) <= cap['cutoff']:
                continue
            out += cap['weight']*self.radial_point(cap['w'], p, r)*(cap['M']-self.cap_self_probability(cap, r))
        return out

    def J_equivalent_point(self, r):
        """The amount stage 1 calls the pointed self reduction J, from my kernel:
        J = sum beta M w(r) (untruncated linear) - sum (2 h_cap + k_cap(x,x))."""
        p = self.piece_of(r)
        lin = sum((cap['weight']*cap['M']*self.radial_point(cap['w'], p, r) for cap in self.caps), Q(0))
        return lin-self.cap_pair_point(r)

    def q0_point(self, r):
        """q0(r) = f + 2K(x,0) + K(x,x) without the gauge, as an arb ball.
        Exact rational part returned too: q0 = f + exact."""
        r = Q(r)
        exact = 2*self.h_point(r)+self.pair_kernel_point(r, r, 1)
        exact += self.cap_self_diag_point(r)
        # tents, allocations, graphs and conditional rays have zero diagonal
        return self.f_point(r)+to_arb(exact), exact

    def cap_self_diag_point(self, r):
        """k_caps(x,x) = -sum beta w~(r) P(x|x)."""
        p = self.piece_of(r)
        out = Q(0)
        for cap in self.caps:
            if cap['cutoff'] is not None and Q(r) <= cap['cutoff']:
                continue
            out -= cap['weight']*self.radial_point(cap['w'], p, r)*self.cap_self_probability(cap, r)
        return out

    def gauge_w_point(self, r):
        """w(r) = -q0(r)/(1 - psi(r)) as an arb ball."""
        q0, _ = self.q0_point(r)
        psi = self.gauge_psi(self.piece_of(r))
        return -q0/to_arb(1-psi)

    # ------------------------------------------------------------------
    # nonpolynomial pair kernels at a point (x != y, both in the support)
    def cap_probability(self, cap, r, t):
        """P(x|y) for |x| = r and cos(x,y) = t. Exact Fraction when rational,
        else an arb ball. Probe direction v with v.yhat = sigma, uniform."""
        d, a2, sigma = cap['d'], cap['a2'], cap['sigma']
        r, t = Q(r), Q(t)
        c = (r*r+d*d-a2)/(2*d*r)      # event: v.xhat >= c
        A = c-sigma*t
        S = (1-sigma*sigma)*(1-t*t)   # v.xhat = sigma t + sqrt(S) U, U ~ Unif[-1,1]
        if S == 0:
            return Q(int(A <= 0))
        if A >= 0 and A*A >= S:
            return Q(0)
        if A <= 0 and A*A >= S:
            return Q(1)
        root = Q(math.isqrt(S.numerator), math.isqrt(S.denominator))
        if root*root == S:
            return (1-A/root)/2
        return (1-to_arb(A)/to_arb(S).sqrt())/2

    def cap_weight_point(self, cap, r):
        p = self.piece_of(r)
        if cap['cutoff'] is not None and Q(r) <= cap['cutoff']:
            return Q(0)
        return self.radial_point(cap['w'], p, r)

    def caps_k(self, r, s, t):
        """Sum over caps of -k_cap(x,y) >= 0 (the subtracted correction)."""
        total = Q(0)
        for cap in self.caps:
            half = cap['weight']/2
            total = _add(total, _mul(half*self.cap_weight_point(cap, s), self.cap_probability(cap, r, t)))
            total = _add(total, _mul(half*self.cap_weight_point(cap, r), self.cap_probability(cap, s, t)))
        return total

    def tents_k(self, r, s, t):
        """Sum over tents of -k_tent(x,y) = beta (w(s) phi(x + d yhat) + w(r) phi(y + d xhat))/2."""
        r, s, t = Q(r), Q(s), Q(t)
        p, q = self.piece_of(r), self.piece_of(s)
        total = Q(0)
        for tent in self.tents:
            d = tent['d']
            px = max(Q(0), 1-2*(r*r+d*d+2*d*r*t))  # |x + d yhat|^2 = r^2 + 2 d r t + d^2
            py = max(Q(0), 1-2*(s*s+d*d+2*d*s*t))
            total += tent['weight']*(self.radial_point(tent['w'], q, s)*px
                                     + self.radial_point(tent['w'], p, r)*py)/2
        return total

    def alloc_capacity(self, al, sq):
        """Capacity as a function of the squared distance |x+y|^2 = sq."""
        sq = Q(sq)
        if al['profile'] == 'two_level':
            if sq < al['inner']**2:
                return Q(1)
            if sq <= al['outer']**2:
                return Q(1, 2)
            return Q(0)
        if sq >= Q(1, 2):
            return Q(0)
        if sq >= Q(1, 4):
            return 1-2*sq
        root = Q(math.isqrt(sq.numerator), math.isqrt(sq.denominator))
        if root*root == sq:
            return min(Q(1), 2*(1-root)**2)
        val = 2*(1-to_arb(sq).sqrt())**2
        # min(1, val): val >= 1 iff sqrt(sq) <= 1 - 1/sqrt(2) iff sq <= (1-1/sqrt2)^2
        # (1-1/sqrt2)^2 = 3/2 - sqrt2; decide exactly: sq <= 3/2 - sqrt2 iff sqrt2 <= 3/2 - sq
        if Q(3, 2)-sq >= 0 and (Q(3, 2)-sq)**2 >= 2:
            return Q(1)
        return val

    def allocs_k(self, r, s, t):
        """Sum of -k_alloc(x,y) = beta (u(r) + u(s)) Cap(|x + y|^2)/2."""
        r, s, t = Q(r), Q(s), Q(t)
        p, q = self.piece_of(r), self.piece_of(s)
        sq = r*r+s*s+2*r*s*t
        total = Q(0)
        for al in self.allocs:
            u = self.radial_point(al['w'], p, r)+self.radial_point(al['w'], q, s)
            total = _add(total, _mul(al['weight']*u/2, self.alloc_capacity(al, sq)))
        return total

    def graphs_k(self, r, s, t):
        """Sum of -k_graph(x,y) = delta * 1[(x,y) selected]."""
        r, s, t = Q(r), Q(s), Q(t)
        inwin = lambda x, w: (w[0] < x if w[2] else w[0] <= x) and x <= w[1]
        total = Q(0)
        for g in self.graphs:
            w1, w2 = g['win']
            if g['mixed']:
                radial = (inwin(r, w1) and inwin(s, w2)) or (inwin(s, w1) and inwin(r, w2))
            else:
                radial = inwin(r, w1) and inwin(s, w2)
            if radial and g['t'][0] <= t <= g['t'][1]:
                total += g['delta']
        return total

    def gauge_k(self, r, s):
        """-k_gauge(x,y) = (w(r) psi(s) + w(s) psi(r))/2, an arb ball (or 0)."""
        pr, ps = self.gauge_psi(self.piece_of(r)), self.gauge_psi(self.piece_of(s))
        total = Q(0)
        if ps:
            total = _add(total, self.gauge_w_point(r)*to_arb(ps/2))
        if pr:
            total = _add(total, self.gauge_w_point(s)*to_arb(pr/2))
        return total

    # ------------------------------------------------------------------
    # cyclic sum, polynomial part, homogenized as in stage 1's branch coordinates
    def gamma_poly(self, kind, regions):
        """factor * (polynomial components of Gamma) as a fmpq_mpoly in (x0,x1,x2).

        full:    x = side lengths (a, b, c) = (|x|, |y|, |x-y|) on pieces regions;
                 factor (2abc)^L; Gamma = k(a,b,t_c) + k(a,c,t_b) + k(b,c,t_a).
        outside: x = (a, b, z), third side c^2 = (a+b)^2 - z((a+b)^2 - R^2) with
                 R = 7/4, so t = -1 + z((a+b)^2-R^2)/(2ab); factor (2ab)^L;
                 c > R gives Gamma = k(a,b,t).
        """
        L = self.L
        ctx3 = flint.fmpq_mpoly_ctx.get(('x0', 'x1', 'x2'), 'lex')
        ctxk = flint.fmpq_mpoly_ctx.get(('r', 's', 't'), 'lex')
        X = ctx3.gens()
        total = ctx3.from_dict({})
        if kind == 'full':
            pairs = [(0, 1, 2), (0, 2, 1), (1, 2, 0)]
        else:
            pairs = [(0, 1, None)]
        for i, j, k in pairs:
            kp = self.pair_kernel_poly(regions[i], regions[j], ctxk)
            groups = {}
            for (a, b, m), c in kp.to_dict().items():
                assert m <= L
                groups.setdefault(m, {})
                e = [0, 0, 0]
                e[i] += a
                e[j] += b
                groups[m][tuple(e)] = c
            if kind == 'full':
                N = X[i]**2+X[j]**2-X[k]**2
                tail = X[k]**L
            else:
                R2 = fq(self.support**2)
                N = -2*X[0]*X[1]+X[2]*((X[0]+X[1])**2-R2)
                tail = ctx3.from_dict({(0, 0, 0): 1})
            two = 2*X[i]*X[j]
            for m, terms in groups.items():
                Sm = ctx3.from_dict(terms)
                total += Sm*(N**m)*(two**(L-m))*tail
        return total

    def factor_point(self, kind, point):
        a, b, c = map(Q, point)
        return (2*a*b*c)**self.L if kind == 'full' else (2*a*b)**self.L

    def outside_cosine(self, a, b, z):
        return -1+Q(z)*((Q(a)+Q(b))**2-self.support**2)/(2*Q(a)*Q(b))


def _add(a, b):
    if isinstance(a, flint.arb) or isinstance(b, flint.arb):
        return to_arb(a)+to_arb(b)
    return a+b


def _mul(a, b):
    if isinstance(a, flint.arb) or isinstance(b, flint.arb):
        return to_arb(a)*to_arb(b)
    return a*b


flint.ctx.prec = PREC
