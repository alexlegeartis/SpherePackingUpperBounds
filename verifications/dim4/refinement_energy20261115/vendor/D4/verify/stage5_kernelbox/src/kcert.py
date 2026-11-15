"""Lemma 1044 (kernelbox): box lower bounds for -K(r,s,t) by a centered form.

Each quantity carries its value and gradient at the box centre and enclosures
of its value, gradient and Hessian over the box (class T). The range lies in
the value enclosure and in the second-order Taylor bound with the Hessian
enclosure as Lagrange remainder; every operation used below is C^2 on the box.

Kinks (relu, clip, step capacities, cap and graph indicators) are replaced on a
box where they are undecided by a CONSTANT LOWER BOUND of that term over the
box (zero derivatives). This is sound because -K is nondecreasing in each such term:
each enters Pos with a nonnegative coefficient, or enters q0 with a sign that
makes -K nondecreasing (derivation in kernelbox.md, section 2). So the triple
for -K is then the centered form of a C^2 function that is <= -K on the box.

The formula is that of kformula.py, which validate_kformula.py checks against
raw_v7 at 72 rational points. Run from any folder.
"""
from fractions import Fraction as Q
from math import comb
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kformula as kf
from flint import arb

A = kf.A
ZERO = arb(0)


def lo_of(x):
    return x.lower()


def hi_of(x):
    return x.upper()


Z3 = (ZERO, ZERO, ZERO)
ZH = ((ZERO,) * 3,) * 3


class T:
    """Second-order centered form. c, gc: value and gradient at the box centre;
    v, gv, H: value, gradient and Hessian enclosed over the box. On the box,
    f in c + gc.d + (1/2) d^T H d with |d_i| <= rad_i (Taylor with Lagrange
    remainder), and f in v."""
    __slots__ = ("c", "gc", "v", "gv", "H")
    R = None  # box half-widths (arb, >= 0), set per box

    def __init__(self, c, gc, v, gv, H):
        self.c, self.gc, self.v, self.gv, self.H = c, gc, v, gv, H

    @staticmethod
    def const(x):
        return T(x, Z3, x, Z3, ZH)

    def rng(self):
        R = T.R
        lo = self.c; hi = self.c
        for i in range(3):
            if R[i] == 0:
                continue
            a = abs(self.gc[i]).upper() * R[i]
            lo = lo - a; hi = hi + a
            hii = self.H[i][i]
            lo = lo + (hii.lower() if hii.lower() < 0 else ZERO) * R[i] * R[i] / 2
            hi = hi + (hii.upper() if hii.upper() > 0 else ZERO) * R[i] * R[i] / 2
            for j in range(i + 1, 3):
                if R[j] == 0:
                    continue
                b = abs(self.H[i][j]).upper() * R[i] * R[j]
                lo = lo - b; hi = hi + b
        lo, hi = lo.lower(), hi.upper()
        l1, u1 = self.v.lower(), self.v.upper()
        lo = l1 if l1 > lo else lo   # uncertain comparison keeps the centred bound (valid)
        hi = u1 if u1 < hi else hi
        return lo, hi

    def lower(self):
        return self.rng()[0]

    def __add__(self, o):
        if isinstance(o, T):
            return T(self.c + o.c, tuple(a + b for a, b in zip(self.gc, o.gc)), self.v + o.v,
                     tuple(a + b for a, b in zip(self.gv, o.gv)),
                     tuple(tuple(a + b for a, b in zip(ra, rb)) for ra, rb in zip(self.H, o.H)))
        o = _arb(o)
        return T(self.c + o, self.gc, self.v + o, self.gv, self.H)
    __radd__ = __add__

    def scale(self, k):
        return T(self.c * k, tuple(a * k for a in self.gc), self.v * k, tuple(a * k for a in self.gv),
                 tuple(tuple(a * k for a in row) for row in self.H))

    def __neg__(self):
        return self.scale(arb(-1))

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, T):
            return self.scale(_arb(o))
        a, b = self, o
        H = tuple(tuple(a.H[i][j] * b.v + a.v * b.H[i][j] + a.gv[i] * b.gv[j] + a.gv[j] * b.gv[i]
                        for j in range(3)) for i in range(3))
        return T(a.c * b.c, tuple(x * b.c + a.c * y for x, y in zip(a.gc, b.gc)), a.v * b.v,
                 tuple(x * b.v + a.v * y for x, y in zip(a.gv, b.gv)), H)
    __rmul__ = __mul__

    def _unary(self, f, f1, f2, fc, f1c):
        """f, f1, f2 on the box range; fc, f1c at the centre."""
        H = tuple(tuple(f1 * self.H[i][j] + f2 * self.gv[i] * self.gv[j] for j in range(3)) for i in range(3))
        return T(fc, tuple(f1c * x for x in self.gc), f, tuple(f1 * x for x in self.gv), H)

    def inv(self):
        assert self.v > 0 or self.v < 0
        u, uc = self.v, self.c
        return self._unary(1 / u, -1 / (u * u), 2 / (u * u * u), 1 / uc, -1 / (uc * uc))

    def __truediv__(self, o):
        if isinstance(o, T):
            return self * o.inv()
        return self.scale(1 / _arb(o))

    def __rtruediv__(self, o):
        return self.inv().scale(_arb(o))

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = T.const(arb(1))
        for _ in range(n):
            out = out * self
        return out

    def exp(self):
        e, ec = self.v.exp(), self.c.exp()
        return self._unary(e, e, e, ec, ec)

    def sqrt(self):
        assert self.v > 0
        sv, sc = self.v.sqrt(), self.c.sqrt()
        return self._unary(sv, 1 / (2 * sv), -1 / (4 * sv * self.v), sc, 1 / (2 * sc))


def _arb(o):
    if isinstance(o, arb):
        return o
    if isinstance(o, (int, Q)):
        return A(o)
    raise TypeError(type(o))


# ---- kinks: decided -> exact branch; undecided -> constant lower bound ----
def relu_lb(x):
    lo, hi = x.rng()
    if lo >= 0:
        return x
    if hi <= 0:
        return T.const(ZERO)
    return T.const(ZERO)


def clip01_lb(x):
    lo, hi = x.rng()
    if lo >= 0 and hi <= 1:
        return x
    if hi <= 0:
        return T.const(ZERO)
    if lo >= 1:
        return T.const(arb(1))
    if lo <= 0:
        return T.const(ZERO)
    return T.const(lo)  # 0 < lo < 1: clip(x) >= lo on the box


def ind_lb(L):
    """Lower bound of the indicator [L <= 0]."""
    lo, hi = L.rng()
    if hi <= 0:
        return T.const(arb(1))
    return T.const(ZERO)  # decided 0, or undecided (lower bound 0)


def cap_prob_lb(row, r, t):
    d, a2, sg = Q(row["d"]), Q(row["a2"]), Q(row["sigma"])
    if sg * sg == 1:
        return ind_lb(r * r + A(d * d - a2) - 2 * A(sg * d) * r * t)
    D = (1 - t * t) * A(1 - sg * sg)
    if not (D.lower() > 0):
        return T.const(ZERO)
    c = (r * r + A(d * d - a2)) / (2 * A(d) * r)
    return clip01_lb((1 - (c - A(sg) * t) / D.sqrt()) * A(Q(1, 2)))


I2 = Q(8453, 20000) ** 2
O2 = Q(11547, 20000) ** 2
S0 = None


def capacity_lb(profile, sq):
    """Both capacities are nonincreasing in sq; undecided pieces use the value at sup sq."""
    lo, hi = sq.rng()
    if profile == "two_level":
        if hi < A(I2):
            return T.const(arb(1))
        if hi <= A(O2):
            return T.const(arb(0.5))
        return T.const(ZERO)
    assert profile == "classical_improvement"
    if lo >= A(Q(1, 2)):
        return T.const(ZERO)
    if lo >= A(Q(1, 4)) and hi <= A(Q(1, 2)):
        return 1 - 2 * sq
    # s in [s0, 1/4]: 2(1-sqrt s)^2 <= 1 iff sqrt s >= 1 - 1/sqrt2
    s0 = (1 - 1 / arb(2).sqrt()) ** 2
    if lo > s0 and hi <= A(Q(1, 4)):
        return 2 * (1 - sq.sqrt()) ** 2
    if hi < s0:
        return T.const(arb(1))
    # undecided: value at sup sq is a lower bound
    h = hi
    if h >= A(Q(1, 2)):
        val = ZERO
    elif h >= A(Q(1, 4)):
        val = (1 - 2 * h).lower()
    else:
        hh = h if h > 0 else ZERO
        v = 2 * (1 - hh.sqrt()) ** 2
        val = v.lower() if v.lower() < 1 else arb(1)
        if not (val < 1 or val >= 1):
            val = val.lower()
    if not (val >= 0):
        val = ZERO
    return T.const(val)


def graph_lb(ra, rb, t):
    total = T.const(ZERO)
    for w, lr, hr, lt, ht in kf.GRAPH_SPECS:
        inr = lambda R: (lr <= R[0] and R[1] <= hr, R[1] < lr or R[0] > hr)
        ya, na = inr(ra); yb, nb = inr(rb)
        if na or nb:
            continue
        assert ya and yb, "box straddles a graph radius"
        lo, hi = t.rng()
        if lo >= A(lt) and hi <= A(ht):
            total = total + A(w)
    return total


# ---- smooth pieces ----
def features(r, region):
    a, b = kf.INTERVALS[region]
    tau = (r - A(a)) * A(1 / (b - a))
    om = 1 - tau
    tp = [T.const(arb(1))]; op = [T.const(arb(1))]
    for _ in range(kf.DEG):
        tp.append(tp[-1] * tau); op.append(op[-1] * om)
    return {region * (kf.DEG + 1) + j: comb(kf.DEG, j) * tp[j] * op[kf.DEG - j] for j in range(kf.DEG + 1)}


_GA = {}


def gram_arb(d, region_a, region_b):
    key = (d, region_a, region_b)
    if key not in _GA:
        ia = range(region_a * 7, region_a * 7 + 7); ib = range(region_b * 7, region_b * 7 + 7)
        M = kf.GRAM[d]
        _GA[key] = [[A(M.get((i, j), 0)) for j in ib] for i in ia]
    return _GA[key]


def quad(d, fa, rega, fb, regb):
    G = gram_arb(d, rega, regb)
    out = T.const(ZERO)
    for i in range(7):
        row = sum((G[i][j] * fb[regb * 7 + j] for j in range(7) if G[i][j] != 0), T.const(ZERO))
        out = out + fa[rega * 7 + i] * row
    return out


def zonals(t, top):
    U = [T.const(arb(1)), 2 * t]
    for _ in range(2, top + 1):
        U.append(2 * t * U[-1] - U[-2])
    return [U[d] * A(Q(1, d + 1)) for d in range(top + 1)]


_RW = {}


def rw(row, rowid, f, region, ra):
    if "pivot_lower_cutoff" in row:
        c = Q(row["pivot_lower_cutoff"])
        assert ra[1] <= c, "pivot cutoff region not handled"
        return None
    key = (rowid, region)
    if key not in _RW:
        _RW[key] = [A(row["w"][region * 7 + j]) for j in range(7)]
    ws = _RW[key]
    return sum((ws[j] * f[region * 7 + j] for j in range(7)), T.const(ZERO))


def scalar_f(r):
    x = A(kf.scalar["alpha"]) * r * r
    poly = T.const(ZERO)
    for c in reversed(kf.scalar["F"]):
        poly = poly * x + A(c)
    return (-x).exp() * poly


def q0(r, f, reg, ra):
    out = scalar_f(r)
    for w, v in kf.E0RAYS:
        vu = sum((v[i] * kf.U0[i] for i in range(49)), Q(0))
        vf = sum((A(v[reg * 7 + j]) * f[reg * 7 + j] for j in range(7)), T.const(ZERO))
        out = out + A(w) * ((vf + A(vu)) ** 2 - A(vu) ** 2)
    for d in range(1, len(kf.GRAM)):
        out = out + quad(d, f, reg, f, reg)
    a3, a5 = kf.shell_vectors(ra)
    out = out + A(kf.count_energy(a3, a5) - kf.count_energy([0] * 3, [0] * 5))
    for row in kf.origin["appended_spherical_local_energy_rays"]:
        if kf.ind_le(*ra, Q(row["R"])):
            out = out + A(row["weight"]) * (A(row["B"]) - sum((A(c) for c in row["scalar_coefficients"]), ZERO))
    for k, row in enumerate(kf.origin["pointed_caps"]):
        w_ = rw(row, ("cap", k), f, reg, ra)
        if w_ is None:
            continue
        d, a2, sg = Q(row["d"]), Q(row["a2"]), Q(row["sigma"])
        P = ind_lb(r * r - 2 * A(sg * d) * r + A(d * d - a2))
        # q0 enters -K with a negative sign, so an exact P is required here: assert decided
        lo, hi = (r * r - 2 * A(sg * d) * r + A(d * d - a2)).rng()
        assert lo > 0 or hi <= 0, "self cap undecided"
        out = out + A(row["weight"]) * w_ * (A(row["M"]) - P)
    for k, row in enumerate(kf.origin["pointed_tents"]):
        w_ = rw(row, ("tent", k), f, reg, ra)
        lo, hi = (1 - 2 * (r + A(row["d"])) ** 2).rng()
        assert hi <= 0
        out = out + A(row["weight"]) * w_ * A(row["capacity"])
    for k, row in enumerate(kf.origin["independent_antipodal_allocations"]):
        w_ = rw(row, ("alloc", k), f, reg, ra)
        lo, hi = (4 * r * r).rng()
        assert lo >= A(Q(1, 2))  # capacity(4 m) = 0
        out = out + A(row["weight"]) * w_
    for row in kf.origin["conditional_radial_rays"]:
        i = kf.ind_le(*ra, Q(row["inner_radius"])); o = kf.ind_le(*ra, Q(row["outer_radius"]))
        out = out + A(Q(row["weight"]) * (23 - i) * (o - i))
    return out


def minusK_lb(r, s, t, ra, rb):
    """Centered-form triple of a C^1 lower bound of -K on the box."""
    rega, regb = kf.region_of(*ra), kf.region_of(*rb)
    fa, fb = features(r, rega), features(s, regb)
    Z = zonals(t, len(kf.GRAM) - 1)
    kq = sum((Z[d] * quad(d, fa, rega, fb, regb) for d in range(len(kf.GRAM))), T.const(ZERO))
    ksph = T.const(ZERO)
    for row in kf.origin["appended_spherical_local_energy_rays"]:
        R = Q(row["R"])
        if kf.ind_le(*ra, R) and kf.ind_le(*rb, R):
            g = sum((A(c) * Z[d] for d, c in enumerate(row["scalar_coefficients"])), T.const(ZERO))
            ksph = ksph - A(row["weight"]) * g
    kc = A(kf.K_count(ra, rb) + kf.K_cond(ra, rb))
    Ga, Gb = kf.ind_le(*ra, kf.GAUGE_R), kf.ind_le(*rb, kf.GAUGE_R)
    kg = T.const(ZERO)
    if Gb:
        kg = kg + q0(r, fa, rega, ra) * A(Q(1, 48) / (1 - Q(Ga, 24)))
    if Ga:
        kg = kg + q0(s, fb, regb, rb) * A(Q(1, 48) / (1 - Q(Gb, 24)))
    out = -(kq + ksph + kg) - kc
    for k, row in enumerate(kf.origin["pointed_caps"]):
        wa = rw(row, ("cap", k), fa, rega, ra); wb = rw(row, ("cap", k), fb, regb, rb)
        if wa is None:
            continue
        out = out + A(Q(row["weight"]) / 2) * (wb * cap_prob_lb(row, r, t) + wa * cap_prob_lb(row, s, t))
    for k, row in enumerate(kf.origin["pointed_tents"]):
        wa = rw(row, ("tent", k), fa, rega, ra); wb = rw(row, ("tent", k), fb, regb, rb)
        d = A(row["d"])
        Ta = relu_lb(1 - 2 * (r * r + d * d + 2 * d * r * t))
        Tb = relu_lb(1 - 2 * (s * s + d * d + 2 * d * s * t))
        out = out + A(Q(row["weight"]) / 2) * (wb * Ta + wa * Tb)
    for k, row in enumerate(kf.origin["independent_antipodal_allocations"]):
        wa = rw(row, ("alloc", k), fa, rega, ra); wb = rw(row, ("alloc", k), fb, regb, rb)
        sq = r * r + s * s + 2 * r * s * t
        out = out + A(Q(row["weight"]) / 2) * (wa + wb) * capacity_lb(row["capacity_profile"], sq)
    out = out + graph_lb(ra, rb, t)
    return out


def box_vars(box):
    """box = ((r0,r1),(s0,s1),(t0,t1)) rationals -> T variables and set T.R."""
    vs, R = [], []
    for i, (a, b) in enumerate(box):
        c = (a + b) / 2; h = (b - a) / 2
        e = [ZERO, ZERO, ZERO]; e[i] = arb(1); e = tuple(e)
        vs.append(T(A(c), e, A(a).union(A(b)), e, ZH))
        R.append(A(h))
    T.R = R
    return vs


EPS = Q(1, 10 ** 40)  # below every gap between distinct radius thresholds


def dec(rint, open_lo):
    """Radius bounds used for threshold decisions. open_lo: the box excludes its
    lower endpoint (a break R_k, which belongs to the lower region by the
    closed-inner convention and is covered by the neighbouring box)."""
    lo, hi = rint
    return (lo + EPS, hi) if open_lo else (lo, hi)


def F_box(box, lam, a, open_r=False, open_s=False):
    """Centered form of a C^2 lower bound of -K - lam q(t) on the box."""
    r, s, t = box_vars(box)
    mk = minusK_lb(r, s, t, dec(box[0], open_r), dec(box[1], open_s))
    q = t * t * (t * t - A(a * a))
    return mk - q.scale(A(lam))


def F_lower(box, lam, a, open_r=False, open_s=False):
    return F_box(box, lam, a, open_r, open_s).lower()


def split_scores(f):
    R = T.R
    out = []
    for i in range(3):
        sc = abs(f.gc[i]).upper() * R[i] + abs(f.H[i][i]).upper() * R[i] * R[i] / 2
        for j in range(3):
            if j != i:
                sc += abs(f.H[i][j]).upper() * R[i] * R[j] / 2
        out.append(float(sc.mid()))
    return out
