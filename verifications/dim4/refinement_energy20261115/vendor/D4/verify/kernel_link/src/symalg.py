"""Exact symbolic values for tracing the two kernel implementations.

A value is a finite sum  sum_m c_m(r,s,t) * m,  where each coefficient c_m is a
rational function over Q in (r,s,t) (FLINT fmpq_mpoly numerator and denominator,
reduced by their gcd, denominator with leading coefficient 1) and each m is a
product of integer powers of atoms. An atom is an opaque function of one value:

  sqrt(x), exp(x), inv(x)               analytic atoms
  relu(x) = max(0, x)
  clip(x) = min(1, max(0, x))
  ind(x)  = 1 if x <= 0 else 0
  band[lo,hi](x) = 1 if lo <= x <= hi else 0
  cap[classical_improvement](x), cap[two_level](x): the allocation capacities
      classical: 0 for x >= 1/2; 1 - 2x for 1/4 <= x < 1/2; min(1, 2(1 - sqrt x)^2) for 0 <= x < 1/4
      two_level: 1 for x < (8453/20000)^2; 1/2 for (8453/20000)^2 <= x <= (11547/20000)^2; 0 above

Atoms are identified by a structural key: (kind, parameters, key of the argument).
Keys of coefficients are canonical (sorted exponent/coefficient lists), so two
values with equal keys are the same function wherever every atom is defined.
Equal keys are a sufficient condition for equality, never a necessary one.

Kink constructors take a cell (a closed box in (r,s,t)) and replace a kink by its
exact branch when an Arb enclosure of its argument over the whole cell decides it.
Otherwise they return the atom. Both implementations go through the same
constructors, so equal arguments give equal results.
"""
from fractions import Fraction as Q
from flint import fmpq_mpoly_ctx, fmpq, arb, ctx as flint_ctx

flint_ctx.prec = 256
VARS = ('r', 's', 't')
CTX = fmpq_mpoly_ctx.get(VARS, 'lex')
GR, GS, GT = CTX.gens()
ONE_P = CTX.constant(1)
ZERO_P = CTX.constant(0)


def fq(x):
    x = Q(x)
    return fmpq(x.numerator, x.denominator)


def tofrac(x):
    """fmpq or exact arb point -> Fraction."""
    if isinstance(x, fmpq):
        return Q(int(x.p), int(x.q))
    if isinstance(x, arb):
        assert x.rad() == 0, 'not an exact arb point'
        m, e = x.mid().man_exp()
        m, e = int(m), int(e)
        return Q(m * 2 ** e) if e >= 0 else Q(m, 2 ** (-e))
    return Q(x)


def pconst(x):
    return CTX.constant(fq(x))


def pkey(p):
    return tuple(sorted((k, str(v)) for k, v in p.to_dict().items()))


def pfromkey(key):
    return CTX.from_dict({k: fq(Q(v)) for k, v in key})


class RF:
    """Canonical rational function n/d over Q in (r, s, t)."""
    __slots__ = ('n', 'd')

    def __init__(self, n, d=None):
        if d is None or d.is_one():
            self.n, self.d = n, ONE_P
            return
        assert not d.is_zero(), 'division by zero rational function'
        if n.is_zero():
            self.n, self.d = n, ONE_P
            return
        if d.is_constant():
            self.n, self.d = n / d.leading_coefficient(), ONE_P
            return
        g = n.gcd(d)
        if not g.is_one():
            n, d = n / g, d / g
        lc = d.leading_coefficient()
        if lc != 1:
            n, d = n / lc, d / lc
        if d.is_one():
            d = ONE_P
        self.n, self.d = n, d

    @staticmethod
    def const(x):
        return RF(pconst(x))

    def is_zero(self):
        return self.n.is_zero()

    def is_poly(self):
        return self.d.is_one()

    def __add__(self, o):
        if self.d.is_one() and o.d.is_one():
            return RF(self.n + o.n)
        if self.d == o.d:
            return RF(self.n + o.n, self.d)
        return RF(self.n * o.d + o.n * self.d, self.d * o.d)

    def __neg__(self):
        r = RF.__new__(RF)
        r.n, r.d = -self.n, self.d
        return r

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        if self.d.is_one() and o.d.is_one():
            return RF(self.n * o.n)
        return RF(self.n * o.n, self.d * o.d)

    def scale(self, q):
        q = Q(q)
        if q == 0:
            return RF(ZERO_P)
        r = RF.__new__(RF)
        r.n, r.d = self.n * fq(q), self.d
        return r

    def inv(self):
        assert not self.n.is_zero()
        return RF(self.d, self.n)

    def key(self):
        return (pkey(self.n), pkey(self.d))

    @staticmethod
    def fromkey(key):
        return RF(pfromkey(key[0]), pfromkey(key[1]))

    def __eq__(self, o):
        return isinstance(o, RF) and self.n == o.n and self.d == o.d

    __hash__ = None


# ---------------------------------------------------------------- atoms
ATOM_ARGS = {}  # atom key -> (kind, params, argument Sym); for evaluation in this process
KINDS = ('sqrt', 'exp', 'inv', 'relu', 'clip', 'ind', 'band', 'cap')


def _mono_mul(a, b):
    d = dict(a)
    for k, e in b:
        d[k] = d.get(k, 0) + e
        if d[k] == 0:
            del d[k]
    return tuple(sorted(d.items()))


class Sym:
    """Exact value: dict monomial -> RF. The empty monomial is the atom-free part."""
    __slots__ = ('terms',)
    BOX = None          # current cell ((rlo,rhi),(slo,shi),(tlo,thi)), Fractions
    GENERIC_T = None    # (tlo, thi) used by lower(): see lower()
    LOWER_LOG = []

    def __init__(self, terms=None):
        self.terms = {m: c for m, c in (terms or {}).items() if not c.is_zero()}

    # -- construction
    @staticmethod
    def const(x):
        if isinstance(x, Sym):
            return x
        return Sym({(): RF.const(Q(x))})

    @staticmethod
    def poly(p):
        return Sym({(): RF(p)})

    @staticmethod
    def rf(rf):
        return Sym({(): rf})

    @staticmethod
    def atom(kind, params, arg):
        key = (kind, params, arg.key())
        ATOM_ARGS.setdefault(key, (kind, params, arg))
        return Sym({((key, 1),): RF.const(1)})

    @staticmethod
    def _wrap(o):
        if isinstance(o, Sym):
            return o
        if isinstance(o, (int, Q)):
            return Sym.const(o)
        if isinstance(o, RF):
            return Sym.rf(o)
        raise TypeError('Sym cannot combine with %r' % type(o))

    # -- structure
    def key(self):
        return tuple(sorted((m, c.key()) for m, c in self.terms.items()))

    @staticmethod
    def fromkey(key):
        return Sym({m: RF.fromkey(c) for m, c in key})

    def free(self):
        """Atom-free part as RF."""
        return self.terms.get((), RF(ZERO_P))

    def is_free(self):
        return all(m == () for m in self.terms)

    def is_const(self):
        if not self.is_free():
            return False
        c = self.free()
        return c.is_poly() and c.n.is_constant()

    def const_value(self):
        assert self.is_const()
        c = self.free().n
        return tofrac(c.leading_coefficient()) if not c.is_zero() else Q(0)

    def is_zero(self):
        return not self.terms

    # -- arithmetic
    def __add__(self, o):
        o = Sym._wrap(o)
        t = dict(self.terms)
        for m, c in o.terms.items():
            t[m] = t[m] + c if m in t else c
        return Sym(t)

    __radd__ = __add__

    def __neg__(self):
        return Sym({m: -c for m, c in self.terms.items()})

    def __sub__(self, o):
        return self + (-Sym._wrap(o))

    def __rsub__(self, o):
        return Sym._wrap(o) + (-self)

    def scale(self, q):
        q = Q(q)
        return Sym({m: c.scale(q) for m, c in self.terms.items()})

    def __mul__(self, o):
        if isinstance(o, (int, Q)):
            return self.scale(o)
        o = Sym._wrap(o)
        t = {}
        for m1, c1 in self.terms.items():
            for m2, c2 in o.terms.items():
                m = _mono_mul(m1, m2)
                c = c1 * c2
                t[m] = t[m] + c if m in t else c
        return Sym(t)

    __rmul__ = __mul__

    def inv(self):
        assert self.terms, 'inverse of zero'
        if self.is_free():
            return Sym.rf(self.free().inv())
        if len(self.terms) == 1:
            (m, c), = self.terms.items()
            return Sym({tuple((k, -e) for k, e in m): c.inv()})
        return Sym.atom('inv', (), self)

    def __truediv__(self, o):
        if isinstance(o, (int, Q)):
            return self.scale(1 / Q(o))
        return self * Sym._wrap(o).inv()

    def __rtruediv__(self, o):
        return self.inv() * o

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = Sym.const(1)
        for _ in range(n):
            out = out * self
        return out

    def exp(self):
        return Sym.atom('exp', (), self)

    def sqrt(self):
        return Sym.atom('sqrt', (), self)

    # -- exact comparisons with numbers (used as 'x == 0' by kformula)
    def __eq__(self, o):
        if isinstance(o, (int, Q)):
            return (self - Sym.const(o)).is_zero()
        if isinstance(o, Sym):
            return (self - o).is_zero()
        return NotImplemented

    __hash__ = None

    # -- enclosures over the current cell
    def rng(self, box=None):
        v = eval_arb(self, box or Sym.BOX)
        if not v.is_finite():
            return None
        return tofrac(v.lower()), tofrac(v.upper())

    def lower(self):
        """Used by kcert.cap_prob_lb as 'D.lower() > 0'. We prove the argument is
        c*(1-t^2), c>0, then use an exact rational subinterval only to select
        the positive branch. Its validity throughout -1<t<1 is algebraic.
        Zero denotes the separate endpoint trace. Every call is logged."""
        Sym.LOWER_LOG.append(self.key())
        if self.is_zero():
            return Q(0)
        assert self.is_free() and self.free().is_poly()
        p = self.free().n.to_dict()
        c = tofrac(p.get((0, 0, 0), fmpq(0)))
        assert c > 0 and p == {(0, 0, 0): fq(c), (0, 0, 2): fq(-c)}, \
            'only c*(1-t^2), c>0 may select the smooth cap branch'
        lo, hi = Sym.GENERIC_T
        assert -1 < lo <= hi < 1
        # Exact rational decision; an Arb ball near -1 can extend across -1
        # because its radius has limited precision even when ctx.prec is high.
        return c * min(1-lo*lo, 1-hi*hi)

    def __repr__(self):
        return 'Sym(%d terms)' % len(self.terms)


# ---------------------------------------------------------------- Arb evaluation
def _ball(lo, hi):
    lo, hi = Q(lo), Q(hi)
    a = arb(lo.numerator) / lo.denominator
    b = arb(hi.numerator) / hi.denominator
    return a.union(b)


def eval_poly_arb(p, vals):
    total = arb(0)
    for (i, j, k), c in p.to_dict().items():
        term = arb(int(c.p)) / int(c.q)
        if i:
            term *= vals[0] ** i
        if j:
            term *= vals[1] ** j
        if k:
            term *= vals[2] ** k
        total += term
    return total


def eval_rf_arb(c, vals):
    n = eval_poly_arb(c.n, vals)
    if c.d.is_one():
        return n
    return n / eval_poly_arb(c.d, vals)


def eval_atom_arb(key, vals):
    kind, params, _ = key
    x = eval_arb(ATOM_ARGS[key][2], None, vals)
    if kind == 'sqrt':
        return x.sqrt() if x > 0 else arb('nan')
    if kind == 'exp':
        return x.exp()
    if kind == 'inv':
        return 1 / x
    lo, hi = x.lower(), x.upper()
    if kind == 'relu':
        if hi <= 0:
            return arb(0)
        if lo >= 0:
            return x
        return arb(0).union(hi)
    if kind == 'clip':
        if hi <= 0:
            return arb(0)
        if lo >= 1:
            return arb(1)
        if lo >= 0 and hi <= 1:
            return x
        return arb(0).union(arb(1))
    if kind == 'ind':
        if hi <= 0:
            return arb(1)
        if lo > 0:
            return arb(0)
        return arb(0).union(arb(1))
    if kind == 'band':
        a, b = map(Q, params)
        if lo >= _ball(a, a) and hi <= _ball(b, b):
            return arb(1)
        if hi < _ball(a, a) or lo > _ball(b, b):
            return arb(0)
        return arb(0).union(arb(1))
    if kind == 'cap':
        return arb(0).union(arb(1))  # enclosure only; never used to decide equality
    raise ValueError(kind)


def eval_arb(sym, box, vals=None):
    if vals is None:
        vals = tuple(_ball(a, b) for a, b in box)
    total = arb(0)
    for m, c in sym.terms.items():
        term = eval_rf_arb(c, vals)
        for k, e in m:
            a = eval_atom_arb(k, vals)
            term *= a ** e if e > 0 else (1 / a) ** (-e)
        total += term
    return total


# ---------------------------------------------------------------- canonical kink arguments
def _content_normalize(p):
    """Divide a polynomial by the absolute value of its leading coefficient."""
    lc = tofrac(p.leading_coefficient())
    return p / fq(abs(lc)), abs(lc)


def positive_on_cell(p):
    lo_hi = Sym.poly(p).rng()
    return lo_hi is not None and lo_hi[0] > 0


def _ind_arg(x):
    """[x <= 0] is invariant under multiplying x by a function > 0. If x is free of
    atoms, with denominator > 0 on the cell, use its numerator divided by |leading
    coefficient|. Returns a canonical Sym argument."""
    if x.is_free():
        c = x.free()
        if not c.is_zero() and positive_on_cell(c.d):
            n, _ = _content_normalize(c.n)
            return Sym.poly(n)
    return x


# ---------------------------------------------------------------- kinks (exact, cell-decided)
def relu(x):
    x = Sym._wrap(x)
    if x.is_const():
        return Sym.const(max(Q(0), x.const_value()))
    rg = x.rng()
    if rg is not None and rg[0] >= 0:
        return x
    if rg is not None and rg[1] <= 0:
        return Sym.const(0)
    if x.is_free() and x.free().is_poly() and not x.free().n.is_zero():
        n, k = _content_normalize(x.free().n)
        return Sym.atom('relu', (), Sym.poly(n)).scale(k)
    return Sym.atom('relu', (), x)


def clip01(x):
    x = Sym._wrap(x)
    if x.is_const():
        return Sym.const(min(Q(1), max(Q(0), x.const_value())))
    rg = x.rng()
    if rg is not None:
        lo, hi = rg
        if lo >= 0 and hi <= 1:
            return x
        if hi <= 0:
            return Sym.const(0)
        if lo >= 1:
            return Sym.const(1)
    return Sym.atom('clip', (), x)


def ind_le0(x):
    """Exact indicator [x <= 0] (closed)."""
    x = Sym._wrap(x)
    if x.is_const():
        return Sym.const(int(x.const_value() <= 0))
    rg = x.rng()
    if rg is not None and rg[1] <= 0:
        return Sym.const(1)
    if rg is not None and rg[0] > 0:
        return Sym.const(0)
    return Sym.atom('ind', (), _ind_arg(x))


def band(x, lo, hi):
    """Exact indicator [lo <= x <= hi] (closed)."""
    x = Sym._wrap(x)
    lo, hi = Q(lo), Q(hi)
    if x.is_const():
        return Sym.const(int(lo <= x.const_value() <= hi))
    rg = x.rng()
    if rg is not None:
        if rg[0] >= lo and rg[1] <= hi:
            return Sym.const(1)
        if rg[1] < lo or rg[0] > hi:
            return Sym.const(0)
    return Sym.atom('band', (str(lo), str(hi)), x)


I2 = Q(8453, 20000) ** 2
O2 = Q(11547, 20000) ** 2


def capacity(profile, x):
    x = Sym._wrap(x)
    if x.is_const():
        v = x.const_value()
        assert v >= 0
        if profile == 'two_level':
            return Sym.const(1 if v < I2 else Q(1,2) if v <= O2 else 0)
        assert profile == 'classical_improvement'
        if v >= Q(1,2):
            return Sym.const(0)
        if v >= Q(1,4):
            return Sym.const(1-2*v)
    rg = x.rng()
    if profile == 'two_level':
        if rg is not None:
            if rg[1] < I2:
                return Sym.const(1)
            if rg[0] >= I2 and rg[1] <= O2:
                return Sym.const(Q(1, 2))
            if rg[0] > O2:
                return Sym.const(0)
        return Sym.atom('cap', ('two_level',), x)
    assert profile == 'classical_improvement', profile
    if rg is not None:
        if rg[0] >= Q(1, 2):
            return Sym.const(0)
        if rg[0] >= Q(1, 4) and rg[1] <= Q(1, 2):
            return 1 - 2 * x
    return Sym.atom('cap', ('classical_improvement',), x)


def var(name):
    return Sym.poly({'r': GR, 's': GS, 't': GT}[name])
