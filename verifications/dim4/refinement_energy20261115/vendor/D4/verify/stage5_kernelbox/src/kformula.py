"""Lemma 1044 (kernelbox, D4-12): explicit arb formula for the raw V7 pair kernel.

K(a,b) is the polarization of the accepted raw V7 anchored energy (kernel_poly.K).
Because E is a quadratic form in the star, K is the sum of the component cross
terms. With r=|a|, s=|b|, t=cos(a,b) (Gegenbauer parameter 1, Z_d=U_d/(d+1)):

  K = Kquad + Ksph + Kcount + Kcond + Kgauge - Pos,
  Kquad  = sum_d Z_d(t) sum_rays (w/DEN) alpha_i^d(r) alpha_i^d(s)  (d=0 is E0),
  Ksph   = -sum_rows weight g_row(t) [r,s <= R_row],
  Kcount, Kcond = polarizations of the count polynomials (depend on shells only),
  Kgauge = -(1/48)(G(s) w(r) + G(r) w(s)), w = -q0/(1-psi), psi = G/24,
  Pos    = caps + tents + allocations + graphs cross terms, each >= 0.

Every function accepts arb balls, so the same code gives point values and box
enclosures. validate() compares point values with raw_v7 via kernel_poly.
Run from any folder; inputs come from paths.py.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import importlib.util, json, sys
from flint import arb, ctx

ctx.prec = 256
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths
SRC = paths.V7
kernel = json.loads((SRC / "kernel.json").read_text())
origin = json.loads((SRC / "origin.json").read_text())
scalar = json.loads((SRC / "scalar.json").read_text())
graphs = json.loads((SRC / "graph_energies.json").read_text())
ACC = json.loads(paths.ACCEPTANCE.read_text())
U_RAW = Q(ACC["exact_source_objective"])
A_CONST = 8 * U_RAW - Q(scalar["F"][0])
DEN = int(kernel["common_denominator"])
DEG = kernel["local_degree"]
BREAKS = list(map(Q, kernel["breaks"]))
INTERVALS = [tuple(map(Q, x)) for x in kernel["radial_basis_intervals"]]
RAD5 = tuple(map(Q, origin["conditional_count_polynomials"][0]["radii"]))
RAD3 = tuple(map(Q, origin["nested_count_quadratic"]["radii"]))
GAUGE_R = Q(origin["full_pair_count_gauge"]["radius"])


def A(q):
    q = Q(q)
    return arb(q.numerator) / arb(q.denominator)


# Per-degree 49x49 Gram matrices M^d = sum_rays (w/DEN) v v^T over features.
def _gram(degree):
    block = kernel["blocks"][degree]
    off = int(degree == 0)
    M = {}
    for ray in block["rays"]:
        v = [Q(x) for x in ray["vector"][off:]]
        w = Q(ray["weight"], DEN)
        nz = [(i, x) for i, x in enumerate(v) if x]
        for i, x in nz:
            for j, y in nz:
                M[i, j] = M.get((i, j), 0) + w * x * y
    return M


GRAM = [_gram(d) for d in range(len(kernel["blocks"]))]
U0 = [Q(x) for x in origin["u"]]
E0RAYS = [(Q(r["weight"], DEN), [Q(x) for x in r["vector"][1:]]) for r in kernel["blocks"][0]["rays"]]


def region_of(r_lo, r_hi):
    """Closed-inner convention: the first region [b_i, b_{i+1}] containing [r_lo, r_hi]."""
    for i in range(len(BREAKS) - 1):
        if BREAKS[i] <= r_lo <= r_hi <= BREAKS[i + 1]:
            return i
    raise ValueError("radius interval straddles a break")


def features(r, region):
    a, b = INTERVALS[region]
    tau = (r - A(a)) / A(b - a)
    om = 1 - tau
    z = {}
    for j in range(DEG + 1):
        z[region * (DEG + 1) + j] = comb(DEG, j) * tau ** j * om ** (DEG - j)
    return z


def zonal_list(t, top):
    """Z_d(t) = U_d(t)/(d+1), d = 0..top."""
    U = [arb(1), 2 * t]
    for _ in range(2, top + 1):
        U.append(2 * t * U[-1] - U[-2])
    return [U[d] / (d + 1) for d in range(top + 1)]


def quad_form(M, fa, fb):
    total = arb(0)
    for (i, j), c in M.items():
        if i in fa and j in fb:
            total += A(c) * fa[i] * fb[j]
    return total


def radial_weight(row, r, region, rq_hi):
    """rw(m); zero when m <= cutoff^2 (rq_hi is an upper bound for r, rational)."""
    if "pivot_lower_cutoff" in row:
        c = Q(row["pivot_lower_cutoff"])
        if rq_hi <= c:
            return arb(0)
        assert False, "box reaches a pivot cutoff region; not needed for radii <= 103/100"
    f = features(r, region)
    return sum((A(row["w"][k]) * v for k, v in f.items()), arb(0))


def cap_prob(row, r, t):
    """P(x|y), |x| = r, cos = t, clipped to [0,1]; sigma = +-1 gives an indicator."""
    d, a2, sg = Q(row["d"]), Q(row["a2"]), Q(row["sigma"])
    if sg * sg == 1:
        # closed event |x - sigma d yhat|^2 <= a2  <=>  r^2 + d^2 - 2 sigma d r t - a2 <= 0
        L = r * r + A(d * d - a2) - 2 * A(sg * d) * r * t
        if L <= 0:
            return arb(1)
        if L > 0:
            return arb(0)
        return arb(0.5, 0.5)
    c = (r * r + A(d * d - a2)) / (2 * A(d) * r)
    D = (A(1 - sg * sg) * (1 - t * t))
    if not (D > 0):
        # only at t = +-1 exactly; enclose by [0,1]
        return arb(0.5, 0.5)
    val = (1 - (c - A(sg) * t) / D.sqrt()) / 2
    return clip01(val)


def clip01(x):
    lo_b = x.lower(); hi_b = x.upper()
    lo2 = arb(0) if lo_b <= 0 else (arb(1) if lo_b >= 1 else arb(lo_b))
    hi2 = arb(0) if hi_b <= 0 else (arb(1) if hi_b >= 1 else arb(hi_b))
    return hull(lo2, hi2)


def hull(a, b):
    return a.union(b)


def relu(x):
    if x >= 0:
        return x
    if x <= 0:
        return arb(0)
    return hull(arb(0), arb(x.upper()))


def capacity(profile, sq):
    """Allocation capacity of a squared distance sq (arb)."""
    if profile == "two_level":
        I2 = A(Q(8453, 20000) ** 2); O2 = A(Q(11547, 20000) ** 2)
        if sq < I2:
            return arb(1)
        if sq >= I2 and sq <= O2:
            return arb(0.5)
        if sq > O2:
            return arb(0)
        parts = []
        if not (sq >= I2): parts.append(arb(1))
        if not (sq > O2) and not (sq < I2): parts.append(arb(0.5))
        if not (sq <= O2): parts.append(arb(0))
        out = parts[0]
        for p_ in parts[1:]:
            out = hull(out, p_)
        return out
    assert profile == "classical_improvement"
    half, quarter = arb(0.5), arb(0.25)
    vals = []
    # piece s >= 1/2: 0 ; 1/4 <= s < 1/2: 1-2s ; s < 1/4: min(1, 2(1-sqrt s)^2)
    lo, hi = sq.lower(), sq.upper()
    if hi >= half:
        vals.append(arb(0))
    if lo < half and hi >= quarter:
        a_ = max(lo, quarter); b_ = min(hi, half)
        seg = hull(arb(a_), arb(b_))
        vals.append(1 - 2 * seg)
    if lo < quarter:
        a_ = max(lo, arb(0)); b_ = min(hi, quarter)
        seg = hull(arb(a_), arb(b_))
        v = 2 * (1 - seg.sqrt()) ** 2
        vals.append(hull(min_arb(v.lower(), arb(1)), min_arb(v.upper(), arb(1))))
    out = vals[0]
    for v in vals[1:]:
        out = hull(out, v)
    return out


def min_arb(a, b):
    return a if a <= b else (b if b <= a else hull(a, b))


def ind_le(rq_lo, rq_hi, R):
    """Indicator r <= R for a radius interval; must be decided."""
    if rq_hi <= R:
        return 1
    if rq_lo > R:
        return 0
    raise ValueError("radius interval straddles a count radius")


def count_poly(cs, ns):
    pairs = [(i, j) for i in range(len(ns)) for j in range(i, len(ns))]
    return Q(cs[0]) + 2 * sum(Q(cs[i + 1]) * n for i, n in enumerate(ns)) + sum(
        Q(c) * ns[i] * ns[j] * (1 if i == j else 2) for c, (i, j) in zip(cs[1 + len(ns):], pairs))


def count_energy(c3, c5):
    v = count_poly(origin["nested_count_quadratic"]["coefficients"], c3)
    return v + sum(Q(row["weight"]) * count_poly(row["coefficients"], c5)
                   for row in origin["conditional_count_polynomials"])


def shell_vectors(ra):
    lo, hi = ra
    return [ind_le(lo, hi, R) for R in RAD3], [ind_le(lo, hi, R) for R in RAD5]


def K_count(ra, rb):
    a3, a5 = shell_vectors(ra); b3, b5 = shell_vectors(rb)
    s3 = [x + y for x, y in zip(a3, b3)]; s5 = [x + y for x, y in zip(a5, b5)]
    return (count_energy(s3, s5) - count_energy(a3, a5) - count_energy(b3, b5)
            + count_energy([0] * 3, [0] * 5)) / 2


def K_cond(ra, rb):
    total = Q(0)
    for row in origin["conditional_radial_rays"]:
        ri, ro = Q(row["inner_radius"]), Q(row["outer_radius"])
        ia, oa = ind_le(*ra, ri), ind_le(*ra, ro)
        ib, ob = ind_le(*rb, ri), ind_le(*rb, ro)
        f = lambda ni, no: (23 - ni) * (no - ni)
        total += Q(row["weight"]) * Q(f(ia + ib, oa + ob) - f(ia, oa) - f(ib, ob) + f(0, 0), 2)
    return total


def scalar_f(r):
    x = A(scalar["alpha"]) * r * r
    poly = arb(0)
    for c in reversed(scalar["F"]):
        poly = poly * x + A(c)
    return (-x).exp() * poly


def q0(r, ra):
    """Unlifted full pair q0(m) = f(m) + [sum comp({m}) - sum comp({})], m = r^2."""
    reg = region_of(*ra)
    f = features(r, reg)
    out = scalar_f(r)
    for w, v in E0RAYS:
        vu = sum((v[i] * U0[i] for i in range(49)), Q(0))
        vf = sum((A(v[i]) * x for i, x in f.items()), arb(0))
        out += A(w) * ((vf + A(vu)) ** 2 - A(vu) ** 2)
    for d in range(1, len(GRAM)):
        out += quad_form(GRAM[d], f, f)
    a3, a5 = shell_vectors(ra)
    out += A(count_energy(a3, a5) - count_energy([0] * 3, [0] * 5))
    for row in origin["appended_spherical_local_energy_rays"]:
        if ind_le(*ra, Q(row["R"])):
            g1 = sum((A(c) for c in row["scalar_coefficients"]), arb(0))  # Z_d(1) = 1
            out += A(row["weight"]) * (A(row["B"]) - g1)
    for row in origin["pointed_caps"]:
        rw = radial_weight(row, r, reg, ra[1])
        if rw == 0:
            continue
        p = cap_prob_self(row, r)
        out += A(row["weight"]) * rw * (A(row["M"]) - p)
    for row in origin["pointed_tents"]:
        rw = radial_weight(row, r, reg, ra[1])
        d = A(row["d"])
        out += A(row["weight"]) * rw * (A(row["capacity"]) - relu(1 - 2 * (r + d) ** 2))
    for row in origin["independent_antipodal_allocations"]:
        rw = radial_weight(row, r, reg, ra[1])
        out += A(row["weight"]) * rw * (1 - capacity(row["capacity_profile"], 4 * r * r))
    for row in origin["conditional_radial_rays"]:
        i = ind_le(*ra, Q(row["inner_radius"])); o = ind_le(*ra, Q(row["outer_radius"]))
        out += A(row["weight"]) * (23 - i) * (o - i)
    return out - 0  # graphs: h = 2m > 0, no edge


def cap_prob_self(row, r):
    """P(x|x): sigma*sigma<1 gives D2 = 0 too; closed event r(r^2 + d^2 - a2) - 2 sigma d r^2 <= 0."""
    d, a2, sg = Q(row["d"]), Q(row["a2"]), Q(row["sigma"])
    L = r * r - 2 * A(sg * d) * r + A(d * d - a2)
    if L <= 0:
        return arb(1)
    if L > 0:
        return arb(0)
    return arb(0.5, 0.5)


GRAPH_SPECS = [
    (Q(graphs["weights"][2]), Q(33, 28), Q(7, 4), Q(-1), Q(-7, 8)),
    (Q(graphs["outer_graph_weight"]), Q(1), Q(103, 100), Q(-1), Q(-3, 4)),
    (Q(graphs["narrow_graph_weight"]), Q(33, 28), Q(7, 4), Q(-1), Q(-183, 200)),
    (Q(graphs["moderate_graph_weight"]), Q(7, 5), Q(3, 2), Q(-4, 5), Q(-11, 20)),
]


def graph_pos(ra, rb, t):
    """Cross term of the graph energies in -K: weight per qualifying unordered pair."""
    total = arb(0)
    for w, lr, hr, lt, ht in GRAPH_SPECS:
        inr = lambda R: (lr <= R[0] and R[1] <= hr, R[1] < lr or R[0] > hr)
        ya, na = inr(ra); yb, nb = inr(rb)
        if na or nb:
            continue
        assert ya and yb, "box straddles a graph radius"
        # angular condition lt <= t <= ht (h < 0 holds since ht < 0)
        if t >= A(lt) and t <= A(ht):
            total += A(w)
        elif t > A(ht) or t < A(lt):
            pass
        else:
            total += hull(arb(0), A(w))
    return total


def pieces(r, s, t, ra, rb):
    """Return (smooth, pos) with -K = smooth + pos; ra, rb rational radius bounds."""
    rega, regb = region_of(*ra), region_of(*rb)
    fa, fb = features(r, rega), features(s, regb)
    Z = zonal_list(t, len(GRAM) - 1)
    kq = sum((Z[d] * quad_form(GRAM[d], fa, fb) for d in range(len(GRAM))), arb(0))
    ksph = arb(0)
    for row in origin["appended_spherical_local_energy_rays"]:
        R = Q(row["R"])
        if ind_le(*ra, R) and ind_le(*rb, R):
            g = sum((A(c) * Z[d] for d, c in enumerate(row["scalar_coefficients"])), arb(0))
            ksph -= A(row["weight"]) * g
    kc = A(K_count(ra, rb) + K_cond(ra, rb))
    Ga, Gb = ind_le(*ra, GAUGE_R), ind_le(*rb, GAUGE_R)
    wa = -q0(r, ra) / A(1 - Q(Ga, 24)) if Gb else arb(0)
    wb = -q0(s, rb) / A(1 - Q(Gb, 24)) if Ga else arb(0)
    kg = -(Gb * wa + Ga * wb) / 48
    smooth = -(kq + ksph + kc + kg)
    pos = arb(0)
    for row in origin["pointed_caps"]:
        rwb = radial_weight(row, s, regb, rb[1]); rwa = radial_weight(row, r, rega, ra[1])
        if rwa == 0 and rwb == 0:
            continue
        pos += A(row["weight"]) * (rwb * cap_prob(row, r, t) + rwa * cap_prob(row, s, t)) / 2
    for row in origin["pointed_tents"]:
        rwb = radial_weight(row, s, regb, rb[1]); rwa = radial_weight(row, r, rega, ra[1])
        d = A(row["d"])
        Ta = relu(1 - 2 * (r * r + d * d + 2 * d * r * t))
        Tb = relu(1 - 2 * (s * s + d * d + 2 * d * s * t))
        pos += A(row["weight"]) * (rwb * Ta + rwa * Tb) / 2
    for row in origin["independent_antipodal_allocations"]:
        rwb = radial_weight(row, s, regb, rb[1]); rwa = radial_weight(row, r, rega, ra[1])
        sq = r * r + s * s + 2 * r * s * t
        pos += A(row["weight"]) * (rwa + rwb) * capacity(row["capacity_profile"], sq) / 2
    pos += graph_pos(ra, rb, t)
    return smooth, pos


def minusK_point(r, s, t):
    r, s, t = Q(r), Q(s), Q(t)
    sm, po = pieces(A(r), A(s), A(t), (r, r), (s, s))
    return sm + po


if __name__ == "__main__":
    print(minusK_point(1, 1, Q(-864, 1000)))
