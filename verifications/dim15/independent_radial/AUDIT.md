# Independent radial audit: 15D curved damping

**All assigned radial gates PASS.** The exact density expression is strictly
less than **0.033671603**, with independent Arb enclosure

```
[0.0336716025337934690284942002842102922656652523924343590239489291717046249604551321481660001565748511622425177388392198059876768483078707102915740810427031 +/- 7.32e-155].
```

This folder is a portable **radial** audit, not a standalone replay of the
angular geometry. The matching separate angular audit has also reported PASS
and is identified below. The full sphere-packing conclusion requires both
components with the identical angular input. No formal proof-assistant
verification or exhaustive literature-priority claim is made.

## Frozen inputs and independent implementation

The primary input package is `research/dim15_curved_ready_20261111`.
Independent snapshots and SHA-256 values are:

- `inputs/certificate.json`:
  `25ead1a73f05106192623c3e7cdefa8885cbd607e903523d0ec8716319b4c7b9`.
- `inputs/angular_exact.json`:
  `ebeb0fd1b0532445cf274f01f6f01f1996f608fc9794eb09e3cb43f3932f6354`.
- `inputs/discovery_candidate.json`, copied from `discovery/radial_cut12.json`:
  `f7240f69bd06f78d4ba08e3885fe205fef7e2db9f9aacece6a9311915c1f27a2`.

The verifier does not trust stored `pass`, objective or interval metadata.
It uses a closed coefficient formula for the radial export, an independent
monomial Fourier-transform recurrence, complete root isolation, and a
different curvature argument for the complete pair interval. All output is
confined to this independent folder; primary and earlier packages are unchanged.

## Exact parameters, export and normalization

The checked parameters are

```
n = 15, s = 15/8, a = 12, b = 16, physical cutoff X = 12,
h(u) = exp[-12(u/s-1)-16(u/s-1)^2], u >= s,
X = pi*u.
```

Both P and H have degree exactly 63 and contain 64 rational coefficients.
Let alpha=15/2 and E_k(X)=L_k^(alpha-1)(2X)/L_k^(alpha-1)(0).
For each numerical discovery parameter v_k and C, reconstruct exactly
`Fraction(float(value)).limit_denominator(10^16)`. The explicit coefficient
formula used here is

```
[X^j] E_k(X) = (-2)^j * k! / ((k-j)! * j! * (alpha)_j).
```

This is distinct from the producer's three-term recurrence export. Exact
rational arithmetic verifies every coefficient of

```
P = sum_(k=0)^63 v_k E_k,
H = sum_(k=0)^63 (-1)^k v_k E_k,
M = (C/B) * (1000001/1000000).
```

The guard is exactly one part per million. The checked constants include

```
B = 106162452199111134301/19544918695346176,
M = 5484650694015782803573433/1698599235185778148816000000,
MB = 4936679292543458869103/281474976710656000000.
```

The stored normalization agrees exactly with H(0)>0. It is approximately
0.99999999999999646898651121476346, **not assumed to be exactly one**.
The final objective divides by its exact rational value:

```
O = (P(0)+MB)/H(0) = 25.928310066563757... .
```

Equivalently, divide P,H,M by H(0); all homogeneous feasibility conditions
are preserved and the normalized Fourier origin is exactly one.

## Independent Fourier identity and both full rays

With the usual exp(-2*pi*i*x.xi) Fourier convention,

```
Fourier[(pi*|x|^2)^k exp(-pi*|x|^2)]
 = k! L_k^(13/2)(pi*|xi|^2) exp(-pi*|xi|^2).
```

The auditor constructs ordinary L_k^(13/2)(X), with argument X rather than
2X, using its three-term recurrence. Applying this identity to every monomial
of the final P reconstructs all 64 coefficients of H exactly. This is a
separate verification from the alternating Laguerre export identity.

Certified complex-root isolation accounts for all 63 roots, with
multiplicity, of each exact rational polynomial. Every root enclosure whose
imaginary part contains zero has its real enclosure strictly below zero for
H, or strictly below 12 for P. There are no unclassified forbidden roots.
Exact endpoint signs H(0)>0 and P(12)<0 therefore prove

```
H(X)>0 for every X>=0,
P(X)<0 for every X>=12.
```

No grid or approximate root location substitutes for either unbounded ray.

## Complete curved pair interval by a distinct proof

The required pair inequality is

```
P(X) exp(-X) - M(2h-h^2) <= 0,  X>=pi*s.
```

For X in [pi*s,12], set z=X/(pi*s)-1, phi=12z+16z^2, h=exp(-phi).
Since h>0, it is equivalent to the ratio gap

```
F(X)=P(X) exp(phi-X)-2M+M exp(-phi) <= 0.
```

This audit proves that gap, rather than the producer's direct-function
midpoint Taylor expression. Put k=phi' and l=phi'', so

```
k = (12+32z)/(pi*s),  l=32/(pi*s)^2,
F'' = [P''+2(k-1)P'+((k-1)^2+l)P] exp(phi-X)
      + M(k^2-l) exp(-phi).
```

All polynomial evaluations use centered Arb polynomial compositions, with the
sum of absolute nonconstant coefficients bounding the remainder. The whole
interval is parameterized by

```
X(t)=pi*s+(12-pi*s)t,  0<=t<=1.
```

For each exact rational parameter interval, both endpoint values and the
whole-interval second derivative are enclosed at 512 bits. If m is a lower
bound for F'' on [A,B], the convexity of F(X)-(m/2)X^2 gives

```
F(X) <= max(F(A),F(B)) + max(0,-m)*(B-A)^2/8.
```

Each accepted interval has a strictly negative upper bound. Exact binary
bisection and a separate sorted endpoint check prove that all accepted
parameter intervals tile [0,1] with no gaps. The proof accepts **67 intervals**,
maximum depth **12**. The complete cover is in `pair_cover.json`.

For X>=12, the global physical ray gives f<0. Since a,b>0 and z>=0,
0<h<=1, so M(2h-h^2)>0. Thus the pair inequality holds on the entire
remaining ray automatically. The kernel is not truncated at X=12.

## Angular data and origin PSD linkage

The auditor independently reconstructs the normalized dimension-15
Gegenbauer expansion from its recurrence and verifies every angular monomial
coefficient, all harmonic coefficients nonnegative, c_0=1 and Q(1)=B.
An independent spherical-moment calculation also gives c_0=1. The stored
harmonic list extends through degree 19; the actual polynomial degree is 12.

For allowed nonzero arguments,

```
K(x,y)=(M/B)h(|x|^2)h(|y|^2)Q(xhat.yhat),
K(0,0)=MB, K(0,x)=K(x,0)=-Mh(|x|^2), K(x,x)=Mh(|x|^2)^2.
```

The degree-zero nonorigin feature is sqrt(M/B)*h and its origin value is
-sqrt(MB). The exact identity (MB)(M/B)=M^2 checks the rank-one origin Gram.
Nonnegative harmonic coefficients and the addition theorem give all remaining
PSD components.

A smooth global extension is available without changing any required value.
Choose smooth radial chi=0 near zero and chi=1 at radii at least sqrt(s), and
eta=1 near zero and eta=0 at radii at least sqrt(s). The degree-zero feature
is sqrt(M/B)*chi*h-sqrt(MB)*eta; higher harmonic features are multiplied by
chi. Quartic damping at infinity makes all features Schwartz. There is no
large-radius cutoff and no Fourier requirement on K itself.

## Exact odd-dimensional density rounding

The packing-ball radius is sqrt(s)/2. Its volume is

```
V15 = pi^7*s^7*sqrt(s)/(128*15!!)
    = 2*pi^7*s^(15/2)*8!/16!.
```

The equivalence 2*8!/16!=1/(128*15!!) is checked exactly. Using 128 terms
of the alternating arctangent series and Machin's identity
pi=16*atan(1/5)-4*atan(1/239), the code constructs rational pi_lower and
pi_upper. The alternating-series remainders prove their strict bounds.

All quantities are positive. Therefore the entirely rational inequality

```
O^2 * 4*pi_upper^14*s^15*(8!)^2/(16!)^2 < (0.033671603)^2
```

proves O*V15<0.033671603, without a floating square root or Arb-only rounding.
The rational squared upper value is stored in `radial_result.json`; a separate
512-bit Arb computation agrees. The strict target is below the supplied
newer comparator 0.03372585, by at least the relative fraction
54247/33725850. This comparison is not a new literature search.

## Angular dependency and full theorem linkage

The distinct independent angular audit is in
`research/density_curved_profile15_angular_independent_20261111`.
Its full entrypoint is `verify_all.py` (geometry core `independent_check.py`), and the angular input hash is the
same one pinned above. That auditor reports a complete curved-profile proof:
10,744 accepted rectangles, 1,041 sorting exclusions, exact Bernstein cap
partition, a unique angular sign transition, analytic tails and all-scale
ratio monotonicity. This radial replay does not rerun or replace those gates.
See `ANGULAR_DEPENDENCY.json` for the exact scope, checked frozen external
manifest/result hashes, and shared-input contract.

The dilation link is compatible with this precise h. If v>=u>=1 and lambda>=1,

```
d/dlambda [phi(lambda*v)-phi(lambda*u)]
 = (v-u)[a-2b+2b*lambda*(u+v)] >= 44(v-u) >= 0.
```

Thus the ordered radial ratios used in the single-sign-transition cyclic
argument are nonincreasing under every dilation. This is the curved profile
covered by the angular audit, not an inherited scalar-Gaussian theorem.

For completeness, the PSD/cyclic summation mechanism is as follows. For a
finite separated configuration p_1,...,p_N, sum the PSD inequalities
sum_(j,k) K(p_j-p_i,p_k-p_i)>=0 over anchors i. The all-equal terms contribute
N*K(0,0), terms with exactly two equal indices contribute the ordered sum of
T(x)=K(x,x)+K(x,0)+K(0,x), and each unordered distinct triple contributes twice
its cyclic kernel sum. The linked angular theorem makes every such cyclic
sum nonpositive. Hence sum_(i!=j)T(p_j-p_i)>=-N*K(0,0). Combining with the
all-radius pair inequality f+T<=0 yields

```
sum_(i,j) f(p_j-p_i) <= N*(f(0)+K(0,0)).
```

For large finite portions of a periodic packing, Schwartz decay makes the
boundary terms vanish per point. Poisson summation and the proved Fourier
positivity give the lower value center_density*H(0) for that pair sum per
point. Periodic approximation then yields the unrestricted sphere-packing
upper bound O*V15. This uses the cyclic theorem, not a claim that each
individual off-origin kernel value is nonpositive.

## Replay

Run `python verify_all.py` with Python 3.11+ and python-flint 0.9.0, assertions
enabled. No optimizer or access to the primary package is needed. The input
snapshots, independent implementation, full cover and exact arithmetic
suffice for every radial claim here. The separate pinned angular/cyclic
component remains an explicit dependency of the final density theorem.
