# Independent audit: 46D noncompact kernel with an extended pair interval

**PASS. The exact frozen construction proves Delta_46 < 0.0000008777229084.**

The independently computed density enclosure is centered at
8.7772290836824778665501244491821618649481951033444446311185195806370359554659247929958469600595550870502668159237640525740203904264426042549214277969722e-7,
with radius below 5.1e-159 at 512-bit precision.

This is an independent computer-assisted audit of the supplied construction,
not a proof-assistant formalization, external peer review, or literature
priority determination. The pinned public-table comparison is separate from
the certificate.

## Inputs and independence

The audited primary package is
`research/dim46_extended_pair_ready_20261111/`.

The frozen input snapshots in this independent folder are:

- `inputs/certificate.json`: the NEW exact degree-111 P,H,M certificate.
- `inputs/angular46.json`: the target-dimensional degree-19 angular polynomial.
- `inputs/discovery_candidate.json`: the final degree-47 refit parameters.
- `inputs/source48_radial.json`: the degree-111 reference physical polynomial
  and its original Fourier normalization.

SHA-256 values are pinned in `INPUT_SHA256.json` and checked by `verify_all.py`.
The certificate SHA-256 is
`1eb9563cadf892101a24b87e5e9bfb2686b1c216297ebe29521ae32f506c6a01`. The angular file is byte-identical
to the previously audited border_joint angular file, with SHA-256
`d1577cd3b71b0d80ec39cc1589c66bd73c2c99867d420434334721eb6d55070e`. Its full harmonic, Sturm and cyclic proofs were
nonetheless replayed in this new standalone audit.

The new radial fields are NOT the old unchanged P,H,M. They are verified afresh,
including exact reconstruction from the refit parameters. All earlier upper-bound and saturation packages were left unchanged.

No `pass` metadata or stored density interval is accepted as evidence. The
audit independently reconstructs the harmonic expansion by recurrence,
recomputes exact Sturm chains, uses a different interval geometry enclosure,
reconstructs the Fourier polynomial by a Laguerre recurrence, verifies full
signed coupling intervals, and recomputes the objective and exact rounding.
Source packages were not modified.

## Dimension-specific angular facts

The normalized dimension-46 Gegenbauer recurrence reconstructs every stored
monomial coefficient exactly. All 20 harmonic coefficients are nonnegative;
the degree is 19 and the constant coefficient is exactly 1. An independent
spherical-moment calculation also gives b0=1. Their sum and Q(1) both equal

```
B = 1153618388872687865313041135709176932237889
    / 2413832882791589326344197467275264
  = 477919742.1233786...
```

Rational Sturm variations at -1, c=5001/10000 and 1 are (10,10,9).
Together with exact endpoint signs, this proves Q<0 on [-1,c] and exactly one
root on (c,1). A 96-step exact bisection brackets that root near
0.5001000483907604. The Sturm chain for Q+1/100000 proves
Q<-1/100000 on [-1,1/2]. All nonzero remainders end in a nonzero constant.

The addition theorem for normalized Gegenbauer polynomials gives |P_l(t)|<=1
on [-1,1]. Therefore |Q(t)|<=B, and the harmonic expansion yields a PSD angular
kernel in dimension 46.

## Full cyclic geometry, scale extension and tails

Sort any triangle's side lengths as a<=b<=d, with a>=1 after scaling the
minimum distance to one. Set

```
u=b/a, t=(a^2+b^2-d^2)/(2ab), v=d/a=sqrt(1+u^2-2tu).
z1=(1-ut)/v, z2=(u-t)/v.
```

Then -1<=t<=1/2, u>=1, 2tu<=1, and z1<=z2. Dividing the cyclic Gaussian sum
by its positive first weight gives

```
Q(t) + Q(z1) exp[-alpha*a^2*(1-2tu)]
     + Q(z2) exp[-alpha*a^2*(u^2-2tu)].
```

The independent compact proof covers 0<=t<=1/2, 1<=u<=5 and 2tu<=1 at alpha=40
and a=1. It certifies **4,212 rectangles**, discards **247** only by the exact
test 2*t_lower*u_lower>1, and reaches maximum axis depth 15. All terminal
rectangles are accounted for by the exhausted subdivision stack.

The first angular term uses an exact rational Bernstein upper bound. The other
two signed angular terms use centered Arb polynomial enclosures. Unlike the
primary verifier, the square-root radicand is enclosed using its exact
monotone endpoint bounds:

```
1+u_lower^2-2*t_upper*u_lower
 <= v^2 <=
1+u_upper^2-2*t_lower*u_upper.
```

These hold because partial derivatives in u and t have signs + and - on the
compact box. Exponent gaps are clipped to nonnegative intervals only because
the feasible sorted-triangle domain has both gaps nonnegative. No feasible
point is discarded by this clipping. The final signed expression is strictly
negative on every accepted rectangle. Precision is 256 bits.

For scale extension, put A=Q(z1), C=Q(z2), p=40(1-2tu), q=40(u^2-2tu), so
0<=p<=q. The unique sign transition and z1<=z2 exclude A>0>C. For
F(T)=A exp(-pT)+C exp(-qT), the same-sign cases give
F(T)<=max(F(1),0) immediately. If A<=0<=C, the bracket
A+C exp[-(q-p)T] is nonincreasing; wherever F is positive, F'<=0. It therefore
cannot exceed its initial positive part or cross from nonpositive to positive.
This includes p=0 and p=q. Since Q(t)<0, compact negativity extends to every
T=a^2*alpha/40>=1.

The remaining domains have separate analytic bounds:

- If t<=0, both exponent gaps are at least 1.
- If u>=5 and t>=0, a positive first secondary term requires z1>c, hence its
  gap is greater than 2*c*v-1>=2*c*u-1>=4.001. The second gap is at least 24.

Dropping negative secondary terms and using |Q|<=B, both tails are dominated
by -1/100000+2B exp(-40). The exact checked rational inequality
40^30/30! > 2B*100000 proves this is strictly negative. These cases cover every
sorted separated triangle, including collinear limits.

## Exact extended-pair radial export and normalization

First, the final discovery parameters are rationalized with exactly the stored
export rule, Fraction(float(value)).limit_denominator(10^14). Independently
using normalized Laguerre three-term recurrences, the verifier reconstructs

```
P(X) = t * P48(X)/H48(0)
       + sum_(l=0)^47 c_l L_l^(n/2-1)(2X)/L_l^(n/2-1)(0).
M = c_M/B * (10001/10000).
```

Every physical coefficient agrees exactly with the frozen certificate, the
correction has degree 47, and the 100-parts-per-million upward M guard agrees
exactly. Source t is nonnegative, and every added Laguerre coefficient of
degree 32 through 47 is nonzero. This proves exact nonmembership in the old
degree-31 span: changing source t changes the nonzero degree-111 coefficient,
while a difference involving only degrees at most 47 cannot cancel that term;
with t fixed, a nonzero degree-47 correction cannot be represented at degree
31. Conversely, the old span embeds by padding sixteen zero coefficients.

The full dimension-specific Fourier-span export is also reconstructed exactly,
separately from the Fourier transform of the final physical polynomial.

The final proof does not use sampled optimizer feasibility: it starts
from these exact rational polynomials and reruns every global condition.

With X=pi|x|^2, f(x)=P(X)exp(-X), and
h(x)=exp[-(17/8)pi(|x|^2-6)]. The radius exclusion is sqrt(6).

All 112 Fourier coefficients are reconstructed independently using

```
Fourier[(pi|x|^2)^k exp(-pi|x|^2)]
 = k! L_k^(22)(pi|xi|^2) exp(-pi|xi|^2)
```

and the three-term Laguerre recurrence. They equal H exactly. H(0)>0. The exported H(0) is not assumed
to be exactly one: the final objective divides by its exact positive value.
Equivalently, divide P, H and M by H(0), preserving every homogeneous
constraint and making the Fourier origin exactly one.
Complete validated root isolation accounts for degree 111 in each polynomial,
with no forbidden real roots. Exact endpoint signs prove H>0 on [0,infinity)
and P<0 on [40,infinity). The physical endpoint is the exact rational
number 40, not 7*pi. All 111 roots of each polynomial are accounted for,
including multiplicities; a root whose imaginary enclosure contains zero
must have its real enclosure strictly below the relevant endpoint.

On the compact pair interval, the audit proves the full signed expression

```
P(X) exp[(9/8)X-(51/4)pi]
 - M (2-exp[-(17/8)(X-6pi)]) < 0.
```

All **6,664 intervals** pass, maximum depth **28**. Their exact endpoints are
stored in `pair_interval_cover.json.gz` and independently checked to tile the
whole rational enclosure [6*pi_lower,40], without gaps. The lower Machin
endpoint is strictly below 6*pi, so this covers every required point from
6*pi to 40. The independent count differs from the primary count of 6,649
because this audit encloses the full signed expression and uses a different
Machin enclosure. Neither proof relies on a sampled grid.

On each interval [a,b], the polynomial is composed exactly as an Arb polynomial
in (X-midpoint)/radius. The sum of absolute nonconstant coefficients gives a
validated symmetric enclosure for P on the whole interval. That enclosure,
an Arb interval for X, and Arb pi are substituted in the entire signed gap.
The computed upper endpoint is strictly negative on every accepted leaf.
Every bisection replaces a parent by both exact rational halves, and a separate
sorted endpoint check proves that all leaves tile the initial enclosure.

For X>=40, P<0 and 0<h<=1 imply f<0<M(2h-h^2), proving the pair inequality
on the entire remaining ray. The kernel is not truncated at X=40 or 7*pi.
Its Gaussian radial features remain nonzero and rapidly decaying at every
allowed radius.

The alternating-series 128-term Machin bounds also prove 51*pi/4>40, so the
cyclic proof applies to this same Gaussian h.

## PSD extension and exact theorem substitution

For allowed nonzero arguments define
K(x,y)=(M/B)h(x)h(y)Q(xhat dot yhat). The degree-zero nonorigin feature is
sqrt(M/B)h and its origin feature is -sqrt(MB); higher origin features vanish.
Consequently K(0,0)=MB, K(x,0)=-Mh(x), and K(x,x)=Mh(x)^2. The exact identity
(MB)(M/B)=M^2 verifies the rank-one origin Gram. Nonnegative harmonic weights
verify the remaining PSD components. For an explicit smooth extension, choose smooth radial functions chi and eta
with chi=0 near zero and chi=1 at radii at least sqrt(6), while eta=1 near
zero and eta=0 at radii at least sqrt(6). The degree-zero feature can be
sqrt(M/B)*chi*h - sqrt(MB)*eta. All higher harmonic features are multiplied
by chi and vanish at zero. These are smooth Schwartz functions, have the
required origin and allowed-radius values, and define a globally PSD kernel.
No cutoff is applied at large radius.

The pair inequality is exactly
f(x)+K(x,x)+K(x,0)+K(0,x)<=0. The full cyclic inequality is already verified.
For completeness, let p_1,...,p_N be a finite separated configuration.
For every anchor i, PSD gives

```
0 <= sum_(j,k) K(p_j-p_i,p_k-p_i).
```

Sum over i. The term with i=j=k contributes N*K(0,0). Terms with exactly two
indices equal contribute sum_(i != j) T(p_j-p_i), where
T(x)=K(x,x)+K(x,0)+K(0,x). For each unordered distinct triple, the remaining
six terms equal twice its cyclic kernel sum and are nonpositive. Therefore

```
sum_(i != j) T(p_j-p_i) >= -N*K(0,0).
```

The all-radius pair inequality f+T<=0 gives
sum_(i,j) f(p_j-p_i) <= N*(f(0)+K(0,0)). For large finite portions of a
periodic separated packing, Gaussian/Schwartz decay makes boundary errors
vanish per point and permits the periodic limit. Poisson summation and
Fourier positivity bound this same pair sum per point below by the center
density times H(0). Standard periodic approximation extends the bound to
unrestricted packings. This is the cyclic generalization of the continuum
argument, not an invocation of the stronger individual-triple-sign premise
of Theorem 1.1 of arXiv:2206.15373.

Thus the correct new numerator is P(0)+MB. No old b0 normalization remains;
the new b0 is exactly 1. The exact reconstructed objective agrees with the new
primary certificate's objective.

## Exact rounding and comparison

The density is

```
((3*pi/2)^23 / 23!) * (P(0)+MB)/H(0).
```

Substituting the rational upper Machin bound for pi proves, using rational
arithmetic alone, that this positive expression is strictly below
0.0000008777229084. The separate 512-bit Arb enclosure agrees. This is below
the preceding degree-47 frozen upper 0.0000008803541044 and the pinned table comparator
0.0000009295151556. The comparison was supplied with the task and not freshly
audited as an exhaustive literature search.

## Exact tail escape and theorem scope

The old radial template imposed P(X)<=0 for every X>=7*pi. This new witness
does not satisfy that condition. Direct rational Horner evaluation at

```
X0 = 15509/500 = 31.018
```

proves P(X0)>0, with value approximately 3.75545943990906790013e-6. The
full exact rational value is saved in `tail_escape_result.json`; it is
independently recomputed from the final physical coefficients and agrees
exactly with the supplied primary witness. The rational upper Machin bound
also proves 7*pi<X0<40. Since exp(-X0)>0, the actual radial function f is
strictly positive there too. This is genuine nonmembership in the old
forced-negative-tail cone, not just a different parameterization or a changed
claim about the same function.

The correct theorem needs the combined inequality f(x)+T(x)<=0 for every
allowed radius. It does not separately require f<=0 after squared radius 7.
The noncompact Gaussian kernel provides a strictly positive allowance
M(2h-h^2), which absorbs the positive lobe. The full interval proof includes
X0, every neighboring point, and all other points up to 40; the complete
negative ray covers everything beyond 40. An earlier rejected search lobe
or a list of numerical peak samples is not used as evidence.

Consequently the older degree-31 saturation lower certificates, whose
necessary rows included physical negativity beyond 7*pi, do not apply to
this enlarged tail template. The present witness also uses the degree-47
correction rather than the degree-31 correction. It therefore escapes both
restrictions. This audit proves a valid upper-bound witness and its tail
scope escape; it does not prove optimality in the new family or invalidate
any correctly scoped older lower certificate. Angular Q and the positive
Gaussian rate are unchanged. Harmonic-dependent radial matrices are not used.

## Replay

Run `python verify_all.py` with Python 3.11+ and python-flint 0.9.0; do not use
Python -O. The runtime used here is
`research/density_solid_cyclic_20261111/.venv/bin/python`. No optimizer or
network is required. The input snapshots, fresh scripts, detailed JSON results
and aggregate `audit_result.json` are sufficient to replay this audit.
