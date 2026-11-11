# Independent audit: new 47D angular kernel with unchanged radial witness

**PASS. The exact frozen construction proves Delta_47 < 0.0000006117584752.**

The independently computed density enclosure is centered at
6.1175847516131160009888892627056047388644491389730877895442236214498733399560517560699327822704351260638803523947521919509350086577234653728351736132042e-7,
with radius below 4.3e-159 at 512-bit precision.

This is an independent computer-assisted audit of the supplied construction,
not a proof-assistant formalization, external peer review, or literature
priority determination. The pinned public-table comparison is separate from
the certificate.

## Inputs and independence

The audited frozen inputs are:

- New angular polynomial and complete primary certificate from
  `research/dim47_border_joint_ready_20261111/`.
- Old physical polynomial from
  `research/dim47_transfer_ready_20261111/continuum_three_point/adaptive_cyclic40_radial111.json`.
- Old dimension-47 Fourier polynomial and M from
  `research/dim47_transfer_ready_20261111/certificate47_transfer.json`.

Copies are frozen in `inputs/`; their SHA-256 values are pinned in
`verify_all.py`. The new P,H,M fields are compared entry-for-entry with these
old source fields. Their shared canonical SHA-256 is
`942f9cf4152a8a5f84519e8d5025360490ffb69e35919404b1b8550ff2fbb476`.

No `pass` metadata or stored density interval is accepted as evidence. The
audit independently reconstructs the harmonic expansion by recurrence,
recomputes exact Sturm chains, uses a different interval geometry enclosure,
reconstructs the Fourier polynomial by a Laguerre recurrence, verifies full
signed coupling intervals, and recomputes the objective and exact rounding.
Source packages were not modified.

## Dimension-specific angular facts

The normalized dimension-47 Gegenbauer recurrence reconstructs every stored
monomial coefficient exactly. All 20 harmonic coefficients are nonnegative;
the degree is 19 and the constant coefficient is exactly 1. An independent
spherical-moment calculation also gives b0=1. Their sum and Q(1) both equal

```
B = 28480256271912302890787749475800706759447932616462987363
    / 43526477349886437471741192702003770547071090688
  = 654320266.7878339...
```

Rational Sturm variations at -1, c=5001/10000 and 1 are (10,10,9).
Together with exact endpoint signs, this proves Q<0 on [-1,c] and exactly one
root on (c,1). A 96-step exact bisection brackets that root near
0.5001000379330464. The Sturm chain for Q+1/100000 proves
Q<-1/100000 on [-1,1/2]. All nonzero remainders end in a nonzero constant.

The addition theorem for normalized Gegenbauer polynomials gives |P_l(t)|<=1
on [-1,1]. Therefore |Q(t)|<=B, and the harmonic expansion yields a PSD angular
kernel in dimension 47.

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
and a=1. It certifies **4,860 rectangles**, discards **251** only by the exact
test 2*t_lower*u_lower>1, and reaches maximum axis depth 16. All terminal
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

## Unchanged radial witness

With X=pi|x|^2, f(x)=P(X)exp(-X), and
h(x)=exp[-(17/8)pi(|x|^2-6)]. The radius exclusion is sqrt(6).

All 112 Fourier coefficients are reconstructed independently using

```
Fourier[(pi|x|^2)^k exp(-pi|x|^2)]
 = k! L_k^(45/2)(pi|xi|^2) exp(-pi|xi|^2)
```

and the three-term Laguerre recurrence. They equal H exactly. H(0)>0.
Complete validated root isolation accounts for degree 111 in each polynomial,
with no forbidden real roots. Exact endpoint signs prove H>0 on [0,infinity)
and P<0 beyond 7*pi. The physical root exclusion is slightly stronger, starting
at seven times the rational lower Machin bound for pi.

On the compact pair interval, the audit proves the full signed expression

```
P(X) exp[(9/8)X-(51/4)pi]
 - M (2-exp[-(17/8)(X-6pi)]) < 0.
```

All **33 intervals** pass, maximum depth **20**. Their exact endpoints are
stored in `pair_interval_cover.json` and independently checked to tile the
whole rational enclosure of [6pi,7pi], without gaps. For |x|^2>=7, f<0 and
0<h<=1 prove f<=M(2h-h^2) automatically.

The alternating-series 128-term Machin bounds also prove 51*pi/4>40, so the
cyclic proof applies to this unchanged h.

## PSD extension and exact theorem substitution

For allowed nonzero arguments define
K(x,y)=(M/B)h(x)h(y)Q(xhat dot yhat). The degree-zero nonorigin feature is
sqrt(M/B)h and its origin feature is -sqrt(MB); higher origin features vanish.
Consequently K(0,0)=MB, K(x,0)=-Mh(x), and K(x,x)=Mh(x)^2. The exact identity
(MB)(M/B)=M^2 verifies the rank-one origin Gram. Nonnegative harmonic weights
verify the remaining PSD components. Smooth cutoffs, zero near the origin and
one at radii at least sqrt(6), extend the features smoothly and rapidly
decreasingly without changing required values.

The pair inequality is exactly
f(x)+K(x,x)+K(x,0)+K(0,x)<=0. The full cyclic inequality is already verified.
For completeness, the direct periodic PSD summation gives the usual cyclic
continuum conclusion: sum each anchored PSD kernel over packing points, group
all distinct triples into cyclic sums, and leave the repeated-index diagonal
and origin-cross terms. Nonpositive distinct-triple sums imply that the
negative of those repeated-index terms is at most N*K(0,0). Combining with
the pair inequality bounds the f pair sum above by N*(f(0)+K(0,0)). Poisson
summation and Fourier positivity bound it below by N*density_of_centers*H(0).
Gaussian decay makes all sums absolutely convergent, and periodic
approximation extends the bound to unrestricted packings.

Thus the correct new numerator is P(0)+MB. No old b0 normalization remains;
the new b0 is exactly 1. The exact reconstructed objective agrees with the new
primary certificate's objective.

## Exact rounding and comparison

The density is

```
((3*pi/2)^(47/2) / Gamma(49/2)) * (P(0)+MB)/H(0).
```

The odd-dimensional volume identity is checked exactly:

```
V_47^2 = 4*pi^46*6^47*(24!)^2/(48!)^2.
```

It follows from Gamma(49/2)=48!*sqrt(pi)/(4^24*24!). The verifier squares the
positive objective and volume, substitutes the rational upper Machin bound
for pi, and proves the result is strictly below the square of
0.0000006117584752. Positivity makes this a rigorous unsquared upper bound.
The separate 512-bit Arb enclosure agrees. This is below the old frozen upper
0.000000616879891 and the pinned table comparator 0.0000006490757338.
The comparison was supplied with the task and not freshly
audited as an exhaustive literature search.

## Replay

Run `python verify_all.py` with Python 3.11+ and python-flint 0.9.0; do not use
Python -O. The runtime used here is
`research/density_solid_cyclic_20261111/.venv/bin/python`. No optimizer or
network is required. The input snapshots, fresh scripts, detailed JSON results
and aggregate `audit_result.json` are sufficient to replay this audit.
