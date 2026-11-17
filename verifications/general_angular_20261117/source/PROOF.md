# Signed logarithmic radial profiles

These packages instantiate the same unrestricted three-point density theorem
used by the earlier cyclic certificates in this repository
([Cohn–de Laat–Salmon, Section 2.2, equation (2.4), and Theorem 2.4](https://arxiv.org/abs/2206.15373)).
They change the admissible functions and supply complete checks of their
hypotheses; they do not introduce a new density theorem or require external
neighbor counts. The periodic PSD summation uses the cyclic sum condition;
it does not substitute that condition for the individual-sign premise of the
printed Theorem 1.1. Absolute convergence permits grouping into triangle sums;
periodic approximation gives unrestricted packing density.
Discovery uses floating-point optimization. Certification uses
the frozen rational data and outward FLINT/Arb arithmetic only.

## Functions, normalization and positivity

Let the minimum squared distance be `s>0`, put `X=pi |x|^2` and
`v=|x|^2/s`, and use

```
f(x) = P(X) exp(-X),       Fourier(f)(y) = H(pi |y|^2) exp(-pi |y|^2),
h(v) = v^(-a) exp[-b(v-1)],       a+b>0, b>0.
```

The frozen packages use integer `a` and `q=1`. In particular `a=-96,b=120`
gives `h(v)=v^96 exp[-120(v-1)]`. On `v>=1` it is positive and decreasing,
with `h(1)=1`. Its logarithmic slope `v*phi'(v)=a+b*v` is positive and
increasing. The feature construction below supplies a smooth extension at
excluded radii; it does not require `h` itself to be a positive definite
radial function.

Write the angular polynomial as `Q(t)=sum c_l G_l^(n)(t)`, with normalized
Gegenbauer polynomials `G_l^(n)(1)=1`, `c_0=1`, every `c_l>=0`, and
`B=Q(1)=sum c_l`. The nonorigin kernel is

```
K(x,y) = C/B^2 h(|x|^2/s) h(|y|^2/s) Q(<x,y>/(|x||y|)),
K(0,x) = -C/B h(|x|^2/s),          K(0,0)=C>0.
```

The spherical addition formula proves positivity in each nonconstant harmonic.
The constant harmonic, together with the origin, has coefficient matrix
`[[C,-C/B],[-C/B,C/B^2]]`, whose determinant is exactly zero and whose
diagonal entries are positive. This proves the full origin-augmented PSD
condition, including configurations with repeated points. In particular
`K(x,x)=C/B*h(v)^2`.

For a global Schwartz kernel, choose a smooth radial cutoff `chi`, zero in
a neighborhood of the origin and one for `|x|>=sqrt(s)`. Use the degree-zero
feature `sqrt(C)/B*chi*h - sqrt(C)*(1-chi)` and multiply every positive-degree
spherical feature by `chi*h`. The origin value of that feature is
`-sqrt(C)`, and allowed nonorigin values are `sqrt(C)/B*h`. This preserves
the displayed kernel, origin and packing-domain values. The cutoff removes
the angular singularity at zero, and Gaussian decay makes all features
Schwartz. It proves global PSD as well as smoothness.

The primary Fourier check uses the coefficient identity

```
H_j = sum_(k>=j) P_k k! (-1)^j (n/2+j)_(k-j) / ((k-j)! j!).
```

The independent path reconstructs the same transform by the Laguerre
recurrence. Both require exact rational equality, `H(0)>0`, and account for
every complex root, with multiplicities. Roots whose imaginary enclosure
contains zero must have their entire real enclosure strictly below zero.
This and `H(0)>0` prove `H>0` on the entire nonnegative ray. Similarly all
possibly real roots of `P` must lie below the rational cutoff, and `P` is
strictly negative at that cutoff. These are complete ray checks, not checks
on a finite sample or a finite truncation.

## Pair condition

For all `X>=pi*s` the required inequality is

```
P(X) exp(-X) <= C/B * (2h(X/(pi*s))-h(X/(pi*s))^2).
```

Beyond the cutoff, `P<0` and `0<h<=1` suffice. On the whole finite interval,
the primary path bounds `P` by its translated polynomial coefficients and
uses the concavity of `2h-h^2` to obtain a lower bound from endpoint
enclosures. The independent path instead encloses the signed difference,
uses integer powers to evaluate the profile, and verifies that all accepted
rational subintervals form an exact gap-free partition. Rational Machin
bounds enclose `pi`; monotonicity is checked on the slightly enlarged
interval caused by those bounds. Neither path accepts a numerical minimum.

## Triangle condition and all dilations

Sort triangle side lengths and write them as
`sqrt(s*T)*(1,u,w)`, with `T>=1`, `u>=1`,
`t in [-1,1/2]`, `2*t*u<=1`, and `w^2=1+u^2-2tu>=u^2`.
The secondary angles are
`alpha=(1-ut)/w` and `beta=(u-t)/w`, with `beta>=alpha`.
After division by a positive common factor the cyclic inequality becomes

```
Q(t) + (w^2/u^2)^(-a) exp[-b*T*(w^2-u^2)] Q(alpha)
     + (w^2)^(-a) exp[-b*T*(w^2-1)] Q(beta) <= 0.
```

The exact angular checks prove `Q<=-delta` on `[-1,1/2]`, no zero up to
`c=5001/10000`, and exactly one zero on `(c,1)`. Its sign therefore changes
only once, from negative to positive. The powers of `T` cancel from the
two profile ratios. Put `g1=w^2-u^2` and `g2=w^2-1`; then `0<=g1<=g2`.
For `T>=1`, if both secondary terms are positive, both decrease. If both
are negative, they already suffice with `Q(t)<0`. In the remaining case
`Q(alpha)<=0<Q(beta)`, factoring `exp[-b*(T-1)*g1]` shows that the positive
term decreases at least as quickly as the negative one. The sum is at most
the larger of its value at `T=1` and zero, so the verified `T=1` inequality
and `Q(t)<0` imply the inequality at every dilation. The opposite sign
ordering is excluded by the unique transition and `beta>=alpha`.

For `1<=u<=U`, both implementations cover the complete compact domain by
outward interval subdivision. Rectangles are discarded only when their
entirety violates the sorting condition. The primary implementation uses
exact Bernstein upper bounds for the first angular term and interval
logarithms for profile ratios. The independent implementation reconstructs
harmonics by a closed formula, proves the negative cap by exact Bernstein
subdivision, and uses endpoint monotonicity of
`R(g,v)=(1+g/v)^(-a)*exp(-b*g)`: it decreases in `g`, and its derivative
in `v` has the sign of `a`. It also checks the exact spherical mean is one.
Every accepted rectangle has strictly negative upper bound.

For the infinite tail `u>=U`, let `p=max(0,-a)`. A first secondary term can
be positive only if `alpha>c`, which implies `w^2-u^2>2cu-1`.
Also `w/u<=2`. The second secondary term is bounded using `w^2>=u^2`.
Since `|Q|<=B` by the normalized addition formula, their total positive
contribution is at most

```
B * [4^p exp(-b*(2cU-1)) + U^(2p) exp(-b*(U^2-1))].
```

The second summand decreases for `U^2>p/b`. The verifiers check this
condition and prove the bound is below `delta/2`. Primary verification also
uses a single rational Taylor term as a lower bound on the exponential.
Thus the tail is strictly negative, uniformly in the dilation. All packages
record the actual `U` used; the selected signed profiles use `U=8`.

## Strict density and implementation scope

The conclusion is

```
density <= (P(0)+C)/H(0) * (pi*s/4)^(n/2) / Gamma(n/2+1)
        < the registered rational decimal.
```

The last comparison uses rational Machin upper bounds for `pi`; in odd
dimensions the positive quantities are squared and the half-integer gamma
factor is eliminated exactly. The displayed decimal is not a rounded Arb
midpoint. Assertions must be enabled; every executable verifier rejects
`python -O`. Each package's manifest binds its rational inputs, wrapper,
and all shared verification code. `verify.py` runs all four gates and only
then writes a complete successful replay result.

The two paths have different Fourier reconstructions, angular reconstructions,
cap proofs and profile enclosures. They share Python, FLINT/Arb, polynomial
root isolation and the same underlying mathematical theorem. They are
independent implementations of the specified gates, not a separate arithmetic
backend, a formally checked proof, or an external audit.
