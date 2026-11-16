# Global Schur coupling with a radial lattice certificate

This method has no near-frame arithmetic or prior density assumption. Numerical
search is not a proof; only a successfully replayed rational certificate is.

Normalize the lattice minimum to1. Its m=2^n-1 nonzero classes modulo2Λ each
contain at most one opposite pair of squared norm<2. Let A be classes whose
short pair has squared norm in [1,U), where U=the first partition endpoint.
Write N_A for their number and N_j for the number with a short pair in bin j.
Classes with no short pair are also allowed.

## Angular counting

Choose a projective Delsarte polynomial f, normalized Gegenbauer constant1,
with all other even Gegenbauer coefficients nonnegative. Require f(x)<=0 on
[0,U-1]. For distinct A classes, any positive contribution f(|cos angle|)
forces a signed sum w of their representatives to have squared norm<2:
||u±v||²<=2U-2|cos angle|<2. Thus w is the unique short pair of its nonzero
sum class. This is an actual vector sum, not merely a congruence.

For a target bin [a,b), use its closed envelope [a,b] for the following safe cosine bounds:
[max(0,(2-b)/(2U)), min(1-1/(2U),(2U-a)/2)].
Let M_j>=0 bound f on this interval. Empty/no-short target classes contribute
nothing positive.

Fix an orientation of w. Each ordered pair of source classes contributing to
this target has exactly one decomposition w=u+v with that ordered labeling.
The centered points y=u-w/2 occur in opposite pairs, are separated by at
least1, and have squared radii in[1-t/4,U-t/4]. Hence their normalized
absolute inner products are at most c=1-1/[2(U-a/4)]. The monotonicity used
here follows from the squared-radius spread<=U-1<1.

Also |y·w|<= (U-1)/2. Their normalized axial coordinates have square at most
eta²=(U-1)²/[a(4-a)]. Projection to w-perpendicular therefore has absolute
normalized inner products at most c'=(c+eta²)/(1-eta²), in dimension d=n-1.
The code is projective: the point count C_j is twice an upper bound on the
number of opposite pairs. Use the minimum of applicable exact bounds:

- d(1-c'²)/(1-d c'²), if d c'²<1;
- d(d+2)(1-c'²)/(3-(d+2)c'²), if (d+2)c'²<3.

The latter follows from the projective polynomial x²(x²-c'²); its normalized
Gegenbauer coefficients are nonnegative. Floors may be taken before doubling.
The17/19 certificates use these two formulas. Higher-dimensional certificates may also use the independently checked sextic polynomial x^4(x²-c'²), yielding d(d+2)(d+4)(1-c'²)/[3(5-(d+4)c'²)] when (d+4)c'²<5. The verifier checks every applicable denominator before taking the minimum.

Consequently N_A²<=f(1)N_A+sum_j C_j M_j N_j. For every positive rational N0,
N_A²>=2N0 N_A-N0² gives the necessary linear inequality
sum_j r_j N_j<=N0, where
r_j=(2-f(1)/N0)1_{j=0}-C_j M_j/N0.
Negative r_j are allowed; M_j and C_j are nonnegative upper bounds.

## Combination with the radial auxiliary function

Use disjoint half-open squared-norm bins [E_j,E_{j+1}) for class accounting. The short part is [1,2), the long part is [2,E_end), and the nonpositive tail covers E_end itself. Actual short and long squared norms are strictly below their respective bin upper endpoints. Thus the conditional zero-cap test on endpoint bounds a+S<=4 is valid: the actual squared norms have sum strictly less than 4. At actual equality 4, distinct orthogonal pairs can occur, so the test would not be valid for closed-bin counts. Angular and radial polynomial inequalities are checked on the closed containing intervals, which is conservative. Likewise a source-prefix window denotes [1,E_prefix); the recorded closed source_interval is an envelope, not a second assignment of boundary classes.

Let h(r)=p(tau r²)exp(-tau r²/2), with nonnegative Fourier transform and
h<=0 outside the finite partition. Exact interval envelopes bound each class
profile by its ordinary conditioned/unconditioned projective caps.
For nonnegative rational multipliers nu_l, subtract sum_l nu_l r_lj from
class j's payoff, and let kappa bound every resulting class profile, including
empty classes. The global sum then gives
rho*hhat(0)<=h(0)+2m*kappa+2*sum_l nu_l*N0_l.
The zero mod2 class has no positive contribution because its nonzero vectors
have squared norm>=4, beyond the tail cutoff.

The default `python run_all.py` is a complete Python3.10+ standard-library proof. It reconstructs the Fourier polynomial and proves both global signs by exact Bernstein subdivision after rational compactification of the half-lines. It also checks every radial interval, continuous angular inequality, projective cap, class budget and exact density normalization. The optional `python verify_sympy.py certificate.json` is a separate SymPy-based implementation using root counts for global signs; public `python-flint` is optional acceleration. No SymPy dependency is needed by the default standard-library verifier.

Odd-dimensional normalization uses a rational lower bound for pi and an exact squared-density comparison. No floating sign tolerances are used. Do not use `python -O`.

A successful replay validates that specific certificate. It does not establish
optimality, universality of sampled experiments, or priority over other work.

Some exact rational strings exceed Python3.11+ default4300-digit conversion limit. The verifier explicitly disables that limit for these supplied certificate files; all arithmetic remains exact. This is a serialization accommodation, not a floating approximation.
