# Rigorous radial and joint saturation bounds in dimensions 46 and 47

**PASS: six genuine-transcendental radial dual lower certificates, plus replay
of four all-degree angular lower certificates.**

These are lower bounds on the best *certificate objective within a specified
family*. They are not lower bounds on sphere-packing density and do not prove
that a packing is optimal.

| Dimension | Fixed saved-Q radial-family floor | Certified current upper witness | Remaining relative improvement |
|---|---:|---:|---:|
| 46 | >0.0000008868479138 | <0.0000008868492132 | <1.466 ppm |
| 47 | >0.0000006010075028 | <0.0000006010083731 | <1.449 ppm |

The current upper witnesses are the separately verified `dim46_border_refit`
and `dim47_border_refit` packages. The zero-margin duals prove that the frozen
radial span is essentially saturated at each saved angular polynomial.

Combining radial duals at the rigorous angular lower values with the replayed
all-degree angular proofs gives:

| Dimension | Rigorous angular floor L_B | Joint-family density-objective floor | Remaining improvement |
|---|---:|---:|---:|
| 46 | 475115741 | >0.0000008830531394 | <0.428041% |
| 47 | 650216751 | >0.0000005981706612 | <0.472159% |

“Joint” means the scalar separable Gaussian angular family at the fixed damping
alpha=51*pi/4, allowing all finite harmonic degrees, together with the radial
span and sign template described below. It does not include matrix-valued
multi-Gaussian kernels, a different damping, a larger radial span, or removing
the physical-tail negativity requirement.

The stronger rate-independent composition permits **every positive scalar
Gaussian rate** and all finite angular harmonic degrees:

| Dimension | Rate-independent angular floor | Joint all-rate density-objective floor | Remaining improvement |
|---|---:|---:|---:|
| 46 | 473143566 | >0.0000008803784166 | <0.72964% |
| 47 | 646517390 | >0.0000005956056153 | <0.89895% |

The radial span, minimum squared separation 6, and physical-tail negativity
beyond squared radius 7 remain fixed. The bound does not extend to higher
Laguerre correction degrees or multi-feature kernels.

## Exact family being bounded

Let P_s be the frozen degree-111 source-48 physical polynomial divided by its
original positive Fourier constant coefficient. Let H_s be its exact Fourier
polynomial in the target dimension n. Put a=n/2 and

```
ell_k(X) = L_k^(a-1)(2X) / L_k^(a-1)(0),  k=0,...,31.
P(X) = t P_s(X) + sum v_k ell_k(X).
H(X) = t H_s(X) + sum (-1)^k v_k ell_k(X).
```

The constraints are t>=0, free v_0,...,v_31, M>=0, and H(0)=1. The scalar
Gaussian feature is h(u)=exp[-(17/8)*pi*(u-6)]. The sign template requires:

1. H(X)>=0 for every X>=0.
2. P(X)<=0 for every X>=7*pi.
3. P(pi*u)*exp(-pi*u)<=M*(2h(u)-h(u)^2) for every 6<=u<=7.

The objective is P(0)+BM, multiplied by the volume
V_n=pi^(n/2)*(3/2)^(n/2)/Gamma(n/2+1). Existing upper witnesses are in this
family after positive Fourier normalization. The nonnegative t restriction is
part of the family, not an assertion about arbitrary radial functions.

The physical-tail condition is important. The proof uses 18 physical-tail
rows. It does NOT give the same lower bound for the larger cone that permits
positive P beyond 7*pi while imposing only the coupled inequality there.

## Necessary sample rows and elimination

For a fixed exact positive B, write C=BM>=0. Each sampled inequality has the
form A_t*t+sum A_k*v_k+A_C*C<=0, with no artificial margin. The actual rows are:

```
Pair, X=pi*u:
  (P_s(X), ell_0(X),...,ell_31(X),
   -exp(X)*(2h(u)-h(u)^2)/B).
Physical:
  (P_s(X), ell_0(X),...,ell_31(X),0), X>=7*pi.
Fourier:
  (-H_s(X), -ell_0(X),..., -(-1)^31 ell_31(X),0), X>=0.
```

The verifier checks the type of every row and its exact rational domain
parameters. Pair points are rational u in [6,7]. Physical points are either
exact rational X with a validated X>=7*pi comparison or exact rational
interpolations between 7*pi and 200. Fourier points are rational X>=0.
Positive powers-of-two row multipliers improve conditioning and do not change
the halfspaces. All polynomial, pi and exponential values are reevaluated with
Arb; stored floating matrix entries are not used.

Since ell_k(0)=1, Fourier normalization eliminates

```
v_0 = 1-H_s(0)*t-sum_(k=1)^31 (-1)^k*v_k.
```

For z=(t,v_1,...,v_31,C), the transformed row and right side are

```
Atilde = (A_t-A_0*H_s(0),
          A_k-A_0*(-1)^k for k=1,...,31,
          A_C),
Atilde*z <= -A_0.
```

The objective is 1+c^T*z, where

```
c = (P_s(0)-H_s(0), 1-(-1)^k for k=1,...,31, 1).
```

The source Fourier transform is rebuilt exactly by the Gaussian/Laguerre
recurrence. All 32 Laguerre Fourier eigenvalue identities and ell_k(0)=1 are
also checked in exact rational polynomial arithmetic.

## Genuine dual certificate

The zero-margin floating LP is used ONLY to choose 33 rows. Each of the six
selected systems has one contact pair row u=6, 18 physical rows and 14 Fourier
rows. Its t and C primal coordinates are positive, so no boundary slack was
needed.

Let D be the 33-by-33 matrix whose columns are the selected true Atilde rows.
Define y to be the unique exact real solution

```
D*y = -c.
```

This is an inverse-defined certificate with transcendental entries, not a
rational approximation to the dual. At 1024-bit precision, the verifier:

- rigorously encloses every true pi/exponential row;
- proves the determinant enclosure excludes zero;
- uses Arb's validated preconditioned solve, never its approximate solve;
- proves every true y_i is strictly positive;
- encloses 1+sum y_i*A_0i and proves the safely rounded density floor below it.

The defining equation establishes exact stationarity, including zero reduced
costs for the two nonnegative variables t and C. Floating residuals are not a
certificate. No claim is made that every point inside the displayed weight
intervals satisfies stationarity: the exact uniquely defined solution does.

For any feasible z, summing the necessary inequalities with those weights gives

```
-c^T*z <= -sum y_i*A_0i,
1+c^T*z >= 1+sum y_i*A_0i.
```

Multiplying by positive V_n proves the reported density-objective lower bound.
The successful finite system is sufficient; there is no need to certify that
the floating LP found the optimum or to verify unselected rows. The primal
continuum feasible set is contained in the finite sampled relaxation.

All selected pair rows happen to be exactly u=6, where h=1. Thus these radial
lower certificates, for a specified B and the same physical/Fourier template,
do not otherwise depend on kappa. The JOINT angular certificate still uses
the fixed alpha=51*pi/4.

The separate all-rate replay uses cap-only angular certificates and explicitly
checks that every angular dual row is a cap row, and every radial pair row is
at contact. These facts remove all rate dependence. Cap negativity on [-1,1/2]
is necessary for every positive-rate scalar Gaussian cyclic kernel: for t<1/2,
send the shortest-side scale to infinity with u=1 so both secondary
exponentials vanish; at t=1/2 the equilateral cyclic sum is 3Q(1/2). Therefore
the cap-only angular lower bounds apply uniformly over all positive rates.

## Joint angular-radial bridge

The angular worker supplied standalone certificates proving B>L_B for all
finite harmonic degrees in the scalar Gaussian cyclic angular family. The
aggregate replay actually executes the pinned standalone angular verifier,
including its finite reduced-cost checks and exact integral-moment tail. It
does not trust copied `pass` metadata.

Radial feasibility expressed in (P,M) is independent of B. Since M>=0,

```
P(0)+B*M >= P(0)+L_B*M.
```

The separately validated radial dual at B=L_B therefore lower-bounds every
feasible joint member with B>=L_B. The angular extension removes the degree-19
restriction, but does not remove the radial span, physical-tail gate, scalar
separability, or fixed angular damping restrictions.

The exact conservative remaining-improvement fractions are stored in
`saturation_result.json`. They compare the safe lower values against the
separately proved rounded upper witnesses, so the reported percentages are
conservative.

## Replay and artifacts

Run `python verify_all.py` with Python 3.11+ and python-flint 0.9.0. Do not use
Python -O. This complete lower-certificate replay needs no optimizer, NumPy,
SciPy, network, or files outside this folder.

- `radial_model.py`: exact source/basis identities and genuine necessary rows.
- `dual46/47_fixed_selected.json`: selected rows at the saved angular B.
- `dual46/47_joint_selected.json`: selected rows at the rigorous lower L_B.
- `dual46/47_allrate_selected.json`: selected rows at the cap-only angular floors.
- `verify_dual.py`: validated inverse-defined dual proof.
- `verify_angular_certificate.py` and `inputs/angular*_lower_certificate.json`:
  pinned standalone all-degree angular proofs supplied by the angular worker.
- `verify_all.py`: actual joint replay and conservative gap calculations.
- `*_verified.json`, `angular*_replayed.json`, `saturation_result.json`: detailed
  verified intervals and aggregate results.
- `INPUT_SHA256.json`: frozen input and selected-row hashes.
- `test_negative_controls.py`: rejection tests for invalid row types, domains,
  signs, singular active sets, nonpositive B and disabled assertions.
- `discover_dual.py`: optional row-selection experiment, requiring NumPy/SciPy;
  no output from its floating optimization is trusted as a proof.

All work is local to this folder. The earlier upper-bound and border_joint
packages were not changed.
