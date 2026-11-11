# Long-edge lens cuts funded by the physical tail

## Outcome

This is one concrete extension of the current R=3/2 cyclic certificate. It
requires no new spherical-code input, count-radius PSD feature, or objective
term. It permits positive K on selected one-long-edge triangles, paid for by
the already negative physical function f at their long base edge. It is valid
in every dimension, including 5, 6 and 7. No improved density is claimed here.

The extension is distinct from the previously inactive all-close-edge cuts:
their allowances vanish in the one-long-edge domain; these allowances are
supported there. The root PSD matrices and compact cyclic conditions stay
unchanged.

## 1. Exact common-neighbor lens caps

Fix a packing with minimum center separation 1. Fix two centers x,y separated
by d, and let z_1,...,z_M be centers with both endpoint distances at most
r=27/25. Write m=(x+y)/2. The parallelogram identity gives

    |z_i-m|^2 = (|z_i-x|^2+|z_i-y|^2)/2 - d^2/4
               <= r^2-d^2/4 = rho^2.

The exact centroid identity, with zbar the average of the z_i, gives

    M(M-1)/2 <= sum_(i<j) |z_i-z_j|^2
              = M sum_i |z_i-zbar|^2
              <= M sum_i |z_i-m|^2 <= M^2 rho^2.

When rho^2<1/2 this implies

    M <= 1/(1-2rho^2) = 2/(d^2-4r^2+2).

Consequently the following conservative, closed-endpoint bins are exact:

| Base length d | Common-neighbor cap B(d) | rho^2 at lower edge | Real count upper |
|---|---:|---:|---:|
| [9/5,46/25) | 3 | 891/2500 | 1250/359 < 4 |
| [46/25,48/25) | 2 | 8/25 | 25/9 < 3 |
| [48/25,54/25] | 1 | 153/625 | 625/319 < 2 |

The exact numerical slack makes these valid even when both short distances
equal r and the long distance equals a bin's lower endpoint. Beyond 2r=54/25
there are no common neighbors. At d=2r one midpoint common neighbor is possible,
so the final bin includes that endpoint.

No projection, angular monotonicity assumption, dimension-dependent kissing
bound, or orientation factor is used. The inequalities also hold with strict
short-distance cutoffs. The closed convention is used below.

The sharper cap floor(2/(d^2-4r^2+2)) can be used wherever its denominator is
positive, but transitions at exact integers require care: with closed short
balls the integer itself is permitted. The conservative rational bins avoid
this issue. `check_lens_and_tail.py` checks all displayed rational arithmetic.

## 2. General weighted pair-triangle inequality

Let J=[9/5,54/25] and choose any nonnegative function eta(d) supported on J.
It depends ONLY on the base-edge length d. For an unordered pair e={x,y}, let
M_e count common neighbors at distance at most r from both endpoints. Then

    sum_e eta(|e|) M_e <= sum_e eta(|e|) B(|e|).

On the left each unordered triangle is counted once for each eligible base
edge. Since every eligible base has length at least 9/5>r, and its other two
sides must have length at most r, a triangle has at most ONE eligible base.
Thus there is no factor 2 or 3 in this incidence inequality.

The same statement with sum over eligible bases is valid without uniqueness;
uniqueness merely makes the implementation particularly simple here.

For a finite periodic quotient, or a sufficiently large finite patch with
the standard boundary limiting argument, let C_K(a,b,c) be the sum of K at
the three vertices of an unordered triangle. Replace the old C_K<=0 by

    C_K(a,b,c) <= eta(c),

when c is the unique eligible long base and a,b<=r. Everywhere else retain
C_K<=0. Summing gives

    sum_triangles C_K <= sum_long_pairs eta(d) M_e
                     <= sum_long_pairs eta(d) B(d).

The root-PSD expansion and Fourier sum have the normalization

    S_f + S_K = N[f(0)+K(0,0)]
              + 2 sum_pairs [f(d)+2K(x-y,0)+K(x-y,x-y)]
              + 2 sum_triangles C_K.

Thus the complete modified pair condition is the original pair expression
plus eta(d)B(d), bounded above by the same count allowance. Both the pair and
triangle correction have multiplicity two. The objective and published count
terms do not change.

## 3. Tail-funded specialization

Retain the original global sign f(d)<=0 for d>=R=3/2. Set

    eta(d) = -f(d)/B(d), for d in J,
             0, otherwise.

This is nonnegative, depends only on d, and requires no sign assumption on K.
For d in J, all nonzero root-kernel pair terms vanish because d>R, and all
published count allowances vanish. The modified pair expression is exactly

    f(d) + B(d) eta(d) = 0.

Elsewhere the original pair condition remains sufficient. Since an eligible
triangle has exactly one edge outside R, C_K is just the surviving kernel at
the opposite vertex. The new far condition is therefore

    K(a,b,c) + f(c)/B(c) <= 0,

for a,b in [1,27/25] and c in J, and K(a,b,c)<=0 for every other one-long-edge
configuration. In the notation of the discovery code, u=a^2, v=b^2, w=c^2,
and the expression is kernel_row(u,v,w) @ vector + f(sqrt(w))/B(sqrt(w)).
Use the same positive row scaling on BOTH terms.

The compact cyclic inequalities are literally unchanged since J starts above
R. No condition on the sign of K is needed other than the modified cyclic/far
conditions just stated. A smaller allowance theta*(-f)/B, 0<=theta<=1, is also
valid. If multiple lens families are combined, their pair budgets must obey
sum_j B_j(d) eta_j(d)<=-f(d); reusing the full tail independently is invalid.

The original density formula, Schur saturation, Fourier nonnegativity, and
unrestricted-packing passage remain unchanged. If f has Gaussian repairs,
use the complete repaired physical f in this formula and its global tail sign.

## 4. Boundaries and implementation hazards

- eta is zero outside J. At a base-bin boundary use the cap assigned in the
  table; both neighboring cap inequalities are valid at that boundary.
- The closed short-radius convention is valid in the counting proof. A
  continuous K must also satisfy limiting inequalities from the unrelaxed
  side a>r or b>r. A verifier that checks both one-sided combinations is safe.
  Do not omit the stronger boundary face merely because the indicator takes
  the closed-side value.
- The pair tail sign remains globally required, including zeros. The budget
  vanishes wherever f vanishes. The pair inequality f+B eta=0 should not be
  burdened with an artificial strict numerical margin.
- Existing saved far cuts, antipodal conditions, and local Bernstein patches
  must be reconsidered in the eligible region. Keeping their old K<=0 right
  sides is sound but can entirely suppress the improvement being tested.
- The previous purely polynomial far verifier does not automatically certify
  K+f(c)/B: f contains exp(-pi*c^2). Use rigorous polynomial bounds for that
  analytic term, or interval/Taylor subdivision. Alternatively use the
  certified rational step allowances below, keeping far checks polynomial.
  Numerical sampling is not a certificate. The compact cyclic verifier is
  unaffected.
- In the existing local far box coordinates, a=1+x/2, b=1+y/2 and
  t=-1+16z/9. The whole patch x,y in [.04,.16], z in [.12,.19] has a,b<=r
  and c about 1.860..2.042, so the cap is 2 or 1 there. The other old patch
  with a radius near 1.4 is not covered by this initial family.

## 5. Relation to radial localization and limits of the conclusion

This is not merely adding PSD radial features at one root. Its proof couples
the common-neighbor count about BOTH endpoints of an edge to a triangle sum.
For example, the natural root-local realization has a negative off-diagonal
entry for two short neighbors with a long base, while its origin cross terms
at that short root vanish. Its sum on those two neighbors is negative, so
that realization is not PSD and not even pointwise packing-positive at each
root. Positivity is restored only after summing across roots and charging the
long pair at its endpoints.

This establishes a distinct valid aggregate constraint, not a theorem that no
different root-PSD certificate could ever reproduce the same final density.
Whether it improves the optimized density is a numerical and verification
question. It specifically removes a previously mandatory pointwise far sign,
which the all-close common-neighbor experiments did not do.

The full geometric refinement would be an edge-anchored PSD block: for fixed
endpoints, a neighbor is described by its two endpoint distances (or axial
coordinate and transverse radius), with Gegenbauer degree in the normalized
transverse directions in dimension n-1. That requires four-point geometric
constraints. The elementary lens cap here is a much cheaper projected input.

## Diagnostics and scope

`check_lens_and_tail.py` reads the frozen direct 5D/6D/7D source coefficients,
checks the rational geometry, and records physical-tail and allowance values
in `lens_and_tail_results.json`. Tail extrema are float64 diagnostics only.
No integration files are edited, no optimization is run by this branch, and
no improvement to a verified density bound is asserted.

## Certified rational step allowances

For an immediate implementation retaining polynomial far expressions, the
following eta values have been certified against the frozen direct physical
functions. Select the unique interval containing c, with the intervals
half-open on the right except for the final one. Set eta=0 elsewhere. These
are alternatives to the saturated eta=-f/B, not additional allowances.

| Base interval | Cap | 5D eta | 6D eta | 7D eta |
|---|---:|---:|---:|---:|
| [46/25,48/25) | 2 | 1e-7 | 1e-7 | 5e-8 |
| [48/25,49/25) | 1 | 5e-7 | 1e-7 | 5e-8 |
| [49/25,2) | 1 | 3e-7 | 3e-8 | 1e-8 |
| [2,51/25] | 1 | 1e-7 | 5e-9 | 2e-10 |

Every decimal above denotes the exact rational. `certify_step_budgets.py`
reconstructs P in exact rational arithmetic from the source's finite-decimal
Laguerre coefficients, then proves f(c)+B*eta<0 on each CLOSED interval using
256-bit Arb second-order Taylor bounds and dyadic subdivision. Both endpoints
are checked for each neighboring formula, but a formula must be selected only
once at an endpoint; do not add two bin allowances there.

All twelve checks passed, at subdivision depths at most 2. See
`step_budget_verification.json` for exact eta fractions, source hashes and
replay results. This certifies the pair budget only, not any new kernel or
density bound. Far verification still must prove K-eta<=0 on each eligible
region and K<=0 on its complement, split at the short-radius and base-bin
boundaries. Those are polynomial inequalities, so the existing Bernstein
technology applies after adding the constant and making the new splits.
