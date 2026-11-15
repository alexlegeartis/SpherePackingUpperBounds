# A radial nonnegativity credit at the origin

## Statement and exact accounting

Let a valid packing-local energy have a positive semidefinite component

`(u+z(C))^T M (u+z(C))`, with `M>=0`, `M u=v`,
`kappa=v^T u`, and `z(C)=sum_{x in C} B(|x|)`.

The vector `B` is the actual piecewise radial feature vector. Other local
energy terms may be present and are left unchanged. Suppose an exact vector
`b` satisfies `M b=a`, and the radial function `w(r)=a^T B(r)` is nonnegative
throughout every closed radial piece used by the certificate. Define

`g=2*u^T*a-a^T*b`.

Then the following identity holds for every finite packing neighborhood,
including the empty neighborhood:

`E(C)-g = (u-b+z(C))^T M (u-b+z(C)) + 2*sum_{x in C} w(|x|) + unchanged_other_local_energies >= 0`.

The identity uses the exact range equation `M b=a`; no inverse of a singular
matrix is asserted. Set `u_new=u-b`, `v_new=v-a`, and `kappa_new=kappa-g`.
Then `M u_new=v_new` and `v_new^T u_new=kappa_new>=0`.

The new linear local energy has origin constant zero, origin coupling
`h_extra(x)=w(|x|)`, and zero kernel on two nonzero arguments. It exactly
cancels the change `-w` in the PSD origin coupling. Consequently the complete
physical `h`, the nonorigin kernel, the pair expression, and the cyclic
expression are **identically unchanged**. Only the total origin constant
falls by `g`. The density objective falls by `g/8` in the retained Gaussian
normalization. All old continuous sign proofs may therefore be reused.

This identity does not make a rejected source feasible. Its pair and cyclic
violations remain exactly the same. It gives a global bound only when the
source already has all required local, pair, cyclic and Fourier properties.

## Completed baseline application

Source: the retained completed `nested_repair2` certificate, objective
`0.6315414010365717...`, with its inherited all-branch verification.

Put `t=4*r-6` on the outer piece `[3/2,7/4]`, and use zero on the inner
pieces. Choose the exact nonnegative quartic

`w(r)=(1-t)*(t-253/500)^2*(a+b*t)`,

where

`a=1683024786574/10^20`, `b=125668197529249/10^20`.

Both coefficients are positive, so the factorization proves `w>=0` on the
entire closed outer interval. Both possible values at every radial jump are
nonnegative. An exact rank-13 reduction gives a vector `b` (the range-solution
vector, distinguished from the scalar factor coefficient above) satisfying
all 25 equations `M b=a_feature`. The exact completion identity then gives

`density_credit = 0.00000004992482939461474...`,

`new_objective = 0.6315413511117424... < 631541351112/10^12`.

This is a very small genuine improvement of the inherited completed bound.
No new cyclic or pair sign verification is necessary because those functions
are unchanged by identity. The data and fresh exact equations are in
`baseline_origin_credit/ORIGIN_CREDIT_VERIFIED.json`. Its generic status does
not assert a new global sign proof; that dependency is the retained baseline.

## Application to the first robust candidate

The candidate called `robust_alpha001_degree17` has rank 17, rather than the
rank 15 of its main functional basis. Two added Gram directions have radial
degree six. The exact reduction contains the ten degree-four mother-interval
features, five near-piece offsets, and those two extra directions. All 49
range equations are checked exactly.

The optimized piecewise polynomial, after exact rationalization and a
`10^-12` positive constant shift, has strict positivity on every closed
piece. Seven independent rational Sturm calculations show positive endpoints
and zero roots on each piece. The resulting exact density credit is

`0.00000010574488819783295...`.

Its formal objective becomes `0.630061051343187...`, below the rational
threshold `630061051344/10^12`. The source was not globally verified, so this
number remains conditional. The source's pair and cyclic expressions,
including any rejection witnesses, remain unchanged. Data are in
`robust_origin_credit/ORIGIN_CREDIT_VERIFIED.json`.

## Reproduction and reuse

For a source in the supported family, create an isolated exact range record:

`PYTHONPATH=<task>/deps python generic_origin_credit.py SOURCE OUTPUT_DIRECTORY`

Discover a radial polynomial numerically:

`PYTHONPATH=<task>/deps python optimize_generic_origin_credit.py OUTPUT_DIRECTORY`

The independent exact acceptance step uses Python's standard library:

`python certify_origin_credit.py OUTPUT_DIRECTORY`

The baseline factorized polynomial can be reconstructed with
`python certify_origin_credit.py baseline_origin_credit --baseline-factor`.

The exact step does not trust numerical polynomial minima. It rebuilds the
integer Gram matrix, checks the range and completion equations, proves each
radial sign by factorization or Sturm arithmetic, verifies a strictly
positive credit, and checks source file hashes again at completion.
