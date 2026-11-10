# Exact origin-coordinate tightening

This dominates the separately explored integer cut on the frozen 5D witness.
It is ordinary PSD slack removal, not a new multipoint theorem or a claim to
improve dimensions outside the existing registered set.

## Matrix lemma

Partition the exact degree-zero block as Q=[[a,b^T],[b,B]]. The source blocks
are positive definite, so B is positive definite and
alpha=a-b^T B^-1 b=1/(Q^-1)[0,0]>0. Set W=[(B^-1 b)^T;I]. Then

    Q-alpha*e0*e0^T = W*B*W^T >= 0.

This is an exact rational identity, verified by the producer research script and verify_algebra.py along with
exact LDL positivity of Q and B. The maximal subtraction gives a singular
PSD matrix. The mathematical theorem only requires semidefiniteness, not
strict definiteness. An unchanged source verifier insisting on positive
pivots will reject this valid singular block; use the factorization above
instead, or subtract slightly less than alpha to retain positive definiteness.

## Why all continuum constraints are unchanged

The origin feature coordinate equals one exactly at displacement zero and
zero at every nonzero displacement. Consequently subtracting alpha*e0*e0^T
changes K(x,y) only at (x,y)=(0,0). It does not change K(x,0), K(x,x), or
K(x,y) for nonzero x,y, and therefore changes neither the pair inequality
nor any cyclic, far, angular, annular, Fourier or tail condition. Closed
support boundaries are unaffected. No origin limit or continuity is needed.

For a center in a periodic packing, the finite displacement sum includes
zero once. Other period translates of that center have nonzero displacement
and are included as ordinary neighbors when inside support. No repeated
origin occurs. Each local PSD sum therefore improves by alpha, and summing
over N centers improves the total by N*alpha, not 2N*alpha. The factor two
on off-diagonal kernel entries and the cyclic grouping are unchanged.

Equivalently, simply insert the modified PSD kernel into the source theorem:
the objective's sole changed term is K(0,0), reduced by alpha. If the stored
matrix is unnormalized and the theorem's kernel is divided by H0, the
objective decreases by alpha/H0, not alpha or alpha/H0^2.

## Numerical normalization checked

In dimensions 4 and 5 the matrices are raw Q, and H0 is reconstructed from
the first 32 Laguerre coefficients as sum_k (-1)^k*c_k. Minimum separation
is 1; density multipliers are pi^2/32 and pi^2/60 respectively. The published
neighbor-count assumptions and their nonnegative weights remain unchanged.

In dimensions 19, 20, 21, the exact certificate's normalization is H0. Its raw
degree-zero matrix is divided by H0 in the theorem. Minimum separation
squared is certificate.exclusion_squared, namely 19/8,5/2,21/8. The density
multiplier is vol(B_1^n)*(exclusion_squared/4)^(n/2). The code checks these
normalizations and computes both approximate decrements and rational lower
bounds, using pi>223/71 and an exactly verified rational square-root lower
bound in odd dimensions.

## Results and scope

Exact matrix tightening passes for dimensions 4, 5, 19, 20, 21. Density decrements:

- 4: approximately 1.59156666237462e-6.
- 5: approximately 2.14303381036118e-5.
- 19: approximately 7.79535801673426e-9.
- 20: approximately 4.63247944494798e-9.
- 21: approximately 3.83182287152171e-9.

The exact fractions and conservative outward new decimals are in
schur_records.json. Those decimals subtract a rigorous rational lower
bound for the decrement from the previously certified outward upper, so
they sacrifice some digits. The independent Arb objective reconstruction sharpens
them to the strict bounds in README.md. The complete source certificates were replayed successfully in an isolated snapshot. The package verify.py reruns all eight source verification paths and then verifies the exact Schur factorization, independent objective reconstruction and strict outward rounding. External neighbor-count mathematical inputs in dimensions 4 and 5 remain inherited. This is a computational certificate audit, not formalization or a novelty claim.

## Why this simple subtraction does not extend to 48D

The 48D source uses a rank-one degree-zero feature g0, with g0(0)=-B and
g0(x)=h(x) on the admissible radial region, scaled by M/B. Its degree-zero
block is already Schur-saturated. The corresponding two-coordinate block
has entries MB,-Mh,Mh^2/B; subtracting any positive alpha from its origin
entry makes its determinant -alpha*M*h^2/B<0. Thus degree-zero blockwise
origin subtraction has no positive allowance in that representation. This
does not rule out redesigning the 48D kernel or exploiting genuinely stronger
constraints; it only closes this particular postprocessing route.
