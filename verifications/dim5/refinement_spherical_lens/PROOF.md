# Five-dimensional density below 0.509953498

## Exact frozen parameters

Minimum separation is 1 and kernel support radius is R=3/2. Dimension n=5,
Gaussian–Laguerre radial degree 31, Chebyshev radius degree p=4, and angular
Gegenbauer degrees 0 through K=12. The source is `source_exact.json`, SHA-256
28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949.

Count radii are exactly 1.003012, 1.008096, 1.013258, 1.017445, with counts
45, 46, 47, 48 and the nonnegative rational weights in `witness.json`.
These thresholds are each 10^-6 below the published thresholds in
[de Laat, de Oliveira Filho and Vallentin, Table 2](https://arxiv.org/html/1206.2608),
so monotonicity suffices. The allowance is sum_i lambda_i 1_(d<r_i).
At each threshold the smaller allowance is checked explicitly. The published
spherical SDP inputs themselves are external mathematical assumptions.

The physical function is exactly the unchanged direct-5D function. For
X=pi |x|^2, it is f(x)=P(X) exp(-X), where

    P(X) = sum_(k=0)^31 c_k L_k^(3/2)(2X)/L_k^(3/2)(0).

Under Fourier convention exp(-2 pi i x dot xi), its Fourier polynomial H
multiplies the kth coefficient by (-1)^k. A second reconstruction verifies
this identity directly in monomial coordinates. Global exact rational Sturm
checks establish H>0 on [0,infinity), and P<0 for X>=9 pi/4 using a rational
lower enclosure of pi. Normalization is exactly
19999999999999998358193/20000000000000000000000. Both Gaussian correction
coefficients are zero and both correction-island lists are empty.

## PSD kernel and unchanged compact hypotheses

For radii in [1,R], let w contain T_0,...,T_4 evaluated at
(2r-R-1)/(R-1); set these features to zero outside [1,R]. Add one isolated
origin coordinate in harmonic degree zero. The rational block Q_0 has size 6
and each Q_l for l>0 has size 5. The spherical harmonic addition theorem and
positive semidefiniteness of these blocks give a globally PSD kernel K.
All 66 positive pivots of the unsaturated matrices are proved rationally.

For nonzero vectors at distance r from the root, the pair condition is

    f(r) + 2K(x,0) + K(x,x) <= sum_i lambda_i 1_(r<r_i).

It is proved on the complete closed support interval [1,R], split at every
count radius, with explicit smaller-allowance endpoint tests. Outside support,
its kernel and count terms vanish and the global physical-tail sign suffices,
except for the added lens charge described below.

For a triangle, let C_K be the sum of the three kernels rooted at its
vertices. The complete compact region [1,R]^3 has C_K<=0. Every such triple
satisfies triangle inequalities because R<2. The exact polynomial numerator
(abc)^K C_K is checked by tensor Bernstein initialization and outward-rounded
dyadic subdivision. Positive denominators preserve its sign. With two edges
outside support all root kernels vanish. With exactly one such edge only the
opposite root kernel can survive; that far condition is the modified one.

## Conservative transverse projection and exact code caps

Set r=27/25 and J=[3/2,54/25]. Write a common neighbor of endpoints at
base length d as (ell,z) relative to their midpoint and axial direction.
Both incident lengths are in [1,r], so |ell| <= (r^2-1)/(2d). Put
A=(r^2-1)^2, D=A/d^2, U=r^2-d^2/4 and L=1-d^2/4-D/4. The transverse
vectors z lie in R^4, have squared norms in [L,U], and transverse squared
separation at least 1-D. For d in [3/2,41/25], L>0, and

    normalized dot <= (u+v+D-1)/(2 sqrt(uv)) <= 1-(1-D)/(2U).

The second inequality follows by differentiating in each squared norm:
L-U+1-D=2-r^2-5A/(4d^2)>0. The cosine upper s(d) decreases with d,
since with x=d^2 the derivative numerator of the subtracted quotient is
(x-A)^2+A(4r^2-A)>0. Exact norm and denominator margins, L at 41/25,
and the monotonicity remainder are reconstructed in spherical_lens_caps.py.
Positive transverse norms exclude zero vectors, and s<1 makes normalization
injective. This is the conservative bound 1-(1-D)/(2U), not the subsequently
proposed sharper projection lemma.

At each of d=3/2,31/20,8/5 the embedded rational dimension-4 polynomial
Q(t)=sum_l a_l P_l^(4)(t) has a_0=1 and all a_l>=0. Exact rational Sturm
root counts and negative endpoint signs prove Q<0 throughout [-1,s(d)].
The addition theorem gives M^2<=sum_(i,j)Q(t_ij)<=M Q(1), hence
M<=floor(Q(1)). Independently replayed coefficients, monomial polynomials,
objectives and integer floors give caps 10,9,8. No external spherical-code
optimum or stored success boolean is accepted for these new caps.

For d>=41/25 the unnormalized transverse dot upper U-(1-D)/2 is strictly
negative (exactly -3583/4202500 at 41/25) and decreasing. Strictly obtuse
vectors in R^4 have at most 5 elements: a nonzero affine dependence among
6 vectors, split by coefficient sign, expresses a squared norm as a strictly
negative cross sum. A zero vector can coexist with no second vector under
this inequality. Thus M<=5 independently of rank or normalization.

The midpoint centroid identity separately gives M<=1/(1-2rho^2) where
rho^2=r^2-d^2/4<1/2. Taking integer floors and the better transverse cap:

| Base LENGTH interval | B | Source |
|---|---:|---|
| [1.50,1.55) | 10 | Dimension-4 exact spherical polynomial |
| [1.55,1.60) | 9 | Dimension-4 exact spherical polynomial |
| [1.60,1.64) | 8 | Dimension-4 exact spherical polynomial |
| [1.64,1.70) | 5 | Strictly obtuse transverse vectors |
| [1.70,1.75) | 5 | Strictly obtuse transverse vectors |
| [1.75,1.80) | 5 | Minimum of transverse and centroid caps |
| [1.80,1.84) | 3 | Centroid |
| [1.84,1.92) | 2 | Centroid |
| [1.92,2.16] | 1 | Centroid |

These are base lengths, not squared radii. Both short arms are at most
27/25. construct_step_allowances.py validate_cap_geometry reconstructs every
embedded projection and polynomial, then replay proves the analytic budgets.
See RANKIN_AND_SPHERICAL_PROOF.md for the detailed geometric argument.

The frozen budget has 132 bins of width 1/200 across J. On each CLOSED bin,
the exact physical polynomial and corresponding rational B,eta prove
f(d)+B eta<0 by 256-bit Arb Taylor bounds and subdivision. Construction
minima and stored success flags are ignored. The canonical physical payload
is bound to source coefficients, dimension, basis, variable, exponential
factor, degree and absence of additional physical terms.

Use eta only on its open bin interior. At each internal endpoint use the
minimum of adjacent eta values; at both support endpoints and outside J use
zero. In particular eta(3/2)=0. At a cap transition the new cap is no larger
than the old one, so the smaller allowance remains valid. Only one allowance
applies at each length. All partition edges and boundary values are checked.

## Weighted incidence and tail funding

Every eligible triangle has a unique base d in J, because d>3/2>r and
its other two edges have lengths at most r. If M_e is the number of short
common neighbors of a base pair e, then

    sum_eligible_triangles eta(d) = sum_base_pairs eta(d) M_e
                                <= sum_base_pairs B(d) eta(d).

Permit C_K<=eta(d) on eligible triangles and retain C_K<=0 elsewhere.
In the one-long-edge region this is K(a,b,c)<=eta(c). All compact cyclic
conditions remain unchanged, because the allowance at d=R=3/2 is zero and
positive allowances apply only for d>R.

The root-PSD and pair/Fourier expansions, per period cell, have the identity

    S_f + S_K = N[f(0)+K(0,0)]
              + 2 sum_pairs [f(d)+2K(x-y,0)+K(x-y,x-y)]
              + 2 sum_triangles C_K.

Both unordered pairs and unordered triangles have multiplicity two. The
incidence bound therefore charges exactly B(d)eta(d) to the long-pair
expression, with no missing factor two or three. At d in J the kernel and
near-neighbor pair terms vanish, while f(d)+B(d)eta(d)<=0. Elsewhere the old
pair inequality and physical-tail sign apply. No new objective term arises.

The new far checker proves this complete condition on three symmetric
rectangles formed by splitting both radius coordinates at 27/25. Exact
squared-distance extrema determine every bin intersecting a closed box.
It uses the minimum eta over all those bins, and zero when the box touches
outside-support or outside-short-radius regions. Boxes on the outer side of
a short-radius boundary are initialized closed and retain K<=0 at that
boundary. Thus the stronger outside-side limiting inequality is not lost.
The whole kernel support boundary is also covered by the closed pair and
triangle/far proofs; immediately outside it the kernel vanishes.

The original pointwise far sign is not assumed. The positive allowance is
used in 113 verified far leaves. The required checker is
certify_general_lens_far.py, which replays all embedded cap and budget proofs.

## Schur saturation, density and unrestricted packings

Write Q_0=[[a,b^T],[b,B]] with B positive definite. Exact rational arithmetic
verifies alpha=a-b^T B^-1 b=10^-12. Replacing a by its Schur value yields a
PSD block with factorization [[b^T B^-1],[I]] B [[B^-1 b,I]]. All nonzero
pair, cyclic and far evaluations remain unchanged. Only the origin cost
falls by alpha. With H(0)>0, the density upper bound is

    (pi^2/60) * [P(0)+b^T B^-1 b+sum_i lambda_i B_i]/H(0).

Every numerator and denominator term is reconstructed exactly and compared
with `prerequisites.json`. Arb enclosure and an independent rational upper
bound for pi give density strictly below **0.509953498**.

For periodic packings, enlarge the period lattice beyond the finite
kernel/lens interaction diameter. Count actual physical centers and
pair/triangle orbits under simultaneous translation, including translated
copies of motif centers. The finite incidence and geometric/code inequalities
continue to apply. Schwartz Fourier sums converge, kernel/lens sums are
locally finite, and the PSD/Poisson argument gives the displayed inequality.
The standard large-patch periodic approximation then extends it to all
sphere packings, not only lattice or periodic packings.

## Verification and trust boundary

`verify.py` performs a fresh complete producer replay. The same-source pair
and compact cyclic checkers are unchanged. The new far checker reconstructs
rational polynomial coefficients, proves signs by Bernstein convexity, and
uses directed binary64 rounding at every subdivision operation. Exact
geometry determines only valid exclusions. The 132 five-dimensional tail-budget inequalities
use rigorous Arb intervals, not optimizer tolerances. The included general
budget file also contains 132 frozen 6D and 132 frozen 7D records. All 396
are replayed, but only the 5D record enters this density certificate.

A separate completed audit independently reconstructed radial algebra,
PSD/Schur arithmetic, count endpoints, all 132 five-dimensional budgets,
projection geometry, the three dimension-4 code polynomials and new far
polynomials; its compact cyclic stage was source-reviewed and freshly
replayed rather than reimplemented. See `INDEPENDENT_AUDIT.md`. Its algebra
checker is included in each producer run as `audit_algebra_exact.py`.

Trusted components include the addition theorem, Poisson summation and
periodic approximation, published near-neighbor count results, Python integer
and rational arithmetic, FLINT/Arb, NumPy and correctly rounded IEEE binary64
arithmetic. This is computer-assisted verification, not a proof-assistant
formalization. The numerical optimization and stored numerical success flags
are not premises of the proof.
