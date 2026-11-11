# Five-dimensional density below 0.509996147

## Exact frozen parameters

Minimum separation is 1 and kernel support radius is R=3/2. Dimension n=5,
Gaussian–Laguerre radial degree 31, Chebyshev radius degree p=4, and angular
Gegenbauer degrees 0 through K=12. The source is `source_exact.json`, SHA-256
235920216afbd10b38be5efd3672e9d94b5ec5af8fc54aca3a63d77d9997721f.

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

## Elementary lens caps: all six actual ranges

Set r=27/25 and J=[41/25,54/25]. For endpoints x,y at distance d and common
neighbors z_i with |z_i-x|,|z_i-y|<=r, their midpoint m obeys

    |z_i-m|^2 <= r^2-d^2/4 = rho^2.

Because all center separations are at least 1, the centroid identity gives

    M(M-1)/2 <= sum_(i<j)|z_i-z_j|^2
             = M sum_i |z_i-zbar|^2 <= M^2 rho^2.

Thus M<=1/(1-2rho^2) whenever rho^2<1/2. At each lower endpoint the
rational real count upper is strictly below B+1, yielding by integrality:

| Base LENGTH interval | B |
|---|---:|
| [41/25,17/10) = [1.64,1.70) | 83 |
| [17/10,7/4) = [1.70,1.75) | 8 |
| [7/4,9/5) = [1.75,1.80) | 5 |
| [9/5,46/25) = [1.80,1.84) | 3 |
| [46/25,48/25) = [1.84,1.92) | 2 |
| [48/25,54/25] = [1.92,2.16] | 1 |

These are length thresholds, not squared radii. Both short arms are at most
27/25. Exact cap arithmetic is replayed at the six thresholds and at every
step-bin lower edge. No transverse projection, Rankin argument, or additional
spherical-code bound is used.

The frozen budget has 104 bins of width 1/200 across J. Each bin has an exact
nonnegative rational eta and the corresponding cap B. On each CLOSED bin the
constructor's `--verify` mode reconstructs the exact physical polynomial and
proves f(d)+B eta<0 using 256-bit Arb second-order Taylor enclosures and
subdivision. It ignores stored construction minima and success flags.
The canonical physical payload is bound to the source coefficients, basis,
variable, exponential factor, degree and absence of additional physical terms.

Use each eta only on its open bin interior. At every internal endpoint use
the minimum of its two adjacent eta values; at either support endpoint and
outside J use zero. At a cap transition the smaller adjacent allowance is
valid with the new cap, since that cap is no larger than the preceding cap.
Only one allowance applies at any base length. The verifier checks all
boundary values and the exact allowance-query conventions.

## Weighted incidence and tail funding

Every eligible triangle has a unique base d in J, because d>=41/25>r and
its other two edges have lengths at most r. If M_e is the number of short
common neighbors of a base pair e, then

    sum_eligible_triangles eta(d) = sum_base_pairs eta(d) M_e
                                <= sum_base_pairs B(d) eta(d).

Permit C_K<=eta(d) on eligible triangles and retain C_K<=0 elsewhere.
In the one-long-edge region this is K(a,b,c)<=eta(c). All compact cyclic
conditions remain unchanged, since J starts strictly above R.

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

A recorded exact example at short legs 1,1 and base 1883/1000 has K>0, about
4.303908118922817e-7. This makes the old pointwise far checker invalid for this
source; it cannot be substituted for `certify_wedge_far.py`.

## Schur saturation, density and unrestricted packings

Write Q_0=[[a,b^T],[b,B]] with B positive definite. Exact rational arithmetic
verifies alpha=a-b^T B^-1 b=10^-12. Replacing a by its Schur value yields a
PSD block with factorization [[b^T B^-1],[I]] B [[B^-1 b,I]]. All nonzero
pair, cyclic and far evaluations remain unchanged. Only the origin cost
falls by alpha. With H(0)>0, the density upper bound is

    (pi^2/60) * [P(0)+b^T B^-1 b+sum_i lambda_i B_i]/H(0).

Every numerator and denominator term is reconstructed exactly and compared
with `prerequisites.json`. Arb enclosure and an independent rational upper
bound for pi give density strictly below **0.509996147**.

For periodic packings, enlarge the period lattice beyond the finite
kernel/lens interaction diameter. Count actual physical centers and
pair/triangle orbits under simultaneous translation, including translated
copies of motif centers. The finite incidence and centroid inequalities
continue to apply. Schwartz Fourier sums converge, kernel/lens sums are
locally finite, and the PSD/Poisson argument gives the displayed inequality.
The standard large-patch periodic approximation then extends it to all
sphere packings, not only lattice or periodic packings.

## Verification and trust boundary

`verify.py` performs a fresh complete producer replay. The same-source pair
and compact cyclic checkers are unchanged. The new far checker reconstructs
rational polynomial coefficients, proves signs by Bernstein convexity, and
uses directed binary64 rounding at every subdivision operation. Exact
geometry determines only valid exclusions. The 104 tail-budget inequalities
use rigorous Arb intervals, not optimizer tolerances. The included general
budget file also contains frozen 6D/7D records and replays them, but only its
5D record enters this certificate.

A separate completed audit independently reconstructed radial algebra,
PSD/Schur arithmetic, count endpoints, all tail budgets and the new far
polynomials; its compact cyclic stage was source-reviewed and freshly
replayed rather than reimplemented. See `INDEPENDENT_AUDIT.md`. Its algebra
checker is included in each producer run as `audit_algebra_exact.py`.

Trusted components include the addition theorem, Poisson summation and
periodic approximation, published near-neighbor count results, Python integer
and rational arithmetic, FLINT/Arb, NumPy and correctly rounded IEEE binary64
arithmetic. This is computer-assisted verification, not a proof-assistant
formalization. The numerical optimization and stored numerical success flags
are not premises of the proof.
