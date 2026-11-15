# Spherical neighbor-count theorem and exact verification

## Verified certificates

| Defining coefficient file | Spherical cutoff s | Exact repaired objective, decimal display | Proved count | Euclidean radius |
| --- | --- | --- | --- | --- |
| s0505_boundary_exact.json | 101/200 | 24.90465042056172 | 24 | 201/200 = 1.005 |
| s05075_taper12_exact.json | 203/400 | 24.972382344791278 | 24 | 50379/50000 = 1.00758 |
| s0515_count25_exact.json | 103/200 | 25.915792680106065 | 25 | 50767/50000 = 1.01534 |
| s0525_count26_exact.json | 21/40 | 26.941041293962012 | 26 | 102597/100000 = 1.02597 |

All comparisons with the next integer use exact rational arithmetic. The
last three radii are obtained by flooring the supported square-root radius
to the denominator-100000 grid, and their angular inequalities are checked
as rational inequalities. The second certificate permits every
radius R<=sqrt(200/197). Its verifier accepted 15257
triangle boxes and excluded 12081 boxes by a strictly negative upper bound on
the Gram determinant, with 54675 total nodes and maximum subdivision depth33.

The defining JSON coefficient files are the inputs to certify.py; the .npz
search files are numerical proposals and may continue to change during a
search. Files without a successful exact verification record are not proofs.

## Auxiliary bound on S^3

Let C be a finite set of unit vectors in R^4 with distinct-point inner products at most s<1. Let

    f(t)=sum_(j>=0) a_j U_j(t)/(j+1), a_j>=0, a_0>0.

For every anchor x in C use the same rotation-invariant kernel on its neighbors. Put a=x.y, b=x.z, t=y.z. Its nonorigin formula is

    K_x(y,z)=sum_l b(a)^T B_l b(b) H_l(a,b,t),
    H_l=[(1-a^2)(1-b^2)]^(l/2)
                P_l((t-ab)/sqrt((1-a^2)(1-b^2))),

where P_l is the Legendre polynomial normalized by P_l(1)=1. The expression is polynomial after expanding, including at a=+-1 or b=+-1. Each B_l is PSD. The degree-zero augmented origin block [[kappa,v^T],[v,B_0]] is PSD, with K_x(y,x)=v^Tb(a) and K_x(x,x)=kappa. Restrict the nonorigin features to a<=s, and set them to zero at other nonorigin points. The origin has its separately specified feature.

For a fixed anchor, the transverse vectors y-(x.y)x lie in the three-dimensional space perpendicular to x. The addition theorem for degree-l spherical harmonics on its S^2 gives positivity of H_l. Multiplication by radial features and the PSD block B_l preserves positivity. The origin factorization completes the PSD kernel.

Assume the pair polynomial

    f(a)+2v^Tb(a)+sum_l (1-a^2)^l b(a)^T B_l b(a)

is nonpositive on [-1,s]. Assume also

    K(a,b,t)+K(a,t,b)+K(b,t,a)<=0

whenever a,b,t are in [-1,s] and 1+2abt-a^2-b^2-t^2>=0. These are exactly the feasible Gram matrices of three distinct unit vectors with the prescribed angle constraint.

For N=|C|, PSD summed over anchors and grouped by coincident indices yields

    sum_(x,y in C) f(x.y) <= N[f(1)+kappa].

The positive harmonic expansion of f gives the reverse bound N^2 a_0. Therefore

    N <= [f(1)+kappa]/a_0.

An exact objective strictly below 25 proves N<=24.

## Cardinality-specific positivity

For the purpose of excluding a code of exactly N=25 points, the kernel energy
only has to be nonnegative on sets with exactly m=24 nonorigin neighbors. Let
the radial basis have first entry b_0(a)=1 (as the Chebyshev basis does), and
write t=sum_(y != x) b(x.y). Then t_0=m exactly. The degree-zero anchored energy is

    kappa+2 v^T t+t^T B_0 t.

It is sufficient that the following reduced matrix be PSD:

    R = [[kappa+2m v_0+m^2(B_0)_(0,0), (v_rest+m(B_0)_(0,rest))^T],
         [v_rest+m(B_0)_(0,rest),         (B_0)_(rest,rest)        ]].

The full origin block and B_0 can be indefinite. With t_0=m, their quadratic
energy is exactly (1,t_rest)^T R(1,t_rest)>=0. Higher angular blocks remain PSD.
This is an exact restriction to a cardinality equality, not an approximation
or a PSD claim for the unrestricted kernel.

To check the summation argument explicitly, let E_x be the full anchored
energy and put P(a)=f(a)+2v^Tb(a)+K(a,a,1). The sum over anchors is

    sum_x E_x = N kappa
                  +sum_(x != y)[P(x.y)-f(x.y)]
                  +2 sum_(unordered distinct triples) Gamma(x.y,x.z,y.z).

The pair and cyclic signs give

    0 <= sum_x E_x <= N[kappa+f(1)]-sum_(x,y)f(x.y)
                    <= N[kappa+f(1)-N a_0].

Thus a_0>0 and [kappa+f(1)]/a_0<25 exclude a 25-point code. A larger code
contains a 25-point subset, so it too is excluded. The uniform repair below
adds an actually PSD kernel and therefore preserves this restricted
nonnegativity argument without requiring full PSD of the original kernel.

For this certificate mode, integer Gram rays specify R and all higher blocks.
The scalar entries v_0, (B_0)_(0,0), and (B_0)_(0,rest) are stored as unrestricted
rational numbers. The verifier reconstructs the remaining origin-block entries
from R and checks exactly that the restriction map L, defined by
L(1,t_rest)=(1,m,t_rest), satisfies L^T A L=R.

## Uniform repair of a small cyclic error

Let q(t)=sum_j q_j U_j(t)/(j+1) have q_j>=0, q_0>0, and q(t)<=0 on [-1,s]. On nonorigin neighbors, the kernel

    T(y,z)=q(y.z)/q_0-1

is PSD because its constant harmonic coefficient is zero and all others are nonnegative. Give it zero origin row. On distinct code points it is <=-1 and its diagonal is H=q(1)/q_0-1.

If a candidate has cyclic excess at most 3 epsilon and pair excess at most eta, add epsilon*T and replace the scalar constant a_0=1 by

    a_0_new=1-epsilon*H-eta.

For epsilon>=0, eta>=0, a_0_new>0, the resulting pair and cyclic inequalities are valid. If the original objective is B=f(1)+kappa, the repaired objective is exactly

    (B-epsilon*H-eta)/(1-epsilon*H-eta).

The repair does not need to be re-expanded in the finite stabilizer basis; its PSD and global signs have this separate proof.

## Free-scalar spherical local energy and its four-point transfer

The harmonic coefficient signs of the scalar function are unnecessary for
the first inequality in the auxiliary-bound proof. In particular, let

    g(t)=sum_j c_j U_j(t)/(j+1)

have arbitrary real coefficients. Suppose an anchored kernel has the full
PSD origin block and PSD higher harmonic blocks, and satisfies the pair
and cyclic inequalities above with g in place of f. For every finite code
C with distinct-point inner products at most s, including the empty code,
write N=|C| and B=g(1)+kappa. Expanding the sum of its anchored energies gives

    0 <= sum_x E_x
       = N kappa + sum_(x != y)[P(x.y)-g(x.y)]
                   +2 sum_(unordered distinct triples) Gamma
       <= B N - sum_(x,y in C) g(x.y).

Consequently the right-hand side is a nonnegative quadratic energy on every
such code. This conclusion requires neither a positive constant coefficient
nor nonnegative higher scalar coefficients. Certificates that assert PSD
only after fixing N cannot be used for this unrestricted energy conclusion.

Let R satisfy the Euclidean angular condition from the next section, and
put chi(x)=1 for 1<=|x|<=R and zero otherwise. For every separation-one
packing neighborhood in R^4, the following augmented kernel has a
nonnegative total energy:

    K(0,0)=0,
    K(0,x)=K(x,0)=B chi(x)/2,
    K(x,y)=-chi(x)chi(y) g((x.y)/(|x||y|)).

Indeed its total energy on an anchor and its packing neighbors is exactly
B N-sum_(i,j)g(u_i.u_j), where u_i are the selected normalized directions.
This is positivity on packing neighborhoods; no unrestricted PSD claim is
made for the displayed augmented kernel. The resulting global pair cost is

    2 K(0,x)+K(x,x)=kappa chi(x),

and its global cyclic contribution is the negative sum of g at the active
triangle vertices, each weighted by the indicators of its two incident
edge lengths. Adding this energy to a globally PSD kernel is valid in the
packing summation argument whenever the resulting global pair and cyclic
conditions are verified. The spherical anchored positivity depends on a
packing center and three of its neighbors, so it supplies four-point
information to that global argument.

Uniform repair also works for this free-scalar cone. If the raw spherical
cyclic excess is at most 3 epsilon and the raw pair excess is at most eta,
replace g by g-delta with delta=epsilon H+eta and add epsilon T to its
anchored kernel. Then B is replaced by B-delta and kappa is unchanged.
There is no requirement delta<1 in this energy setting. Its exported
global kernel must use the repaired g, and the resulting global signs
must be checked after that substitution. The exact verifier's
verification_goal="local_energy" branch implements this contract.

The compact numerical model is local_energy_api.py. The file
accepted_local_energy_rays.json exports exact repaired scalar coefficients,
kappa, B, support radius, and source-certificate hashes. Its accepted
count-optimized rays have very small nonconstant scalar coefficients;
allowing free scalar coefficients in a joint global optimization is needed
to explore the additional angular information efficiently.

## From a spherical bound to a Euclidean shell count

Suppose points around a separation-one packing center have radii A,B in [1,R]. Their normalized directions have inner product at most (A^2+B^2-1)/(2AB). For 1<=R<(1+sqrt(5))/2 the maximum on the radius rectangle is 1-1/(2R^2): fixing either radius gives no interior maximum, and the corner comparison follows from (R-1)(R^2-R-1)<=0.

Consequently a spherical count of 24 at cutoff s proves at most 24 packing neighbors in [1,R] whenever 1-1/(2R^2)<=s. In particular s=101/200 would support R=201/200, since its required cosine is 20401/40401<101/200.

## Exact certificate contract

Numerical SDP results are proposals only. Rationalization expresses every matrix as a sum of integer outer products with a common positive squared denominator, giving exact PSD. Scalar harmonic coefficients are nonnegative rationals for a count certificate; they may be signed rationals in the explicitly marked local-energy mode. The verifier expands the pair and cyclic polynomials with integer coefficients after clearing positive denominators.

All scalar and pair signs are certified on their entire closed intervals with outward integer Bernstein bounds. For triangles the verifier subdivides the full cube [-1,s]^3. A box is accepted when its cyclic upper coefficients are nonpositive, or excluded only when the upper coefficients of its Gram determinant are strictly negative. Remaining boxes are subdivided. Explicit overflow guards protect integer averaging. UNKNOWN and VIOLATION fail verification.

The first sampled s=.505 candidate was rejected by this exact process at the rational admissible triple (-12499/12800,-377/1280,523/6400). Direct evaluation of the Gram rays agreed with the polynomial expansion at that witness. It was added to the numerical exchange search. The rejection is retained as evidence that sampling alone is not the certificate.
