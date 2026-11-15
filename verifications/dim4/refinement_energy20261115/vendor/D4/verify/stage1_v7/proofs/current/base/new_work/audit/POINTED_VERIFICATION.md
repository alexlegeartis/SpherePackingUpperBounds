# Independent audit of pointed ball-count repairs

This note audits `new_work/pointed/interval_certificate.py`,
`new_work/pointed/build_repair.py`, and their integration into
`new_work/counts/verify_augmented_certificate.py`. It records a mathematical
justification of the new local energy and its exact interval proof. A
particular packing bound is accepted only after the complete global
verifier reports that all conditions passed.

## Ball count and the pointed energy

Let the packing have minimum center separation one. A closed ball of
squared radius `a2 < M/[2(M+1)]` contains at most M centers. Indeed, if it
contained n=M+1 centers x_i, then

    n(n-1)/2 <= sum_(i<j)|x_i-x_j|^2
               = n sum_i |x_i-mean(x)|^2 <= n^2 a2,

contradicting the strict radius condition. This is the exact condition
checked by `validate_caps`.

Fix a packing center at the origin and let S be its selected nonorigin
neighbors with radii in the compact global kernel support [1,T]. For each
y in S, choose v uniformly on the latitude `v.yhat=sigma` of S^3, where
`-1<=sigma<=0`; sigma=-1 is the point mass at -yhat. Let P(x|y) be the
probability that x belongs to the ball centered at d v with squared radius
a2. For a nonnegative radial weight w, the quantity

    M sum_y w(y) - sum_(x,y in S) w(y) P(x|y)

is nonnegative: for every y and every v the selected neighbors in the ball
are a subset of a separation-one configuration, so their count is at most
M. Omitting packing points outside [1,T] preserves this bound.

Write chi(x) for the active-annulus indicator. The corresponding symmetric
augmented kernel has

    Kp(0,0)=0,
    Kp(0,x)=M w(x)/2,
    Kp(x,y)=-chi(x)chi(y)[w(y)P(x|y)+w(x)P(y|x)]/2.

Thus it is explicitly zero whenever either nonorigin argument is outside
the active annulus. Merely giving w compact radial support would not be
sufficient: the other term in the symmetric formula could remain nonzero.

Its total energy is the displayed nonnegative quantity. This is positivity
on packing neighborhoods; unrestricted positive definiteness is neither
claimed nor needed in the packing summation argument. Nonnegative linear
combinations of these energies remain valid.

For r>0, write `c(r)=(r^2+d^2-a2)/(2dr)`. The assumptions `d>0` and
`d^2>a2` imply c(r)>0. Since sigma<=0, P(x|x)=0. Thus Kp(x,x)=0 and the
global pair contribution is exactly M w(x). Every nonorigin kernel value
Kp(x,y) is nonpositive.

## Conditional probability on S^3

Put t=xhat.yhat. For -1<sigma<=0, write

    v=sigma yhat + sqrt(1-sigma^2) u,

where u is uniform on the S^2 in the orthogonal three-dimensional space.
The coordinate of u along any fixed unit direction is uniform on [-1,1].
Consequently, with A=c(r)-sigma t and
`D=sqrt((1-sigma^2)(1-t^2))`, the probability is

    clip_[0,1]((1-A/D)/2)

when D>0. When D=0 it is one if A<=0, and zero if A>0. The weak inequality
at A=0 is required because the ball is closed. The sigma=-1 branch reduces
directly to the indicator of `c(r)+t<=0`.

For a radius interval [rlo,rhi], c is convex because d^2-a2>0, so its maximum
is attained at an endpoint. Its stated lower bound

    rlo/(2d) + (d^2-a2)/(2d rhi)

is valid term by term. Since -sigma>=0, interval bounds for A use tlo in
the lower bound and thi in the upper bound. Bounds for 1-t^2 use its two
endpoint minima and its maximum at zero when the interval crosses zero.
The integer square-root routine encloses D from below and above.

For A>=0, the probability decreases as A grows and increases as D grows.
For A<=0, it decreases as A grows and decreases as D grows. The four
sign-dependent choices in `probability_bounds` follow precisely these
monotonicities. If the A interval straddles zero, the lower and upper
bounds use the corresponding separate cases. Its explicit D=0 branches
give conservative enclosures, including the boundary probability one.

`audit_pointed.py` independently compares these interval results against
exact radical inequalities at 16,282 rational points in 1,800 rectangles,
including sigma=-1, sigma=0, t=+-1, and the degenerate equality case. It
uses rational squares and comparisons, with no floating-point or
trigonometric oracle. These finite checks supplement the continuum
monotonicity proof above.

## Exterior interval comparison

For an exterior branch the third side lies outside [1,T]. The core's
coordinates are radii a,b and z in [0,1], with

    t=-1+z h(a,b),  h=((a+b)^2-T^2)/(2ab).

On 1<=a,b<=T<=2, h is nonnegative and increases in each radius, since
its derivative in a has numerator a^2-b^2+T^2>0; the b derivative is
analogous. Thus the corner enclosure for t is valid and remains within
[-1,1/2]. The radial weight uses nonnegative Bernstein coefficients on a
basis interval containing the entire branch. Each Bernstein basis
function is unimodal, with its maximum at j/degree, so the stated endpoint
minimum and possible interior maximum enclose it.

The pointed correction is

    C = sum_caps beta/2 [w(b)P(a|b)+w(a)P(b|a)] >= 0.

The routine bounds it using products of nonnegative interval endpoints.
The polynomial core represents `(2ab)^L (Gamma_base-eta)`. If U is an
outward Bernstein upper bound for that polynomial and C_lower bounds C,
then

    U <= (2 a_lower b_lower)^L C_lower

implies the desired `Gamma_base-eta-C<=0` throughout the box. Both factors
on the right are nonnegative and bounded from below. Conversely, an exact
lower polynomial bound at a corner greater than the corresponding factor
times C_upper certifies a violation at that corner.

The initial Bernstein coefficients are rounded outward into signed
integers. The chosen precision bounds their absolute values below 2^59,
so the core's averaging splits do not overflow int64. The lower and upper
arrays are split with outward rounding. UNKNOWN, a depth limit, a node
limit, or a positive exact witness does not count as acceptance.

## Origin shift and the 11-to-15 embedding

Let M0 be the degree-zero nonorigin Gram matrix and let v_old be the origin
coupling whose global pair inequality is to be retained. For pointed rays
with weights beta_j, set

    h = sum_j beta_j M_j w_j/2,
    v_new = v_old-h.

If `M0 u_new=v_new`, choose `kappa_new=v_new.u_new`. The augmented origin
block is PSD: its quadratic energy is

    (u_new + sum_x b(x))^T M0 (u_new + sum_x b(x)).

This argument also holds for a singular M0, provided v_new lies in its
range. The builder explicitly works in the 11 functional coordinates of
the hybrid model, where it solves the rational linear system. It then
embeds into the 15 piecewise coordinates. The full verifier checks
`M0 u15=v15` and `v15.u15=kappa_new>=0` directly, so the range condition is
checked exactly even though the expanded matrix is rank deficient.

The embedding E maps an 11-vector v to

    (v_0+v_10,...,v_4+v_10, v_0,...,v_4, v_5,...,v_9).

Every nonorigin Gram ray is checked to lie in the image of E. The
builder's dual lift satisfies E^T u15=u11 exactly. Hence the embedded
matrix and coupling have the same quadratic mass as their functional
counterparts. This audit found no missing term in that transformation.

The added pointed pair term 2h cancels the origin shift 2(v_new-v_old)
exactly, while its nonorigin diagonal is zero. The global verifier checks
this coefficient identity and then recomputes the original pair
inequality using v_old. For full triangles, ignoring the nonpositive
pointed correction is conservative. For exterior triangles it checks the
correction explicitly by the interval method above. No new origin
constant comes from the pointed energy, so the global objective changes
only by `(kappa_new-kappa_old)/8` when the other certificate components
are fixed. The final verifier recomputes the complete objective rather
than relying on the builder's incremental formula.

## Audit conclusion

The pointed energy, conditional probability bounds, interval comparison,
and origin-shift integration are mathematically consistent. In particular,
the exterior correction has the required factor 1/2, its radial weights
are paired with the correct conditional probabilities, and the origin
shift is verified in the actual expanded Gram coordinates. Acceptance of
any numerical proposal still requires successful exact verification of
the scalar, all pairs, all full triangles, and all exterior branches.

The complete replay for `new_work/pointed/probe4_repair1` has now passed.
Its `VERIFIED.json` records objective 0.631993119255666 (decimal display),
strict threshold `3949957/6250000 = 0.63199312`, and 17 pointed caps. This
accepted threshold remains above the requested 0.63 target.

The separate general helper also replayed the two exterior branches that
needed pointed repair: [1,1] passed in 6421 nodes and [1,2] in 4211 nodes.
The accepted v1 helper and its completed global certificate were unchanged.
