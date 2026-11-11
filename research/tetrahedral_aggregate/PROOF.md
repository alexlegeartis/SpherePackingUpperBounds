# Four-point edge-kernel proof

## Definitions and exact input conditions

Consider a finite packing in R^22 whose distinct centers have distance at least 1. Fix q equal to 21/20 or 11/10. A close edge has squared length at most q. Let E and T count unordered close edges and all-close triangles.

For a close base edge AB, project each common neighbor X perpendicular to AB and normalize its transverse vector. Every incident squared distance lies in [1,q], and the transverse squared radius is at least 1-q/4>0.

The exact angular polynomial is

    Q(t) = sum_l c_l P_l^21(t),

where P_l^21 is the normalized Gegenbauer polynomial for vectors in the 21-dimensional transverse space, P_l^21(1)=1. Every c_l is nonnegative, c_0=1, and B=Q(1)=sum c_l>0. The independent verifiers reconstruct Q by an explicit rising-factorial expansion rather than trusting the stored power coefficients or producer recurrence.

The sharp projection lemma proves that if the opposite edge is long, its projected cosine is at most

    F = 1 - q/[2(q-1/4)].

Thus F=11/32 for q=21/20 and F=6/17 for q=11/10. Exact rational Sturm signs prove Q<0 on [-1,F], including the boundary.

The full compact proof establishes

    sum over all six base edges Q(t_e) < 0

for every six-tuple of squared edge lengths in [1,q]^6. In particular it covers every all-close tetrahedron. No geometric exclusions are used.

## PSD edge form

For any M>0, take a kernel on an origin plus the normalized projected common neighbors with

    J(0,0)=BM,
    J(0,u)=J(u,0)=-M,
    J(u,v)=(M/B)Q(u dot v).

The constant angular component is a rank-one extension with origin feature sqrt(BM) and neighbor feature -sqrt(M/B). All positive-degree components are PSD by the spherical addition theorem and c_l>=0. Hence every finite sum of J is nonnegative. Its neighbor diagonal is M because Q(1)=B.

For N common neighbors of one close edge,

    0 <= BM - MN + (2M/B) sum over unordered neighbor pairs Q(t).

Sum this inequality over all close edges. Each all-close triangle contributes exactly three to sum N. Each tetrahedron contributes once per base edge for which the other two vertices are common neighbors.

- If every edge is close, the sum over all six bases is negative by the full compact verification.
- If exactly one edge is long, only its opposite edge has both other vertices as common neighbors. Its Q value is nonpositive by the far sign.
- If at least two edges are long, no supported base remains: a supported base and its four incident edges would already require five close edges.

Therefore the entire tetrahedron contribution is nonpositive, and

    0 <= BME - 3MT,
    3T <= BE.

There is no missing factor two. This is an aggregate inequality, not a pointwise bound on each edge's common-neighbor count; B must not be floored. A proved upper rounding of B may be substituted.

Equality q is included in close edges. The far sign also includes the limiting endpoint, so support boundaries are covered conservatively. Compact verification includes the closed six-edge cube. The finite-set proof requires no limiting or periodic approximation.

## Exact compact interval verification

For squared base a and incident squared lengths u,v,w,z, write

    R=(u+v)/2-a/4-(u-v)^2/(4a),
    S=(w+z)/2-a/4-(w-z)^2/(4a).

If c is the opposite squared edge, the projected cosine is

    t = [R+S+(u-v-w+z)^2/(4a)-c]/[2 sqrt(RS)].

The independent vertex-coordinate version uses the equivalent formula with longitudinal coordinates (a+u-v)/(2 sqrt(a)) and (a+w-z)/(2 sqrt(a)). Interval positivity of denominators is checked independently on every box.

Automatic differentiation encloses every partial derivative of the sum over six edges. For a box with center m and coordinate radii r_i, the multivariate mean-value theorem gives the upper enclosure

    F(m) + sum_i [-r_i,r_i] partial_i F(box).

A box is accepted only when this enclosure, or a direct function enclosure, is strictly negative. Otherwise exact rational bisection preserves complete closed coverage. Nonfinite intervals, positive centers, excessive depth and uncompleted coverage fail the replay. Interval squares crossing zero use explicit multiplication where needed.

For q=11/10, fixed-sign partial derivatives also allow replacing a coordinate by its maximizing endpoint. This bounds the maximum over the whole box by the maximum on the retained face. The producer fixes all certified coordinates simultaneously; the independent implementation fixes one at a time and recomputes derivatives. Both are valid since the signs hold uniformly on the original box and every intermediate face.

The q=21/20 producer and independent midpoint proofs each use 3,090 leaves, depth 17. The independent vertex-origin proof uses 21,735 leaves, depth 21. The q=11/10 producer uses 9,954 leaves, 309 face reductions, depth 19; the independent implementation uses 9,954 leaves, 407 sequential reductions, depth 22. There are no excluded or unresolved boxes.

## Strict separation from affine and normal-harmonic envelopes

An affine envelope would require Q(t)<=a(t-1/3) throughout the sharp compact interval, with a>=0. For q=21/20 that interval is [17/59,3/8]. Exact endpoint evaluations force a>=107.34533753... and a<=75.76211606..., which is impossible. For q=11/10 the interval is [7/29,7/17], forcing a>=93.43885873... and a<=33.54012259.... The runner independently reconstructs and checks the exact rational contradictions.

NORMAL_HARMONIC_ENVELOPE_REDUCTION.md proves that every finite nonnegative combination of the unweighted face-normal harmonic envelopes reduces to this same affine family for 1<=q<=6/5, once the far sign is imposed. Therefore both full-geometry witnesses strictly enlarge that feasible family.

These are computer-assisted auxiliary proofs, not external peer review, formal proof-assistant verification, a claim of literature novelty, or a new density upper bound.
