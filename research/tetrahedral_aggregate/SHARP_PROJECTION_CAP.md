# Sharp common-neighbor projection cap for a short edge

This is a geometric lemma, not a packing-density certificate. The derivative identities and rational lower bounds below were checked independently by symbolic expansion.

## Statement

Let 1 <= q <= 6/5. Let A,B,X,Y be Euclidean points with squared base length |A-B|^2 in [1,q] and all four incident squared lengths |A-X|^2, |B-X|^2, |A-Y|^2, |B-Y|^2 in [1,q]. Suppose |X-Y|^2 >= c, where 1 <= c <= q. Project X and Y orthogonally onto the subspace perpendicular to AB, using the base midpoint as longitudinal origin. Their transverse vectors have nonzero lengths. The cosine t between them satisfies

    t <= 1 - c / (2(q-1/4)).

In ambient dimension at least three, this bound is attained when the base has squared length 1, both neighbors lie in its midpoint hyperplane, their four incident squared distances equal q, and their mutual squared distance equals c.

Consequently, the compact-tetrahedron cap is 1-1/[2(q-1/4)], and the one-long-edge cap is 1-q/[2(q-1/4)]. If the long edge is required to be strictly greater than q, the latter is a supremum. At q=6/5 the two caps are 9/19 and 7/19.

## Coordinates and monotonicity

Put a=|A-B|^2 and b=sqrt(a). Write z for a neighbor's signed longitudinal coordinate relative to the base midpoint, and R for its transverse squared radius. Its two incident squared distances are

    R + (z+b/2)^2,    R + (z-b/2)^2.

Their range gives |z| <= (q-1)/(2b). For fixed x=|z|, the maximal permissible transverse squared radius is

    U - bx - x^2,    U=q-a/4.

All actual and maximal permissible transverse squared radii lie in

    L0=1-q/4 <= R <= U0=q-1/4.

For two neighbors, let delta be their longitudinal difference. Then delta^2 <= (q-1)^2. The separation condition bounds their cosine by

    F(R1,R2) = (R1+R2-C)/(2 sqrt(R1 R2)),
    C=c-delta^2.

The derivative with respect to R1 has the sign of R1-R2+C, and likewise for R2. Since

    C >= 1-(q-1)^2,
    U0-L0 = 5(q-1)/4,

both signs are positive: their uniform lower bound is 71/100. Thus each radius may be replaced by its individual maximum. For fixed absolute longitudinal coordinates x,y, opposite longitudinal signs maximize delta^2, giving delta^2=(x+y)^2.

## Balancing the longitudinal coordinates

Set p=x+y and r=xy. The resulting cosine bound is N/(2 sqrt(J)), where

    N = 2U-c-bp+2r,
    J = U^2-Ubp-Up^2+(2U+a+bp)r+r^2.

Here J is the product of the two positive maximal transverse squared radii. At fixed p, the derivative in r has the sign of

    4J-N J_r
      = 2U(c-a)+ac+(-4U+c+a)bp
        +(-4U+a)p^2+2r(a+c+2bp).

Use U<=19/20, a,c in [1,6/5], p<=1/5, bp<=1/5 and r>=0. The expression is bounded below by

    1 - 19/50 - 9/25 - 14/125 = 37/250 > 0.

Therefore the cosine increases with r. For fixed feasible p, r is maximized at x=y=p/2; this balanced choice remains feasible because each coordinate is bounded by (q-1)/(2b).

Write the common value as z. The cosine becomes

    1 - (c-4z^2)/(2(U-bz-z^2)).

The derivative of the subtracted ratio has numerator

    cb + (2c-8q+2a)z + 4bz^2.

Because z<=1/10, this numerator is at least

    1 - (28/5)(1/10) = 11/25 > 0.

Thus the cosine is maximized by z=0. The remaining expression is 1-c/(2U), which is maximized at a=1 because U=q-a/4 decreases with a. This proves the bound.

## Equality and use in the tetrahedral count lemma

Take a base of length 1. Choose two transverse vectors of squared length q-1/4 with cosine 1-c/[2(q-1/4)]. Put both longitudinal coordinates at zero. All four incident squared lengths are q, and their mutual squared length is c. The required cosine lies strictly between -1 and 1, so this construction exists in a two-dimensional transverse subspace.

The compact lower cosine bound 1-q/[2(1-q/4)] from the earlier argument remains valid. Hence the tetrahedral affine-envelope certificate may use the compact interval

    [1-q/(2(1-q/4)), 1-1/(2(q-1/4))]

and the far interval

    [-1, 1-q/(2(q-1/4))].

This sharpens only the geometric ranges. A usable common-neighbor aggregate bound still requires a Gegenbauer-positive polynomial with certified signs/envelope on these intervals, and a density result requires the subsequent full packing certificate.
