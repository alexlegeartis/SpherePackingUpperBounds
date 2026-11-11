# Reduction of tetrahedral normal-harmonic envelopes

## Result

For 1 <= q <= 6/5, replacing the affine compact envelope by any finite nonnegative combination of the unweighted face-normal harmonic envelopes below does **not** enlarge the feasible set of angular polynomials Q, when the same far sign is imposed. The entire family reduces exactly to the original nonnegative affine envelope.

This is a structural simplification, not a new packing-density bound. It explains why higher-degree envelope searches can return the same objective with an apparently additional degree-five component.

## Validity of the proposed harmonic constraints

A tetrahedron has four outward unit face normals n_1,...,n_4 in its three-dimensional affine span. Its six internal dihedral cosines are t_e=-n_i·n_j over the six pairs of adjacent faces. For normalized three-dimensional Gegenbauer polynomials, which are ordinary Legendre polynomials P_l with P_l(1)=1, the addition theorem gives

    0 <= sum_(i,j) P_l(n_i·n_j)
       = 4 + 2 sum_(i<j) P_l(n_i·n_j).

Using parity, define

    E_l(t) = (-1)^(l+1) P_l(t) - 1/3.

Then sum_e E_l(t_e) <= 0. Degenerate tetrahedra follow by continuity when the supported faces have nonzero transverse radius. Therefore

    Q(t) <= sum_l a_l E_l(t),    a_l >= 0,

on an interval containing all six compact cosines is sufficient for the required nonpositive compact tetrahedral sum.

The sharp compact interval from SHARP_PROJECTION_CAP.md lies inside [0,1/2] for 1<=q<=6/5. The far sign is imposed on [-1,c_far], with

    c_far = 1 - q/[2(q-1/4)] >= 1/3.

These two facts drive the reduction.

## Only degrees one and two survive on the compact range

The first two envelopes are

    E_1(t) = t-1/3,
    E_2(t) = 1/6 - 3t^2/2.

On 0<=t<=1/2:

- E_3<=-1/3, because t(5t^2-3)<=0.
- Writing x=t^2, E_4=(-105x^2+90x-17)/24 is increasing on 0<=x<=1/4 and has maximum -17/384.
- E_6=(-693x^3+945x^2-315x-1)/48 <= -1/48, because 945x^2-315x=315x(3x-1)<=0.
- Degree five is dominated by degree two:

    (5/9)E_2(t)-E_5(t)
      = (3t-1)^2 [92+147t-126t^2-189t^3]/216 >= 0.

  The bracket is at least 92-126/4-189/8=295/8>0 on this interval.

All degrees l>=7 are strictly negative there. A self-contained bound follows from the Legendre differential equation. Set

    u(theta)=sqrt(sin theta) P_l(cos theta),
    A(theta)=(l+1/2)^2 + 1/(4 sin^2 theta).

Then u''+A u=0. On theta in [pi/3,pi/2], A'<=0, so the Sonin quantity

    S(theta)=u(theta)^2 + u'(theta)^2/A(theta)

has derivative -A' u'^2/A^2>=0. Thus u(theta)^2<=S(pi/2).

For even l, S(pi/2)=P_l(0)^2. For odd l,

    S(pi/2)=l^2 P_(l-1)(0)^2 / [l(l+1)+1/2].

For l=7 this is 1225/14464<2/21. The exact formula |P_(2k)(0)|=binomial(2k,k)/4^k decreases with k, since successive values have ratio (2k+1)/(2k+2). Hence every even l>=8 satisfies P_l(0)^2<=(35/128)^2<2/21. For odd l>=9 the displayed ratio is strictly smaller than P_(l-1)(0)^2 and has the same bound. Since sin(theta)>=sqrt(3)/2>6/7, it follows that

    |P_l(t)|^2 < (2/21)/(6/7) = 1/9

for all l>=7 and 0<=t<=1/2. Hence E_l(t)<0 for all these degrees.

Consequently, dropping the strictly negative envelopes and replacing a_5 E_5 by (5a_5/9)E_2 only raises the envelope. Every feasible Q remains feasible with just E_1 and E_2, after replacing a_2 by a_2+5a_5/9.

## Reduction to the original affine envelope

The exact identity

    E_2 = -E_1 - (3/2)E_1^2

gives

    a_1 E_1 + a_2 E_2
      = (a_1-a_2)E_1 - (3a_2/2)E_1^2.

If a_1>=a_2, this is bounded above by the original affine envelope (a_1-a_2)E_1, whose slope is nonnegative.

If a_1<a_2, the expression is nonpositive whenever t>=1/3. For t<=1/3, the far condition already gives Q(t)<=0 because c_far>=1/3. Thus Q is nonpositive throughout the compact range as well, and the original affine envelope with slope zero suffices.

The reverse inclusion is immediate by taking a_1 equal to the original affine slope and every other coefficient zero. Therefore the projected feasible-Q sets are identical. This conclusion holds regardless of the degree used to represent Q.

## Scope

The proof concerns finite nonnegative combinations of these unweighted, one-variable face-normal harmonic envelopes, together with the stated far sign and q range. It does not rule out weighted face-normal identities, constraints depending jointly on multiple dihedral angles or edge lengths, a different support range, or other four-point constructions. No optimizer status or sampling assumption enters the reduction.
