# Attacking the ten remaining dimensions

**Only 58D supplies a new certified density improvement in this continuation.**
The other computations certify restrictions on particular routes or provide
exact counterexamples. The former conditional 92D candidate is rejected.
These results are excluded from the SOTA table.
All parents and targets remain at most96D.

**55D correction:** the former minimum 10/local-optimality claim below
is withdrawn. Its Gram contains `(-3Q,3)` of norm 6. A
[new exact correction and tilted 55D construction](../tilted55_20261116/PROOF.md)
gives center density 91,274.808704… without claiming a world record.

| Target | Result of this continuation | Consequence for the research |
| --- | --- | --- |
| 48 | The classical harmonic-theta argument proves local optimality of the extremal even unimodular parent | A small lattice Gram deformation cannot improve it |
| 54 | Reconstructed the parent from polynomial points; supplied8,748 minimum vectors, full perfection rank 1,485 and positive eutaxy | Local lattice Gram deformations are blocked |
| 55 | The former Gram has an explicit norm 6 vector, rejecting its minimum 10 claim | The former 55D local-optimality conclusion is withdrawn; the Q-hole geometry alone remains valid |
| 56 | Reconstructed MW56; a 2,436-ray orbit has full rank 1,596 and exact positive eutaxy | Local lattice Gram deformations are blocked; its cosets successfully yield58D |
| **58** | **Explicit minimum-10 lattice, center density161,833.931161…** | **New certified best-known bound; +6.9277% over the stronger classical section** |
| 63 | A universal ceiling rules out every inherited-minimum-only multi-translate codimension-one Ne64 guarantee | Using only the parent minimum cannot beat Vardy's reported Elkies guarantee |
| 64 | A ceiling covers every finite projection cluster for the Gamma72 / 4E8 normal | That inherited-minimum-only family cannot beat Ne64; other normals remain open |
| 72 | The classical harmonic-theta argument proves local optimality of the extremal even unimodular parent | A small lattice Gram deformation cannot improve it |
| 80 | Reconstructed MW80 from 160 polynomial sections; a 4,920-ray orbit has full rank 3,240 and exact positive eutaxy | Local lattice Gram deformations are blocked |
| 92 | Three additional exact shortest vectors reject the proposed deformation and complete the local-optimality proof | The former 31.1436% candidate is rejected; a different construction is required |

“Local optimality” here concerns lattice Gram deformations near the stated
parent. It is not a claim of global optimality among all sphere packings.
A different lattice, a distant deformation, or a periodic construction
can still improve these dimensions.

## Exact certificates for 54D–56D and 80D

The source-to-Gram reconstruction in 54D uses the curve
`y²=x³+x+t^28` over F_(3^6), with field polynomial
`X^6-X^4+X²-X-1`. For a rational x-coordinate N/D, the explicit tripling
formula is

```
x(3P)=x(P)^9+t^84-t^28.
```

Since3P is narrow, Leterrier's Lemma3.2.15, Lemma3.2.17 and
Corollary3.2.18 give

```
q(P)=max(deg(N),deg(D)+28/3).
```

After cancellation the displayed degrees cannot tie. The group-addition
formula and this identity compute the height Gram of 124 saved polynomial
points. Their integral coordinate matrix spans Z^54, and the resulting
Gram has determinant 3^25. Thus the published narrow minimum 10 applies
to this exact reconstructed parent.

The point Q=(0,t^14) has norm 28/3. Its coordinates q satisfy3q integral,
q nonintegral, and the directly computed pairings to the polynomial
points. Adding(q,sqrt(2/3)) gives the former 55D Gram of determinant
2·3^24, whose actual minimum is6. Its full Gram is retained as a
negative certificate; it is not a minimum 10 lattice.

For the 55D hole, translations of the base variable produce closest points
explicitly. For b in F_(3^6), set

```
z_b(t)=d+b^27 t-b^81 t³+b^243 t^9, d³+d=b^28, d in F_27.
```

Then `(z_b(t),(t+b)^14)` lies on the curve and has canonical norm 28/3.
Together with the three x-translations by constants c with c³+c=0,
these give2,187 closest vectors in the same Q coset. The supplied vectors
are checked directly by their norms and integral differences from q;
regenerating the discovery group is unnecessary.

The finite checks are:

* The 54D minimum vectors have dyad rank 1,485 modulo 65521 and moment
  `sum(v^T v) G54=1620 I`.
* The 2,187 closest vectors span all 54 dimensions, sum to zero, and have
  moment `sum(w^T w) G54=378 I`.
* In the former 55D Gram, weight 1/2187 on the supplied horizontal norm10 vectors and
  weight 1/1458 on the supplied lifted closest vectors give exactly
  `G55^-1`. All weights are strictly positive. Their dyad rank is1,540.
  **These are not minimum vectors of G55:** its minimum is6, so this
  calculation does not prove local lattice optimality.
* The 56D and 80D minimal orbits have full dyad rank and equal-weight
  moments435 I and 861 I, respectively. Their parent Grams are bound to
  actual unitary polynomial sections as in the main proof.

For54D,56D and 80D, full rank modulo a prime proves full rank over Q. Voronoi's criterion
then gives local lattice optimality from perfection and strict eutaxy.
No complete classification of all minimum vectors is needed.

## The classical local obstruction in 48D and 72D

For an extremal even unimodular lattice of rank 48 and minimum 6, a harmonic
theta series of degree 2 or4 begins at q^3. Its modular weight is26 or28,
so the level-one valence bound forces it to vanish. In rank 72, minimum 8,
the corresponding weights are38 and 40 and the initial order is4; the
same argument applies. The minimum shells are therefore spherical
4-designs. Their second and fourth moments imply eutaxy and perfection,
which rule out an improving small lattice Gram deformation.

## A route ceiling in 63D

Let the parent be Ne64, center density3^16 and minimum 12. For a primitive
dual normal of squared norm q and M>=2 distinct integer labels on its
projection line, the squared diameter is at least(M-1)^2/q. The antipode
ratio is consequently bounded by

```
2M/(M-1) * sqrt(x) * (1-x)^(63/2),
x=(M-1)^2/(12q).
```

The continuous maximum occurs at x=1/64, and M/(M-1)<=2. Hence every
such inherited-minimum guarantee satisfies

```
delta63² <=3^32/4*(63/64)^63,
delta63 <=13,105,940.958282… .
```

An exact integer-power check puts this below2^23.92, the lower endpoint
of the Elkies guarantee reported by Vardy(1999). This bounds the inherited-minimum formula for every normal and every
multi-point cluster. An antipode with a stronger actual minimum would
require a separate minimum proof and is not ruled out.

## A full fixed-normal family ceiling in 64D

For a 4E8 normal in Gamma72, the projection label lattice is E8 scaled
by 1/2. Its squared norms are multiples of1/2 and its minimum is 1/2.
The bare inherited-minimum section density is 2^24. A cluster of M
points with squared diameter beta gives the guarantee

```
2^24 * M * (1-beta/8)^32.
```

All possible diameters below 8 are covered as follows:

* beta=0 implies M=1.
* beta=1/2 gives an equilateral simplex, so M<=9.
* beta=1 gives a two-distance set. The polynomials
  `(||x-p||²-1/2)(||x-p||²-1)` are independent by evaluation at the
  cluster points. They lie in the span of 1, the 8 linear monomials,
  the 36 quadratic monomials, the 8 functions `||x||²*x_j`, and `||x||^4`.
  Thus M<=54.
* beta=3/2 gives M<=2^8 by injection into the label lattice modulo 2.
* For 2<=beta<9/2, injection modulo 3 gives M<=3^8.
* For 9/2<=beta<8, injection modulo 4 gives M<=4^8.

The residue arguments use the fact that a nonzero difference in j times
the label lattice has norm at least j²/2. Taking the worst endpoint in
each interval gives the exact ratio bounds

```
1, 9*(15/16)^32, 54*(7/8)^32,
256*(13/16)^32, 6561*(3/4)^32, 65536*(7/16)^32.
```

The largest is the nine-point value, already attained by the archived
finite cluster. Therefore **every inherited-minimum-only guarantee in
this fixed-normal family** is at most 19,144,465.698823…, below Ne64's
43,046,721. This does not bound an antipode's unknown actual minimum or
constructions using other normals.

## 92D: exact counterexample rejects the deformation

**The former conditional candidate is now rejected.** The data contain two integer
symmetric forms H0,H1. With G the verified MW92 Gram, both satisfy

```
A_i=G^-1 H_i integral, A_i²=3I, trace(A_i)=0,
v H_i v^T=0 for all 6,486 rays in the supplied minimum-vector orbit.
```

The exact dyad rank of this orbit is4,276, leaving two directions in the
4,278-dimensional symmetric-form space. Both null forms are independent
and verified exactly. This gives the explicit candidate

```
G_candidate=G+H0/16,
det(G_candidate)/det(G)=(253/256)^46.
```

The candidate is positive definite. Since3/256<1/81, it satisfies
`G_candidate>(8/9)G`. Every parent vector of norm at least18 therefore
still has norm greater than 16. The original 6,486 rays retain norm 16.
However, an integral involution
`J=(H0 G^-1-H1 G^-1)/3` splits the lattice into two eigenspaces.
Reduction of the image of `I-J` found the saved integer vector w in
[counterexample92.json](frontier/counterexample92.json). Exactly,

```
w G w^T=16, w H0 w^T=-24,
w (G+H0/16) w^T=29/2 <16.
```

Thus the smaller determinant does not supply the stated density. The
earlier finite BKZ check missed this witness; it was never a minimum
certificate. The candidate was never included in the SOTA envelope.

The same computation settles local optimality. Put
`C=(H0 G^-1)(H1 G^-1)/3`. The verifier proves C integral,
`C G C^T=G` and `C²+C+I=0`. The three norm-16 vectors w,wC,wC²
have null-form values spanning R² and summing to zero. Therefore their
positive cone is all R², and they raise the perfection rank to 4,278.

Here is the local determinant argument without assuming that any shell
is complete. Normalize a perturbation to keep minimum 16. For an
admissible Gram increment D, the original minimum
vectors obey `v D v^T>=0`. Their exact positive-eutaxy identity gives
`trace(G^-1 D)>=0`. Equality forces every original dyad pairing to
vanish, so D lies in the span of H0,H1. The three additional vectors
then force D=0, by the positive-spanning null values. On the compact
unit slice of this closed admissible cone, trace is therefore strictly
positive. The first derivative of log determinant dominates its bounded
quadratic remainder for sufficiently small increments. Hence G is
strictly locally optimal among lattice Gram deformations.

The [next continuation](../creative_sections_20261116/) supplies stronger
certified bounds in 57D, 78D and 79D; it makes no new 92D density claim.

## Replay

```sh
python research/cm58_20261116/verify_frontier.py
```

This checks the finite-field reconstructions, exact matrix identities,
positive eutaxy, modular ranks, the 63D ceiling and the rejected 92D
deformation. The rank checks take longer than the main58D saved-report audit.
The result explicitly distinguishes the certified density, local
obstructions and rejected candidates.
