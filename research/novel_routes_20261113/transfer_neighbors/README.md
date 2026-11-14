# New parent lattices, quotient holes, and a complete layering obstruction

**No new density record was obtained.** This directory pursues a different
construction route from the earlier36D deformation and P48 section searches:
lamination over seven32D catalogue lattices of minimum6 and determinant2^16.
It contains a new exact finite result and a resulting family-wide obstruction,
plus explicitly labelled floating discovery searches.

## Primary inputs and benchmark

The seven Gram matrices are archived verbatim from the Nebe–Sloane catalogue's
[32D level-two data](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/ExtLatDat/2_dim32.dat).
The source is `sources/2_dim32.dat`; its hash and flags are in
`source_metadata.json`. All seven have published minimum6 and determinant65536.
The bundle marks indices1 and2 as **Modular=false**. They are therefore not
described here as seven2-modular lattices. The finite computations verify the
required level-two structure directly for every matrix.

The minimum6 is inherited from the catalogue and is not re-proved here.
Positive definiteness, determinants, and all saved changes of basis are checked
exactly. The32D center density is3^16/2^24. The33D and34D benchmark is
3^(33/2)/2^25, as in
[Henry Cohn's table](https://cohn.mit.edu/sphere-packing/).
[SPLAG's third-edition preface](https://neilsloane.com/doc/splag3.pdf)
describes Elkies's lamination over putative norm-four deep holes of these32D
parents. A later source lookup did not identify a replacement33D/34D record;
that lookup is not an exhaustive priority guarantee.

## Why holes farther than4 would matter

Let L be one of these32D parents, so covol(L)=256 and min(L)=6. For any real
shift x, consider the33D periodic point set

\[
(L\times 2h\mathbb Z)\ \cup\ ((L+x)\times(h+2h\mathbb Z)).
\]

Its squared minimum is at least
min(6,4h^2,h^2+d(x,L)^2). Unlike a lattice extension, this construction does
not require2x to belong to L. If4<d(x,L)^2<=9/2, choose
h^2=6-d(x,L)^2. Its center density divided by the33D benchmark is

\[
\sqrt{\frac{2}{6-d(x,L)^2}}>1.
\]

Consequently a rigorously certified hole farther than4 would immediately give
a new33D bound. No such hole was found.

## Exact result: all65536 dual quotient classes are covered

For every parent, the quotient A=L*/L has65536 classes and exponent2. The saved
certificate gives one vector in every class, with the following **witness**
squared norms:

| Squared norm | Number of classes |
|---:|---:|
|0|1|
|3|32640|
|4|32895|

This proves d(x,L)^2<=4 for every x in L*. It does not prove an upper bound4
for the full covering radius, nor does it cover every point of(1/2)L.

The proof is finite and does not depend on the completeness of any enumeration.
Write G for the primal Gram and D=2G^(-1). Both D and its diagonal are integral
and its diagonal is even. For a saved unimodular T, H=TDT^t is an integral
Gram. A coefficient row w represents the dual vector with primal coordinates
wTG^(-1), and has squared norm wHw^t/2. Its class is determined by the parity
of wTD. Sixteen selected parity coordinates provide its stored label.

Each `dual_quotient_i.npz` contains G,D,H,T, the label coordinates, and65536
integer witness rows. The independent verifier checks the exact matrix
identities, det(T)=+/-1, all witness norms, and every distinct class label.
Since |L*/L|=det(G)=65536, all classes are accounted for. Across the seven
parents these compressed certificates occupy5,369,414bytes.

The discovery method enumerated a pool of short vectors in2L*, then formed
short pair sums and differences. The labels of a pair sum are the XOR of the
individual labels. Random pairing followed by directed searches for uncovered
labels filled every class. Only the final exact witnesses are needed for the
proof.

## Stronger consequence: no33D gain from arbitrary layers with dual offsets

The finite certificates exclude more than two-layer constructions. Fix one
of the seven parents. Consider **any packing formed from full translates of L
in parallel32D layers**, with every translation offset in L*, and arbitrary
real heights. The layers may be nonuniform, nonperiodic, or share a height when
their point sets are disjoint. Allow its overall squared minimum mu to be any
positive value. Necessarily mu<=6 because every full layer contains L.

The quadratic form q(a)=||a||^2 mod2 is well-defined on A. The verifier checks
this quadratic form on all65536 classes and verifies that its alternating
polar form is nondegenerate of rank16. The norm-three witnesses lie in every
q=1 class; the other classes have witnesses of norm at most4.

Any set of classes with q(a-b)=0 for every pair has at most256 elements. After
translation it contains0; its span is isotropic for the nondegenerate
alternating form, hence has dimension at most8.

Let I be the upper linear intensity of layers, counting different translates
at the same height separately. Three necessary bounds follow:

1. A repeated offset has vertical separation at least sqrt(mu), so
   I<=65536/sqrt(mu).
2. For mu>3, every open interval of length sqrt(mu-3) contains at most256
   layers. Otherwise two offsets have q-difference1, whose norm-three witness
   would produce a pair closer than sqrt(mu). Thus
   I<=256/sqrt(mu-3).
3. For mu>4, every two distinct layers have vertical separation at least
   sqrt(mu-4), because every offset difference has a witness of squared norm
   at most4. Thus I<=1/sqrt(mu-4).

The center density satisfies

\[
\delta_{33}^2\le\frac{\mu^{33}I^2}{2^{82}}.
\]

The minimum of these three intensity bounds changes at

\[
\mu_1=\frac{3\cdot65536}{65535}=\frac{65536}{21845},\qquad
\mu_2=4+\frac1{65535}=\frac{262141}{65535}.
\]

The three squared-density bounds are respectively

\[
\frac{\mu^{32}}{2^{50}},\quad
\frac{\mu^{33}}{2^{66}(\mu-3)},\quad
\frac{\mu^{33}}{2^{82}(\mu-4)}.
\]

The first is increasing. For mu^33/(mu-a), the derivative has the sign of
32mu-33a, so every interior critical point is a minimum. It suffices to check
mu1,mu2,6. Their squared-density upper bounds, divided by the squared33D
benchmark3^33/2^50, are respectively

\[
0.333496134735\ldots,\quad0.202556611151\ldots,\quad1.
\]

The comparisons are performed as exact rational inequalities. Therefore no
packing in this entire specified family can improve the33D benchmark, even
by lowering its minimum or using nonlinear, irregular, or many-layer offsets.
This does not exclude arbitrary real offsets, changes to the32D period lattice,
deletion of selected points within layers, or constructions outside this family.

`layer33_theorem_checked.json` records the exact endpoint comparisons. The
derivation was independently audited by the coordinating research agent.

## Floating searches and rejected apparent gains

The initial half-lattice search made17507 CVP calls across seven parents.
1548 calls hit their node cap. Capped calls only provide upper witnesses for
the distance and cannot certify a large hole. Poor bases created apparent
scores4.5 and5; explicit exact witnesses reduce these to3.5 and3 respectively.
Both rejecting vectors are retained in `capped_candidate_rejections.json`.
After basis reduction, an additional2371-call run on parent4 had no caps and
found only squared distance4. A464-call preliminary run is retained separately.

The continuous Voronoi LP/CVP searches used107 starts and2914 CVP calls.
Best floating squared distances were:

| Catalogue index | Best continuous score | Starts | CVP calls |
|---:|---:|---:|---:|
|0|3.8673469388|14|380|
|1|4.0000000000|14|433|
|2|3.8632576846|16|435|
|3|4.0000000000|16|417|
|4|3.8672000000|14|416|
|5|4.0000000000|16|414|
|6|3.8533057851|17|419|

These are discovery searches, not covering-radius proofs. Some starts were
biased near previously found norm-four half-lattice holes. The scripts save
all improving numerical points. No score exceeded4 by the candidate threshold.

## A separate36D antipode calculation

For a one-dimensional normal projection of the existing36D lattice, take m>=2
distinct projected levels and use the standard inherited minimum guarantee
8-beta. Writing q for the squared dual normal and using diameter at least
(m-1)/sqrt(q), this density guarantee obeys

\[
\delta_{35}\le\delta_{36}\frac23
\left(\frac{35}{36}\right)^{35/2}
=2.1506693474\ldots<2\sqrt2.
\]

This bounds that standard guarantee only. It does not bound an independently
proved larger minimum of the projected point set, and it does not cover
singleton sections. `verify_transfer.py` checks the exact squared inequality.

## Replay

Exact replay needs only NumPy and SymPy and does not download anything:

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/transfer_neighbors/verify_transfer.py
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/transfer_neighbors/verify_dual_quotient.py
```

Discovery additionally needs SciPy, fpylll/cysignals, and the repository's
`research/density_breakthroughs_20261113/bin/shell_discovery`. The source
`half_holes.cpp` compiles with a standard C++17 compiler. `prepare_inputs.py`
reconstructs the source inputs and reduced bases. `dual_quotient.py i` recreates
the finite witness for parent i; basis reduction and floating enumeration can
change its particular witnesses without changing the exact verification test.
Large raw short-vector pools are ignored and are not needed for the result.
