# A 604-edge antipodal local energy for the broad annulus

**Subsequent application status:** the four-graph proposal at E weight
\(4\times10^{-7}\) was rejected after its full run certified 110 of 112
branches. Both remaining packets contain actual positive exterior values,
including strict radial interior witnesses. The 604-edge local theorem
below remains valid; the proposed global application is rejected. The
broader graph with a degree-nine projection proof is a separate theorem.

## Result and scope

The following local packing theorem is proved with exact arithmetic. Let
\(C\subset\mathbb R^4\) be a finite set of neighbors of a packing center,
with pairwise distances at least one. Retain the neighbors whose radii lie
in

\[
I=\left[\frac{33}{28},\frac74\right].
\]

Connect two retained neighbors when the cosine of their directions lies in

\[
T=\left[-1,-\frac{183}{200}\right].
\]

The resulting graph has at most **604 ordered edges**. Every selected edge
has length strictly larger than \(7/4\). Consequently, for any
\(\delta_E\ge0\),

\[
\mathcal E_E(C)=\delta_E\bigl(604-e_E(C)\bigr)\ge0
\]

is a valid local energy whose pair term and kernel diagonal vanish. Here
\(e_E(C)\) counts both orientations of each ordinary graph edge.

This theorem is an unconditional geometric statement. Applying it to a
particular scalar function and kernel still requires the complete global
certificate. The proposed four-graph candidate has a conditional density
threshold of \(0.63013902009\); the local result and finite point checks in
this note do not by themselves certify that threshold.

## 1. The radial extension has no cost in directional separation

Write \(x=ru\) and \(y=sv\), with unit vectors \(u,v\). Packing separation
gives

\[
u\cdot v\le h(r,s):=\frac{r^2+s^2-1}{2rs}.
\]

For fixed \(s>0\),

\[
\frac{\partial h}{\partial r}
=\frac{r^2-s^2+1}{2r^2s}.
\]

An interior critical point is a minimum, so the maximum over a radial
rectangle occurs at one of its four corners. Taking
\(L=33/28\) and \(R=7/4\), those corner values are exactly

\[
h(L,L)=\frac{697}{1089},\qquad
h(L,R)=h(R,L)=h(R,R)=\frac{41}{49}.
\]

Thus all distinct directions in the broad annulus satisfy

\[
u\cdot v\le s_0:=\frac{41}{49}<1.
\]

The choice \(33/28\) comes from \(R-1/R\). For a fixed upper radius
\(R=7/4\), this is the smallest lower endpoint that keeps the mixed radial
corner at the same angular cutoff as the outer corner. Extending the
eligible vertices from the outer shell \([3/2,7/4]\) down to \(33/28\)
therefore preserves the spherical count used below. It also includes the
middle-shell antipodal counterexample near radii 1.415 and 1.419.

## 2. A fresh exact count of 151 directions

On \(S^3\), use the normalized zonal polynomials

\[
Z_\ell(t)=\frac{U_\ell(t)}{\ell+1},\qquad Z_\ell(1)=1,
\]

where \(U_\ell\) denotes the Chebyshev polynomial of the second kind. For
every finite list of unit vectors, the addition formula gives

\[
\sum_{i,j} Z_\ell(u_i\cdot u_j)\ge0.
\]

The certificate supplies an explicit rational polynomial of degree 20,

\[
p(t)=\sum_{\ell=0}^{20}a_\ell Z_\ell(t),
\]

with

\[
a_0=\frac{999999}{1000000}>0,\qquad a_\ell\ge0\quad(\ell\ge1),
\]

and proves on the entire closed interval that

\[
p(t)\le0\qquad\text{for }-1\le t\le\frac{41}{49}.
\]

The proof reconstructs the monomial polynomial from its rational harmonic
coefficients, converts it to Bernstein coefficients, and subdivides using
exact rational de Casteljau steps. The complete sign proof uses 93 nodes,
47 accepted leaves, and maximum depth 12. No pre-existing numerical count
or previously asserted neighbor count is used.

For \(n\) retained directions, the addition formula and off-diagonal sign
give

\[
a_0n^2\le\sum_{i,j}p(u_i\cdot u_j)\le n\,p(1).
\]

If \(n>0\), it follows that

\[
n\le\frac{p(1)}{a_0}
=\frac{151196972764393}{999999000000}<152.
\]

Since \(n\) is an integer, **\(n\le151\)**. The case \(n=0\) is immediate.
The exact polynomial coefficients are included in the graph certificate,
and the graph checker freshly proves the count as part of its own replay.

## 3. The antipodal graph has degree at most four

Fix a graph vertex with direction \(u\), and set \(p=-u\). If \(v\) is a
selected neighbor of \(u\), then

\[
v\cdot p\ge a:=\frac{183}{200}.
\]

Project its direction into the three-dimensional space perpendicular to
\(p\):

\[
w_v=v-(v\cdot p)p\in p^\perp\cong\mathbb R^3.
\]

For distinct selected neighbors \(v_i,v_j\), the preceding directional
separation yields

\[
\begin{aligned}
w_{v_i}\cdot w_{v_j}
&=v_i\cdot v_j-(v_i\cdot p)(v_j\cdot p)\\
&\le\frac{41}{49}-\left(\frac{183}{200}\right)^2\\
&=-\frac{961}{1960000}<0.
\end{aligned}
\]

At most four vectors in \(\mathbb R^3\) can have pairwise strictly negative
inner products. Here is a direct proof of the fact being used. If there
were five such vectors, they would have a nontrivial affine dependence
\(\sum_i\lambda_iw_i=0\), \(\sum_i\lambda_i=0\). The positive and negative
coefficients give two nonempty positive combinations with the same value:

\[
w=\sum_{\lambda_i>0}\lambda_iw_i
=\sum_{\lambda_j<0}(-\lambda_j)w_j.
\]

Taking the inner product of these two representations makes \(\|w\|^2\)
a strictly negative sum of cross inner products, a contradiction.

Therefore each graph vertex has degree at most four. Together with the
fresh count, this proves

\[
e_E(C)=\sum_{v\in C_I}\deg(v)\le4|C_I|\le4\cdot151=604.
\]

The projection argument needs the strict inequality \(a^2>41/49\), rather
than any special property of the fraction \(183/200\). This allowed the
original exploratory band ending at \(-11/12\) to be widened at the same
capacity. The wider rational endpoint gives more room around the observed
positive pocket.

## 4. The correction occurs only on exterior triangles

For a selected pair of radii \(r,s\ge L\) and cosine \(t\le-a<0\),

\[
|x-y|^2=r^2+s^2-2rst
\ge 2L^2(1+a)
=\frac{417087}{78400}.
\]

The exact gap above the squared support radius is

\[
\frac{417087}{78400}-\left(\frac74\right)^2
=\frac{176987}{78400}>0.
\]

Hence a selected graph edge cannot occur in a full triangle whose three
side lengths are at most \(7/4\). In an exterior configuration, only the
anchor having both supported radial edges can contribute. There is at
most one selected anchor for E.

Define the symmetric kernel

\[
K_E(x,y)=-\delta_E\,
\mathbf 1_{\{|x|,|y|\in I,\ \widehat x\cdot\widehat y\in T\}}
\quad(x\ne y),
\]

with zero diagonal and zero linear term. The sum of this kernel over
ordered distinct neighbors is \(-\delta_E e_E(C)\), and the origin
constant \(604\delta_E\) gives the proved nonnegative local energy.

The graph adds no pair cost, and leaves the complete existing pair gauge
unchanged. Its entire objective charge is the origin constant, divided by
eight in the current density normalization.

## 5. Exact charge and integration with A, B, and D

The selected E weight is

\[
\delta_E=\frac{1}{2500000}=4\times10^{-7}.
\]

Its origin increment and density increment are, respectively,

\[
604\delta_E=\frac{151}{625000}=0.0002416,
\qquad
\frac{604\delta_E}{8}
=\frac{151}{5000000}=0.0000302.
\]

The four independent graph energies have the following charges.

| Graph | Ordered capacity | Exact weight | Density increment |
|---|---:|---:|---:|
| Wider A | 862 | \(3/10000000\) | 0.000032325 |
| B | 1102 | \(1/10000000\) | 0.000013775 |
| D | 6252 | \(1/500000000\) | 0.000001563 |
| Widened E | 604 | \(1/2500000\) | 0.000030200 |
| Total | — | — | **0.000077863** |

The unchanged robust base objective is strictly below
\(63006115709/10^{11}\). Adding these exact costs gives the conditional
threshold

\[
\frac{63006115709}{10^{11}}+\frac{77863}{10^9}
=\frac{63013902009}{10^{11}}
=0.63013902009.
\]

The final integration explicitly supplies the A/B/E weight tuple
\((3/10^7,1/10^7,1/2500000)\), with capacities \((862,1102,604)\).
The frozen helper retains its exploratory default E weight of
\(3/12500000\); the explicit weight argument overrides that default. D retains
its separate helper and weight. Every full-triangle graph decrement is
zero. For exterior boxes, the helper subtracts a graph weight only if the
whole radial and cosine intervals lie inside that graph's closed support
(with A's middle radial class open at its lower endpoint).

A box that straddles a support boundary receives zero lower credit from
that graph until subdivision resolves it. The safe exterior upper
correction is the sum of the graph weights. At an exact point, the
membership predicates are exact, so the correction upper bound can equal
the exact lower bound. This last identity prevents a coarse upper bound
from masking a positive exact counterexample during diagnostic rejection.

## 6. Point replay, the rejected smaller weight, and the remaining gate

**The smaller E weight \(2.4\times10^{-7}\) was rejected.** After the
initial 21-point screen, a further middle/middle diagonal peak was found near
\(r=s=1.4142120740794129,z=0\), with numerical value about
\(3.4704299\times10^{-7}\). A separate exact witness at
\(r=s=282842415029/200000000000,z=0\) has positive value about
\(3.45423919000603\times10^{-7}\), which exceeds that smaller correction.
The strict angular perturbation \(z=10^{-10}\) is also positive by more
than the smaller correction. This is why the proposed integration uses the larger weight
\(4\times10^{-7}\). The exact graph theorem is independent of this
choice of weight.

The complete exact evaluator was run independently on the frozen base
candidate. It includes the direct rational Gram terms, all existing local
corrections, and the complete pair gauge, with no nonnegativity clamps.
The exploratory graph corrections at E weight \(2.4\times10^{-7}\)
made all 21 initial diagnostic points strictly negative. Their finite-point
conclusions remain valid, but the subsequent extra witness rejects that
smaller-weight candidate. The enlarged final screen uses
\(\delta_E=4\times10^{-7}\) and contains **242 exact point enclosures**.
Every corrected upper bound is strictly negative; the largest is displayed
approximately as \(-5.359490962379267\times10^{-10}\).

The following values use the larger weight. They are computed from the
retained exact base enclosures and the exact graph membership, including
the later rejecting point.

| Diagnostic configuration | Base cyclic value, displayed approximately | With final E weight \(4\times10^{-7}\) |
|---|---:|---:|
| Outer first exact discovery, \(r=s=25/16,z=3/64\) | \(+1.662849940636316\times10^{-7}\) | \(-2.337150059363684\times10^{-7}\) |
| Rational replay of the discovered outer maximum | \(+2.2013279690239963\times10^{-7}\) | \(-1.798672030976004\times10^{-7}\) |
| Middle antipodal discovery, \(z=0\) | \(+3.201572637488607\times10^{-8}\) | \(-3.679842736251139\times10^{-7}\) |
| Difficult middle/outer box midpoint | \(-1.508072251962956\times10^{-7}\) | \(-5.508072251962956\times10^{-7}\) |
| Later middle diagonal witness | \(+3.45423919000603\times10^{-7}\) | \(-5.45760809993970\times10^{-8}\) |

The exact endpoint replay also includes points just outside the widened
band. There E contributes zero, and the tested values near radius 1.57
remain about \(4.7\times10^{-8}\) below zero. These checks support the
choice of the widened band; they do not prove all of its complement is
negative.

For a global bound, all 112 continuous branches of the final frozen
candidate must still be accepted: 84 full branches and 28 exterior
branches, including all one-sided radial closures and angular boundaries.
The completed local theorem, the exact objective computation, and the
finite point screen are separate from that remaining requirement.

## Reproducible sources

The exact graph proof is self-contained apart from the shared standard-
library polynomial utilities in `verify_exterior_graphs.py`.

- `antipodal_wide_exterior_graph_certificate.json`: complete rational graph
  definition and spherical count polynomial.
- `verify_antipodal_wide_exterior_graph.py`: fresh count replay, exact
  projected-degree inequality, and selected-distance proof.
- `antipodal_wide_exterior_graph_verification.json`: successful exact proof
  report, including all constants and sign-proof statistics.
- `wide_ab_wide_e_graph_repair_bounds.py`: canonical A/B/E interval helper.
- `audit_antipodal_wide_graph.py` and
  `antipodal_wide_graph_preparation_audit.json`: independent fresh A/B/D/E
  proof replay, import binding, endpoint controls, and exact charge audit.
- `probe_antipodal_wide_graph.py` and
  `antipodal_wide_graph_point_replay.summary.json`: exact discovery replay
  and angular boundary diagnostics, with the full rational enclosures in
  the corresponding `.exact.json` file.

The updated final point screen is recorded in
`current/four_graphs_e4_peak_screen.json`; the smaller-weight rejection is
recorded in `current/robust_diagonal_antipodal_points.json`.

All preceding frozen graph, candidate, and archive sources were left
unchanged. The E polynomial and geometric proof do not depend on the weight;
the larger weight is supplied explicitly by the final coordinator.
