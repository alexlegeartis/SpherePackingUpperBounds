# A broader antipodal graph with 1,358 ordered edges

## Proven local theorem

Let \(C\subset\mathbb R^4\) be a finite set of neighbors of a packing
center, with pairwise distances at least one. Retain the points whose
radii about that center lie in

\[
I=\left[\frac{33}{28},\frac74\right].
\]

Join two retained points when their directional cosine is at most
\(-7/8\). The resulting graph has at most **1,358 ordered edges**. Every
selected edge has length greater than \(7/4\). Thus, for any
\(\delta\ge0\),

\[
\mathcal E(C)=\delta\bigl(1358-e(C)\bigr)\ge0
\]

is a valid exact local energy with zero linear term and zero kernel
diagonal. Its entire density cost in the current normalization is
\(1358\delta/8\).

This geometric theorem is independent of the numerical kernel. The
single-graph proposal at \(\delta=10^{-7}\) was subsequently rejected by
a larger positive antipodal configuration. A two-graph refinement can
combine this theorem with the earlier 604-edge theorem. No proposed
weights in this note are a global density certificate.

**Later integration status.** The frozen dual graph route has passed an
independent exact integration audit with broad weight
\(b=6\times10^{-8}\) and narrow weight \(n=2.1\times10^{-7}\). Its
exact density charge is \(651/25000000=0.00002604\). Both local graph
proofs remain valid, but the unmodified round-six scalar function failed
its forbidden-halfline gate. A separately proved slower Gaussian local
capacity is available for a scalar repair; that repair and the complete
global verification are separate requirements. See
`dual_antipodal_integration_audit.json` and
`GAUSSIAN_THREE_QUARTERS_LOCAL_CAPACITY.md` for the current local results.

## 1. Fresh annular count: at most 151 vertices

For \(x=ru\), \(y=sv\), where \(u,v\) are unit vectors, packing separation
implies

\[
u\cdot v\le\frac{r^2+s^2-1}{2rs}.
\]

For a fixed positive \(s\), the derivative of the right side in \(r\)
has numerator \(r^2-s^2+1\); an interior critical point is a minimum.
Therefore the maximum on the radial rectangle occurs at a corner. For
\(r,s\in I\), the four corner values are

\[
\frac{697}{1089},\quad\frac{41}{49},\quad
\frac{41}{49},\quad\frac{41}{49}.
\]

All distinct retained directions consequently have inner product at most
\(\rho=41/49\). The supplied rational degree-20 polynomial in the
normalized \(S^3\) zonal basis has nonnegative nonconstant coefficients,
constant coefficient \(a_0=999999/1000000>0\), and is nonpositive on the
whole interval \([-1,41/49]\). Its complete rational Bernstein proof
uses 93 subdivision nodes, 47 leaves, and maximum depth 12.

The addition formula gives

\[
|C_I|\le\frac{p(1)}{a_0}
=\frac{151196972764393}{999999000000}<152.
\]

Hence \(|C_I|\le151\). The complete polynomial is included in the new
graph certificate and its sign is proved afresh by the new checker. No
retained numerical or previously asserted neighbor count is assumed.

## 2. Normalize the perpendicular projections

Fix a graph vertex with direction \(u\), write \(p=-u\), and put
\(\beta=7/8\). Its selected neighbor directions \(v_i\) satisfy

\[
h_i=v_i\cdot p\ge\beta.
\]

If \(h_i=1\) for one neighbor, then \(v_i=p\). Another selected neighbor
would have \(v_i\cdot v_j=h_j\ge\beta>\rho\), contradicting directional
separation. That case has degree at most one.

Otherwise the unit perpendicular projections

\[
w_i=\frac{v_i-h_i p}{\sqrt{1-h_i^2}}
\]

lie on a two-sphere in \(p^\perp\cong\mathbb R^3\). Their mutual inner
products satisfy

\[
w_i\cdot w_j\le
F(h_i,h_j):=
\frac{\rho-h_i h_j}{\sqrt{(1-h_i^2)(1-h_j^2)}}.
\]

The partial derivative is

\[
\frac{\partial F}{\partial h_i}
=\frac{\rho h_i-h_j}
{(1-h_i^2)^{3/2}\sqrt{1-h_j^2}}<0,
\]

since \(\rho h_i-h_j\le\rho-\beta<0\). The same holds in the other
variable, so the maximum occurs at \(h_i=h_j=\beta\). Therefore

\[
w_i\cdot w_j\le
\kappa:=\frac{\rho-\beta^2}{1-\beta^2}
=\frac{223}{735}.
\]

The exact monotonicity margin is \(\beta-\rho=15/392>0\). The preceding
derivative computation justifies the projected cutoff for the entire
latitude interval, not just for neighbors on its boundary.

## 3. A factored quartic proves degree at most nine

Consider

\[
q(t)=(t+1)\left(t-\frac{223}{735}\right)
             \left(t+\frac23\right)^2.
\]

This factorization immediately proves \(q(t)\le0\) throughout
\([-1,223/735]\). In Legendre polynomials normalized by \(P_\ell(1)=1\),
its exact expansion is

\[
q(t)=\frac{62}{147}P_0(t)
+\frac{37144}{33075}P_1(t)
+\frac{1214}{945}P_2(t)
+\frac{2984}{3675}P_3(t)
+\frac{8}{35}P_4(t).
\]

Every coefficient is positive. The exact checker reconstructs the
Legendre polynomials from their recurrence and verifies this identity
coefficient by coefficient.

For \(m\) projected neighbors, the \(S^2\) addition formula gives

\[
\frac{62}{147}m^2
\le\sum_{i,j}q(w_i\cdot w_j)
\le m\,q(1).
\]

Since

\[
q(1)=\frac{5120}{1323},\qquad
\frac{q(1)}{62/147}=\frac{2560}{279}<10,
\]

integrality implies **\(m\le9\)**. This count requires no numerical
optimization, approximate roots, or additional interval subdivision.

Summing graph degrees now gives \(e(C)\le151\cdot9=1359\). Each ordinary
graph edge contributes two ordered edges, so \(e(C)\) is even and

\[
\boxed{e(C)\le1358.}
\]

## 4. Exterior support and exact local-energy cost

For a selected pair,

\[
|x-y|^2\ge
2\left(\frac{33}{28}\right)^2\left(1+\frac78\right)
=\frac{16335}{3136}.
\]

Its exact gap above the squared support radius is

\[
\frac{16335}{3136}-\left(\frac74\right)^2
=\frac{6731}{3136}>0.
\]

The graph kernel equals \(-\delta\) on selected distinct ordered pairs
and zero elsewhere. Its diagonal and linear term vanish. Adding the
origin constant \(1358\delta\) gives the nonnegative local energy.
Neither the complete pair inequality nor its existing pair gauge changes.

All full-triangle graph decrements vanish. On an exterior box the helper
subtracts a graph weight only if both complete radial intervals and the
complete computed cosine interval lie inside the graph support. A
straddling box receives zero lower credit until subdivision resolves it.
At an exact point, the membership predicate is exact.

At \(\delta=10^{-7}\), the density charge is exactly

\[
\frac{1358\delta}{8}
=\frac{679}{40000000}=0.000016975.
\]

## 5. Rejection diagnostics and two-graph allocation

The earlier fallback ended with 110 certified branches and two actual
`VIOLATION` packets. One has strict interior radii
\(r=s=8983761/6400000\), \(z=51/64\), and directional cosine near
\(-0.02552\). It is already negative on the stronger bare round-six
candidate. The other survives at

\[
r=\frac{1499999999}{10^9},\quad s=\frac{103}{64},\quad
z=\frac{19}{256},
\]

with cosine approximately \(-0.8984558698322078\). Its complete bare
round-six value is about \(2.0282658732479322\times10^{-8}\), with exact
rational enclosure

\[
\frac{10141}{500000000000}
\le\Gamma\le
\frac{20283}{1000000000000}.
\]

The new graph at weight \(10^{-7}\) makes that point strictly negative.
However, a later strict interior antipodal point near
\((r,s,z)=(1.4187922484720006,1.546674626372355,0.01787352828670068)\)
has complete bare value about \(2.53545328798179\times10^{-7}\).
Its cosine is about \(-0.9766584144\), so a single weight of \(10^{-7}\)
is insufficient. That application is rejected; the geometric theorem
remains valid.

The new larger peak lies inside both this broad graph and the earlier
604-edge graph, whose angular interval is \([-1,-183/200]\). One can
therefore use two independent weights: \(b\) on the broad graph and
\(n\) on the narrow graph. The exact density charge is

\[
\frac{1358b+604n}{8}.
\]

Both weights apply inside the narrower band; only \(b\) applies on the
additional broad band. Independent nonnegative local energies justify
their sum without changing either theorem.

The reoptimized broad-only point at
\((r,s,z)=(1.5,1.607710943647149,0.07405172899659528)\) has exact bare
value approximately \(4.80235207066555\times10^{-8}\), and remains
positive after moving the first radius inward by \(10^{-9}\). Thus
\(b=4\times10^{-8}\) is rejected. At \(b=5\times10^{-8}\), this tested
point has negative margin only about \(1.98\times10^{-9}\).

For illustration, \(b=5\times10^{-8}\), \(n=2.25\times10^{-7}\) costs
exactly \(0.000025475\), and the sum \(b+n=2.75\times10^{-7}\) exceeds
the discovered inner peak. Adding that cost to the stored round-six
objective gives approximately \(0.6299998737751055\). These are
conditional budget calculations and finite point checks. The maximum on
every continuous branch still requires the full exact verifier.

The subsequently selected stronger broad weight
\(b=6\times10^{-8}\), with \(n=2.1\times10^{-7}\), gives total
decrement \(2.7\times10^{-7}\) in the narrow band. A distinct narrow
basin at
\((1.5,1.5877395337381568,0.05031723440882012)\) has complete exact
mixed round-six value about \(2.461757772931045\times10^{-7}\).
After this graph correction its value is about
\(-2.38242227068955\times10^{-8}\). Its strict inward radial
perturbation is already negative. These checks are recorded in
`round6_second_narrow56_exact_summary.json`; they do not replace a whole
branch proof.

## Exact files and evidence

- `antipodal_cap9_exterior_graph_certificate.json`: full graph definition,
  fresh annular-count polynomial, and five exact Legendre coefficients.
- `verify_antipodal_cap9_exterior_graph.py` and
  `antipodal_cap9_exterior_graph_verification.json`: complete independent
  proof and successful report.
- `wide_ab_cap9_e_graph_repair_bounds.py`: canonical interval API, with
  capacities \((862,1102,1358)\) and exploratory default weights
  \((0,0,10^{-7})\); final applications supply their weights explicitly.
- `cap9_graph_preparation_audit.json`: fresh proof replay, source binding,
  20 endpoint controls, and the initial exact witness repair.
- `fallback_failures_vs_round6_exact_summary.json`: exact comparison of
  both failed fallback packets on the fallback and bare round-six data.
- `round6_newest_antipodal56_summary.json` and
  `round6_broad_only56_peak_summary.json`: later exact failures that
  invalidate the smaller exploratory weights.
- `dual_antipodal_integration_audit.json`: exact reconstruction of the
  frozen dual coordinator and helper, fresh local graph proofs, all
  source bindings, origin charging, overlap, and 28 endpoint/straddling
  controls.
- `round6_second_narrow56_exact_summary.json`: direct exact evaluation
  of the later second narrow basin on both native and mixed candidates.

No earlier mathematical source or candidate input was modified.
