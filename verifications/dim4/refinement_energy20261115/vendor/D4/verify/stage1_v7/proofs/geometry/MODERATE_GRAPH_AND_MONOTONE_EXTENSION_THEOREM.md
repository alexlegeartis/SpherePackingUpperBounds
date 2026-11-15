# A moderate-angle graph energy and exact branch transfer

## Status and purpose

The graph theorem below and its application to `current/candidates/round9_six_graphs_degree18` are **completed**. The fresh exterior `(5,5)` proof certified its whole closed domain in 69,385 nodes, and the exact source-bound monotone composer accepted all 112 branches. Its final record is `audit/round9_six_graphs_monotone_origin_credit_global_verified.json`, with status `EXACT_GLOBAL_MONOTONE_GRAPH_EXTENSION_AND_ORIGIN_CREDIT_VERIFIED`. The preceding V6 sweep remains a rejected overall run: only its 111 individually certified branches are transferred by the theorem and exact source checks below, and its failed packet is excluded.

The added graph has ordered-edge capacity **1880**, covers radii `[7/5,3/2]` and directional cosines `[-4/5,-11/20]`, and costs exactly `47/20000000 = 0.00000235` in the density objective at weight `1/100000000`. Together with the retained graph charges and exact origin credit, the verified strict upper bound for arbitrary four-dimensional sphere packings is **Delta_4 < 78749932299/125000000000 = 0.629999458392**. The exact objective is approximately `0.6299994583912497039091`, strictly below the stated rational threshold; the decimal is not being used to decide the inequality.

## 1. An unrestricted spherical edge bound

Consider any finite unit-separated set of displacement vectors in four dimensions whose radii lie in

\[
I=[7/5,3/2].
\]

Write their unit directions as \(u_1,\ldots,u_N\). Distinct directions satisfy

\[
u_i\cdot u_j
\le \frac{r_i^2+r_j^2-1}{2r_ir_j}
\le\frac79=\rho.
\]

For fixed \(s\), the derivative numerator of the middle expression with respect to \(r\) is \(r^2-s^2+1\), an increasing function. Every interior critical point is therefore a minimum. Applying the endpoint argument in both variables proves that the radial square's maximum occurs at a corner; its largest corner value is \(7/9\).

Select an unordered pair precisely when

\[
-\frac45\le u_i\cdot u_j\le-\frac{11}{20},
\]

including both endpoints, and let \(e\) be twice the number of selected unordered pairs. Thus \(e\) counts **ordered distinct pairs** and is an even nonnegative integer.

Let \(Z_k(t)=U_k(t)/(k+1)\) be the normalized zonal polynomial on \(S^3\), where

\[
U_0=1,\qquad U_1=2t,\qquad U_{k+1}=2tU_k-U_{k-1}.
\]

The frozen rational degree-16 polynomial is

\[
P(t)=\sum_{k=0}^{16}a_kZ_k(t).
\]

Its complete coefficients and independently reconstructed power coefficients are in `geometry/moderate_exterior_graph_certificate.json`. Every \(a_k\) is nonnegative, and

\[
a_0=\frac{109916026941}{500000000000}>0,
\qquad
P(1)=\frac{40670769968163}{1000000000000}.
\]

The exact checker proves

\[
P(t)\le0\quad(-1\le t\le7/9),
\qquad
P(t)\le-1\quad(-4/5\le t\le-11/20).
\]

These are whole-interval assertions, including every endpoint. The checker reconstructs the zonal expansion coefficient by coefficient and proves the three closed interval pieces by 69 exact rational Bernstein-subdivision nodes, with maximum depth seven. A separate audit proves the same strict signs through primitive pseudo-remainder Sturm sequences: each relevant polynomial has zero roots on its interval and positive signs at both endpoints.

The spherical addition formula gives

\[
a_0N^2\le\sum_{i,j=1}^N P(u_i\cdot u_j)
\le NP(1)-e.
\]

The right inequality retains the diagonal \(NP(1)\), charges at least one for each selected ordered pair, and uses nonpositivity on all other distinct pairs. Consequently

\[
e\le NP(1)-a_0N^2.
\]

This proves an edge bound without assuming any separate bound on the number of vertices. The concave quadratic's real center is

\[
\frac{P(1)}{2a_0}
=\frac{13556923322721}{146554702588},
\]

which lies between 92 and 93. Its maximum over all nonnegative integer cardinalities is attained at 93 and equals

\[
\frac{1881054173013741}{1000000000000}
=1881.054173013741<1882.
\]

Since \(e\) is even, **\(e\le1880\)**. The independent audit also enumerates all integer cardinalities from zero through 400, checks the same maximum, and verifies that the quadratic is negative and decreasing beyond 400. No inherited spherical count, latitude bound or cap-projection monotonicity assumption enters this argument.

## 2. The exact local energy and its origin cost

For any finite unit-separated packing neighborhood \(C\), apply the preceding graph only to its vertices with radii in \(I\). Let \(e_I(C)\) be its selected ordered-edge count. For every exact weight \(\delta\ge0\),

\[
E_\delta(C)=\delta(1880-e_I(C))\ge0.
\]

Its origin constant is \(1880\delta\), its origin-to-neighbor coupling is zero, and its nonorigin self-diagonal is zero. Its kernel on a selected distinct pair is \(-\delta\), and it is zero on every other nonorigin pair. The empty neighborhood also satisfies the inequality.

The selected cosine is negative, so squared distance is increasing in both positive radii. Every selected edge therefore obeys

\[
|x-y|^2\ge2(7/5)^2\left(1+\frac{11}{20}\right)
=\frac{1519}{250}
>\frac{49}{16}=R^2,
\qquad R=7/4.
\]

Thus a selected edge is strictly exterior to the supported full-triangle domain. On any triangle with all three sides at most \(R\), the graph contributes zero at each anchor. On an exterior triangle with exactly one side greater than \(R\), at most the opposite anchor can contribute; at either other anchor one incident displacement has radius greater than \(R\), outside \(I\). The cyclic decrement is therefore either zero or \(\delta\), and is always nonnegative.

Both radial endpoints and both angular endpoints are included in the graph. The separate helper computes a conservative exact lower bound on that decrement for every closed coordinate box. Boxes straddling a threshold may have lower bound zero. A uniform upper bound is \(\delta\) on exterior boxes; at a singleton the exact closed-membership rule gives the exact value for both bounds. All full-branch bounds are zero. The minimum-distance gap also shows that no selected edge is lost or ambiguously added at the closed full/exterior interface \(c=R\).

The graph changes neither the physical pair function

\[
q(x)=f(x)+2h(x)+K(x,x)
\]

nor the already proved complete-pair gauge. It has no linear or diagonal term with which those quantities could interact. The added origin charge must nevertheless be included in the objective. With the retained normalization \(\alpha=2\), \(G(0)=1\), its density cost is

\[
\frac{1880\delta}{8}.
\]

At \(\delta=10^{-8}\), this is exactly \(47/20000000\). The previous three active graphs remain present with unchanged weights, so the total graph density cost becomes

\[
\frac{4299}{200000000}+\frac{47}{20000000}
=\frac{4769}{200000000}.
\]

The existing origin-credit identity still applies to the same degree-zero Gram block, linear coupling and radial coefficients. Its local range equations and radial nonnegativity must be replayed against the same frozen scalar, origin and kernel bytes. Applying that credit and adding this graph commute: each graph origin charge is paid once, the credit is taken once, and the physical pair and cyclic functions have precisely the stated changes.

## 3. Monotone transfer of completed branch proofs

**Lemma.** Let an old source have a valid local-energy decomposition, scalar and complete pair proofs, and an exact cyclic proof on a closed branch \(\Omega\). Add an independently nonnegative local energy whose nonorigin kernel has nonpositive cyclic contribution on \(\Omega\), with its origin charge included. If every old physical coefficient and every other local-energy weight is unchanged, the old cyclic proof remains valid for the new source on \(\Omega\).

**Proof.** At every configuration \(T\in\Omega\), write \(D(T)\ge0\) for the added decrement. The exact physical identity is

\[
\Gamma_{\mathrm{new}}(T)
=\Gamma_{\mathrm{old}}(T)-D(T).
\]

The old branch proof supplies \(\Gamma_{\mathrm{old}}(T)\le0\) for every \(T\in\Omega\). Therefore \(\Gamma_{\mathrm{new}}(T)\le0\) everywhere on that same closed branch. No sampling or extrapolation is involved. The new local energy is independently nonnegative, and its zero coupling and self-diagonal preserve the old pair proof. Including its origin constant gives the correct new objective. \(\square\)

The old source need not have had a successful **global** run for this branchwise implication to hold. A failed branch does not revoke other complete, source-bound branch proofs. It also cannot be inherited: that branch must receive a new exact proof or another separately justified argument.

The completed coverage is the old V6 source's **111 certified closed branches** together with the fresh proof of its remaining exterior `(5,5)` branch after the moderate correction. The normal V7 coordinator continues to define all 112 branches. A separate selected-branch run is partial and must never emit a global certificate by itself. The accepted final composition independently verifies:

- The old certified packets' mathematical source bindings and frozen input hashes, including their original graph weights and the complete pair/dependency records.
- Equality of the old scalar, origin, kernel and all old local-energy inputs in the new candidate; the only added physical term is the separately proved moderate graph, with its full origin charge.
- The fresh selected-branch proof's source bindings, exact identities, complete closed domain, graph weight and `CERTIFIED` status.
- An exact partition of the declared 84 full and 28 exterior branch identities into the 111 retained branches and the one newly proved branch, with no missing or duplicate labels.
- The complete objective identity, and the separately replayed source-bound origin credit if the credited threshold is claimed.

All these proof gates passed for the frozen V7 source. The exact coverage record contains 84 full and 28 exterior branches, 672 supplementary exact arithmetic crosschecks, and 560 signed spherical checks. The complete sign proofs use 391,420 accepted subdivision nodes and 195,766 leaves in aggregate. The finite-neighborhood summation, periodic limit and cube-periodization argument therefore applies unchanged to arbitrary upper-density packings. A point screen, a resource-limited `UNKNOWN`, a failed packet, a selected-stage status alone or an unbound source record would not justify this conclusion.

## 4. Frozen mathematical evidence

- `geometry/moderate_exterior_graph_certificate.json`: all exact harmonic and power coefficients.
- `geometry/verify_moderate_exterior_graph.py`: standalone standard-library graph proof.
- `geometry/moderate_exterior_graph_verification.json`: complete interval proofs and unrestricted integer-cardinality maximum.
- `geometry/moderate_exterior_graph_bounds.py`: closed-membership bounds and exact density charge.
- `geometry/SIX_GRAPHS_SMALL_CAP_BOUNDS_HELPER_CONSTRUCTION.json`: five edits whose inverse reconstructs the frozen V6 helper byte for byte.
- `geometry/six_graphs_small_cap_bounds_helper_exact_audit.json`: independent Sturm/cardinality checks and 1,406 singleton plus 270 rational box controls of the actual additive helper code.
- `geometry/ROUND9_LOCAL_ENERGY_TRANSFER_AUDIT.md`: the unchanged finite, periodic, arbitrary-density and complete-pair-gauge arguments; its current status is the completed V7 application, with the earlier V6 rejection retained explicitly.

The exact global completion is recorded in `audit/round9_six_graphs_monotone_origin_credit_global_verified.json`; `current/FINAL_MONOTONE_CERTIFICATE_SUMMARY.json` independently checks its current hashes and reports the aggregate proof statistics. The selected-branch report remains correctly labelled partial on its own; only the completed source-bound composition gives the global conclusion.
