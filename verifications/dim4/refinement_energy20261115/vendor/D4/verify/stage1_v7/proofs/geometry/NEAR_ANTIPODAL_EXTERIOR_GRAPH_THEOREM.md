# An exterior graph on near neighbors: exact capacity 120

## Scope and current status

The following local theorem is proved independently of every candidate coefficient. Its application with weight

\[
\delta=\frac1{5000000}=2\times10^{-7}
\]

adds exactly \(3\times10^{-6}\) to the density objective. It repairs the six exact closure and interior witnesses in `round9_new_small_radius_outside55_exact.json`. A global density conclusion still requires the complete source verification of the candidate with this added energy and the source-bound origin-credit composition. The uncorrected round-9 candidate is rejected by those strict interior witnesses.

## The local graph

Fix a center of a packing with minimum distance one. Let \(C\) consist of its neighbors whose radii belong to the closed interval \([1,103/100]\). Join two distinct vertices if their directional cosine belongs to the closed interval \([-1,-3/4]\). Let \(e(C)\) count each selected edge in both orders. Then

\[
\boxed{e(C)\le120.}
\]

The origin is not a vertex. All radial and angular endpoints are included.

### At most 30 eligible vertices

For distinct vertices \(ru,sv\), separation implies

\[
u\cdot v\le\frac{r^2+s^2-1}{2rs}\le\rho:=\frac{5609}{10609}.
\]

The derivative with respect to \(r\) has numerator \(r^2-s^2+1\), which is at least \(9391/10000>0\) on the whole radial square; the same holds with the variables exchanged. Thus the maximum occurs at \(r=s=103/100\). The four corner values are \(1/2,103/200,103/200,5609/10609\).

Use the polynomial

\[
p(t)=(t-\rho)(t+1)(t+1/2)^2t^2.
\]

Its factorization proves \(p(t)\le0\) on the entire closed interval \([-1,\rho]\). In the normalized spherical zonals \(Z_k(t)=U_k(t)/(k+1)\) on \(S^3\), its coefficients are

\[
\left(
\frac{46955}{678976},\frac{43173}{169744},
\frac{73173}{169744},\frac{5625}{10609},
\frac{153045}{339488},\frac{46827}{169744},\frac7{64}
\right).
\]

Every coefficient is positive. The verifier reconstructs the Chebyshev recurrence and checks the polynomial identity coefficient by coefficient. The spherical addition theorem gives, for every possible number \(N\) of eligible vertices,

\[
p_0N^2\le\sum_{i,j}p(u_i\cdot u_j)\le Np(1).
\]

For \(N>0\),

\[
N\le\frac{p(1)}{p_0}=\frac{288000}{9391}<31,
\qquad N\le30.
\]

The empty set also satisfies the claim. No previously assumed kissing or neighbor count is used.

### At most four selected neighbors per vertex

Fix a graph vertex with direction \(u\). For each selected neighbor \(v_i\), write \(h_i=-u\cdot v_i\ge3/4\) and let \(w_i=v_i+h_i u\) be its projection into \(u^\perp\cong\mathbb R^3\). For two such neighbors,

\[
w_i\cdot w_j=v_i\cdot v_j-h_i h_j
\le\rho-\frac9{16}
=-\frac{5737}{169744}<0.
\]

Five vectors in \(\mathbb R^3\) admit an affine dependence. Separating its positive and negative coefficients produces two equal, nonempty positive combinations with common value \(w\). Their cross inner product is strictly negative, because all cross pairs are distinct; it is also \(\lVert w\rVert^2\ge0\), a contradiction. Consequently every graph degree is at most four, and \(e(C)\le30\cdot4=120\).

## Local energy, support and exact charge

For any \(\delta\ge0\), add the energy

\[
E_\delta(C)=\delta(120-e(C))\ge0.
\]

Its origin constant is \(120\delta\), its origin coupling and nonorigin diagonal are zero, and its off-diagonal kernel is \(-\delta\) on a selected ordered pair and zero elsewhere. The count is already ordered, so no second factor of two is introduced in the origin charge.

Every selected pair satisfies

\[
|ru-sv|^2\ge r^2+s^2+\frac32rs\ge\frac72
>\left(\frac74\right)^2.
\]

The selected edge therefore exceeds the kernel support radius. The new correction is identically zero on every full triangle, and on an exterior triangle it subtracts \(\delta\) when the two supported radii and their directional angle meet the graph conditions. The pair expression and the original complete-pair gauge are unchanged. The source-bound origin credit is also unchanged, because it only shifts the original degree-zero Gram component and restores its physical coupling.

In the normalization \(\widehat f(0)=\pi^2/4\), the density charge is

\[
\frac{120\delta}{8}=15\delta=\frac3{1000000}
\quad\text{for }\delta=\frac1{5000000}.
\]

## Closed-box evaluator

An exterior branch uses \(1\le r,s\le7/4\), \(0\le z\le1\), with

\[
t=-1+z\frac{(r+s)^2-49/16}{2rs}.
\]

This expression increases in each coordinate on the full parameter box. Its \(r\)-derivative has the sign of \(z(r^2-s^2+49/16)\ge0\); the \(s\)-derivative is analogous. The helper therefore encloses the angle by the two corresponding corner values. It returns the lower correction \(\delta\) only when the entire closed box is selected, and zero otherwise. Its optional upper correction is zero only when selection is impossible. At a singleton, both bounds coincide with the exact closed indicator. Boxes straddling an endpoint retain the interval \([0,\delta]\).

The certificate, verifier and helper are `near_antipodal_exterior_graph_certificate.json`, `verify_near_antipodal_exterior_graph.py`, and `near_antipodal_exterior_graph_bounds.py`. The fresh local proof is `near_antipodal_exterior_graph_verification.json`; exact endpoint, straddle and repaired-point checks are in `near_antipodal_graph_endpoint_and_point_audit.json`. These local artifacts do not assert completion of the 112 global cyclic branches.
