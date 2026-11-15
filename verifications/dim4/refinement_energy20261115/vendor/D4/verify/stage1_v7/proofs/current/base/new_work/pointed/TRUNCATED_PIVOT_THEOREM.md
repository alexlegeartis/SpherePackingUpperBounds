# A geometric truncation that weakly dominates anchored pointed rays

This is an analytic statement about admissible local energies. It improves
the search cone but does not by itself give a numerical packing bound.
Use the global local-energy theorem in `../LOCAL_ENERGY_THEOREMS.md` and the
pointed construction in `THEOREM.md`. Packing-center separation is one.

## 1. When a pointed ball cannot contain a distinct neighbor

Let the distinguished packing center be 0, let another center be
\(y=s\widehat y\), and use the deterministic pointing direction
\(v=\widehat y\), so \(\sigma=1\). Assume

\[
0<d\le a,\qquad a^2<\frac12.
\]

The probe ball is \(B(d\widehat y,a)\), and it contains the origin.
If \(x\ne0,y\) is a packing center, put \(r=|x|\). Separation from
\(y\) gives

\[
x\cdot\widehat y\le\frac{r^2+s^2-1}{2s}.
\]

Since \(s\ge1>d\) and \(r\ge1\),

\[
\begin{aligned}
|x-d\widehat y|^2
&\ge \left(1-\frac ds\right)r^2+d^2-ds+\frac ds\\
&\ge 1+d^2-ds.
\end{aligned}
\]

Consequently the ball contains no such distinct center whenever

\[
s<T(d,a):=\frac{1+d^2-a^2}{d}.
\]

Moreover \(T(d,a)=d+(1-a^2)/d\) decreases on \(0<d\le a\), because
\(a^2<1/2\). Its minimum there is \(1/a>\sqrt2\). The strict inequality
at the onset matters for a closed probe ball.

## 2. Exact weak dominance

Choose a rational cutoff \(C<T(d,a)\), and replace a nonnegative pivot
weight \(w(s)\) by

\[
\widetilde w(s)=w(s)\,\mathbf1_{s>C}.
\]

This is still an admissible nonnegative weight in the pointed ball-count
energy. On every actual pair of distinct neighbors, removing the weight
for \(s\le C\) removes a term whose probability was already zero.
Therefore the **nonorigin off-diagonal kernel, and hence every full and
exterior cyclic expression, is unchanged** on its feasible domain.

The pair contribution of the old ray is

\[
\beta\,[M-P(x\mid x)]w(|x|).
\]

For the current rays \(M\ge1\), so this expression is nonnegative.
Truncating the weight sets it to zero below the cutoff and leaves it
unchanged above. The pointed origin constant is zero in either case.
Thus replacing any such ray by its truncation weakly improves every pair
constraint, preserves every cyclic constraint, and leaves the objective
unchanged. Reoptimization can use the resulting slack; an improved numerical
objective still requires complete global verification.

No extra PSD feature is needed for this dominance operation in a joint
optimization: the pointed linear coupling is a separate local-energy
component. A new PSD feature would be an optional enlargement of the
search space. A post-hoc origin shift that is required to preserve an old
pair expression would have an additional range requirement.

## 3. A concrete window disjoint from nearby D4 shells

The rational choices

\[
d=\frac{707106}{10^6},\quad
a^2=\frac{499999}{10^6},\quad
C=\frac{707107}{500000}=1.414214
\]

satisfy \(d^2\le a^2<1/2\), \(C^2>2\), and \(a^2C^2<1\).
Thus \(\sqrt2<C<1/a\le T(d,a)\). The origin-containing ball has at
most four nonorigin neighbors by the elementary dimension-four bound.

Use, for example,

\[
\widetilde w(s)=\mathbf1_{C<s\le3/2},
\]

or the quartic inner-annulus Bernstein feature multiplied by the same
indicator. Normalized D4 has possible nonzero squared lengths \(1,2,3\)
within the global radius \(7/4\); none falls in this window. Hence the
whole pointed local energy is exactly zero on that neighborhood. This is
a property of the additional inequality, not an assumption that a packing
being bounded has D4 local structure.

Finite-dual pricing suggested that this window could be useful. Those
prices are search diagnostics only and are recorded separately in
`../copositive/pointed_radial_window_prices.json`.

## 4. Verification without another cyclic radial partition

Store the untruncated nonnegative Bernstein coefficient vector together
with `pivot_lower_cutoff=C`. The validator checks the preceding parameter
conditions and the strict rational inequality \(C<T(d,a)\). Because the
off-diagonal kernels agree on every feasible triangle, the cyclic verifier
may use the untruncated expression and its existing radial partition.

For the pair calculation, define the auxiliary restored origin row using
the untruncated linear contribution,

\[
v_{\rm aux}=v_{\rm PSD}+\frac12\sum_j\beta_jM_jw_j.
\]

For a truncated ray subtract the following nonnegative amount from its
restored pair expression:

\[
J_j(r)=\beta_jw_j(r)
\begin{cases}
M_j,&r\le C_j,\\
P_j(x\mid x),&r>C_j.
\end{cases}
\]

Then the remaining contribution is exactly
\(\beta_j(M_j-P_j(x\mid x))w_j(r)\mathbf1_{r>C_j}\).
Untruncated rays continue to use \(J_j=\beta_jw_jP_j(x\mid x)\).
The exact pair verifier encloses this nonnegative reduction on each
interval and keeps the endpoint convention explicit. The auxiliary
origin row is an algebraic representation of the pair expression; the
actual PSD row still separately satisfies its checked range identity.
