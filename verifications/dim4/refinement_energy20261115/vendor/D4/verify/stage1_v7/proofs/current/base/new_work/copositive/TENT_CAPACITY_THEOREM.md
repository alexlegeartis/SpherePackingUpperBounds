# Pointed tent occupancy energies

This note proves a new family of packing-local energies. It does not assert
that a numerical optimization using the family is a global density bound.
The accepted certificate at `pointed/nested_repair2` is unchanged.

## 1. A dimension-independent weighted occupancy inequality

Let

\[
\phi(z)=(1-2\lVert z\rVert^2)_+.
\]

For every unit-separated set \(X\subset\mathbb R^m\) and every center
\(c\),

\[
\sum_{x\in X}\phi(x-c)\le1.
\tag{1}
\]

Only finitely many points contribute. Write them as \(x_1,\ldots,x_n\)
and put \(z_i=x_i-c\). The pair-distance identity gives

\[
n\sum_i\lVert z_i\rVert^2-\Bigl\lVert\sum_i z_i\Bigr\rVert^2
=\sum_{i<j}\lVert z_i-z_j\rVert^2
\ge\frac{n(n-1)}2.
\]

Therefore

\[
\sum_i\phi(z_i)
=n-2\sum_i\lVert z_i\rVert^2
\le1-\frac2n\Bigl\lVert\sum_i z_i\Bigr\rVert^2\le1.
\]

The case \(n=0\) is immediate. The support used here is the strict ball
\(\lVert z\rVert^2<1/2\). Points on its boundary contribute exactly zero.

If the packing contains the origin and \(c=dv\), \(\lVert v\rVert=1\),
deducting the origin's contribution gives the independent bound

\[
\sum_{x\in X\setminus\{0\}}\phi(x-dv)\le\min(1,2d^2).
\tag{2}
\]

## 2. Sharp anchored capacity in four dimensions

Suppose \(0<d\le1/\sqrt2\). Then the exact capacity for nonorigin
unit-separated points satisfying \(\lVert x\rVert\ge1\) is

\[
C_*(d)=\max\left(0,\max_{1\le n\le4}
\left[-n(1+2d^2)+4d\sqrt{\frac{n(n+1)}2}\right]\right).
\tag{3}
\]

First, at most four nonorigin points have positive tent weight. Relative
to the probe center, every such point has squared norm strictly below
\(1/2\), and the vector representing the packing origin has squared norm
\(d^2\le1/2\). The inner product of every pair of these vectors is strictly
negative, since their corresponding packing points are at distance at least
one. A strictly obtuse set in \(\mathbb R^4\) has at most five members.

For completeness, the latter fact has an elementary linear-algebra proof.
A linear dependence among vectors with strictly negative pairwise inner
products cannot have both positive and negative coefficients: equating the
two positive combinations and taking their inner product would make a
squared norm strictly negative. A nullspace of dimension at least two
would contain a dependence with both signs. Its dimension is consequently
at most one, and the number of vectors is at most the ambient dimension
plus one. This argument applies at the closed endpoint \(d^2=1/2\),
because every nonorigin contributing point lies in the strict support.

Now fix the number \(n\le4\) of contributing nonorigin points and set

\[
Q=\sum_i\lVert x_i\rVert^2\ge n,\qquad S=\sum_i x_i.
\]

The same pair-distance identity gives

\[
\lVert S\rVert^2\le nQ-\frac{n(n-1)}2,
\]

and hence

\[
\sum_i\lVert x_i-dv\rVert^2
\ge Q+nd^2-2d\sqrt{nQ-\frac{n(n-1)}2}.
\tag{4}
\]

The right side is increasing for \(Q\ge n\), because its derivative is
at least

\[
1-d\sqrt{\frac{2n}{n+1}}
\ge1-\sqrt{\frac n{n+1}}>0.
\]

Substituting \(Q=n\) into (4) proves the corresponding branch of (3).
Each positive branch is attainable: take \(n\) unit vectors with pairwise
inner product \(1/2\), and let \(v\) point along their sum. Their Gram
matrix is positive definite of rank \(n\le4\), so the configuration exists
in \(\mathbb R^4\). All its tent values are equal, and their sum is exactly
the displayed branch. Choosing a positive maximizing branch proves
sharpness. If every branch is nonpositive, the empty set attains zero.

This sharpness also survives adding a fixed negative-pole pivot
\(y=-sv\), \(s\ge1\): every simplex point has positive projection onto
\(v\), so its distance from \(y\) is greater than one. The pivot's own
tent value is zero. Thus merely conditioning on that existing pivot cannot
improve (3) for the negative-pole family.

### Exact rational capacity verification

A positive rational number \(C\) is a valid capacity if \(C\ge1\).
For \(d^2\le1/2\), it is also valid if \(C\ge2d^2\), or if all four
rational comparisons

\[
\left[C+n(1+2d^2)\right]^2\ge8d^2n(n+1),
\qquad n=1,2,3,4,
\tag{5}
\]

hold. The quantity inside the square is positive, so no sign ambiguity is
introduced by squaring. The checker uses (5), not floating-point square
roots. It need not establish that the supplied capacity is optimal.

For example, \(d=3/5\) admits

\[
C=\frac{44923461417477}{62500000000000},
\]

and the continuously priced distance \(d=59611859/10^8\) admits

\[
C=\frac{221437889840917466849}{312500000000000000000}.
\]

## 3. Packing-local kernel and its origin accounting

Fix the common annulus \(A=\{x:1\le\lVert x\rVert\le R\}\), a
nonnegative radial function \(w\) supported on \(A\), a positive capacity
\(C\) justified above, and a coefficient \(\beta\ge0\). For nonzero
arguments in \(A\), define

\[
P_d(r,t)=\left(1-2(r^2+d^2)-4drt\right)_+,
\]

\[
K(x,y)=-\frac\beta2
\left[w(\lVert y\rVert)P_d(\lVert x\rVert,\widehat x\cdot\widehat y)
+w(\lVert x\rVert)P_d(\lVert y\rVert,\widehat x\cdot\widehat y)\right].
\tag{6}
\]

Set \(K(x,y)=0\) when either nonzero argument lies outside \(A\). The
origin couplings and origin mass of this local-energy summand are

\[
h(x)=K(0,x)=\frac{\beta C}2w(\lVert x\rVert),\qquad K(0,0)=0.
\tag{7}
\]

For any finite packing neighborhood \(S\subset A\), its energy is exactly

\[
2\sum_{y\in S}h(y)+\sum_{x,y\in S}K(x,y)
=\beta\sum_{y\in S}w(\lVert y\rVert)
\left[C-\sum_{x\in S}\phi(x+d\widehat y)\right]\ge0.
\tag{8}
\]

Each inner sum is bounded by the capacity theorem with probe center
\(-d\widehat y\). All ordered pairs, including the diagonal, are present
in (8). The diagonal is zero because

\[
P_d(r,1)=(1-2(r+d)^2)_+=0\quad(r\ge1).
\]

Thus the pair contribution is exactly \(\beta Cw(r)\), and every
distinct-point kernel value is nonpositive. This is a packing-local
nonnegativity argument; the kernel (6) need not be positive semidefinite.
No lattice assumption enters it.

If the radial function belongs to an existing PSD feature span, the
coupling in (7) can be canceled by shifting its origin vector from \(v\)
to \(v-\beta Cw/2\). This keeps the combined pair expression unchanged.
For an invertible nonorigin Gram block \(G\), the new PSD origin mass is
\((v-\beta Cw/2)^T G^{-1}(v-\beta Cw/2)\). Every other local-energy origin
constant, including the constant term of any integer count polynomial,
must still be added. This construction is neither a claim of zero total
origin cost nor a replacement for the global scalar and cyclic checks.

An optimizer may normalize the column by writing \(\gamma=C\beta\).
The exact record always stores \(\beta\), so its pair contribution is
\(C\beta w\), and its negative cyclic contribution uses \(\beta P_d\).

## 4. Exact cyclic bounds

For the support radius \(R=7/4\) used in the current search, every
negative-pole tent is identically zero on a full triangle, whatever the
probe distance \(d>0\). Indeed, positive tent weight gives
\(\lVert x+d\widehat y\rVert^2<1/2\). Writing
\(q=x\cdot\widehat y\), we must have \(q<0\), and the perpendicular
squared norm is strictly below \(1/2\). Thus

\[
q<-\sqrt{\lVert x\rVert^2-1/2},\qquad
\lVert x-y\rVert^2>
\lVert x\rVert^2+\lVert y\rVert^2+
2\lVert y\rVert\sqrt{\lVert x\rVert^2-1/2}\ge2+\sqrt2.
\]

Since \((7/4)^2<2+\sqrt2\), the third side must be outside the support.
More generally all full corrections vanish when \(R^2\le2+\sqrt2\).
This exact observation also shows a limitation of the new columns: they
cannot directly repair a full-triangle violation. They can relieve the
exterior constraints during joint optimization.

On a full triangle with side lengths \((r,s,z)\), the uncut expression in
the directed tent is

\[
A(r,s,z)=1-2(r^2+d^2)-\frac{2d(r^2+s^2-z^2)}s.
\tag{9}
\]

It decreases in \(r\) and increases in \(z\). As a function of \(s\),
it is monotone or has one maximum at \(s=\sqrt{r^2-z^2}\). Consequently a
whole-box lower bound uses only rational corner evaluations; an upper
bound additionally includes that possible critical point, enclosed with
an integer square-root bound. Applying \(u\mapsto\max(0,u)\) preserves
the enclosure.

The exterior branch uses the normalized coordinate \(0\le z\le1\) and

\[
t=-1+z\frac{(r+s)^2-R^2}{2rs}.
\]

There the uncut tent is

\[
A(r,s,z)=1-2(r-d)^2-\frac{2dz((r+s)^2-R^2)}s.
\tag{10}
\]

On \(1\le r,s\le R\le2\), it decreases in \(s\) and \(z\), and is
concave in \(r\). Its possible interior maximum in \(r\) occurs at the
rational value

\[
r_* = \frac{d(1-z)}{1+dz/s}.
\]

The exact helper `pointed/interval_certificate_tent.py` uses these bounds,
nonnegative Bernstein bounds on \(w\), and the existing outward integer
Bernstein branch subdivision. It can add already-validated pointed-ball
corrections. Its output certifies only the requested cyclic branch; the
global verifier remains responsible for all origin, pair, scalar,
Fourier, count-capacity, and remaining branch conditions.

## 5. Numerical discovery and its limits

The source dual is the unverified `nested_window` numerical candidate with
objective approximately \(0.6311100367738846\). Pricing enumerated every
one-dimensional tent activation interval and each change in the sharp
capacity envelope. On each resulting interval, the reward divided by
capacity is a ratio of quadratics, whose stationary equation is quadratic.

The best negative-pole indicator price occurs near
\(d=0.596118591654659\), where the two-point and three-point capacity
branches cross. Its rationalized column price is approximately
\(-3.50767879821\). The universal-capacity family has its best indicator
distance near \(d=0.830993850223466\), with price about
\(-3.48823786534\). Nonnegative degree-six radial Bernstein columns also
have negative prices. Full numerical records are in
`copositive/tent_optimized_distances.json`.

These prices justify adding the columns to a subsequent optimization.
They do not quantify the eventual density improvement and are not global
feasibility certificates.
