# Pointed ball-count energies and exact origin shifts

This is an application of the global local-energy theorem in
[LOCAL_ENERGY_THEOREMS.md](../LOCAL_ENERGY_THEOREMS.md), equations (1)–(5).
That theorem requires nonnegative energy on every packing neighborhood;
the additional energy below need not be a positive-definite kernel on
arbitrary point sets. The numerical optimizer discovers coefficients.
Only the exact verifier contract at the end establishes a packing bound.

## 1. Two elementary ball-count bounds

Normalize the separation of packing centers to one. If a closed ball of
squared radius \(a^2\) contains \(n\) centers, then

\[
\frac{n(n-1)}2
\le \sum_{i<j}|x_i-x_j|^2
=n\sum_i|x_i-\bar x|^2
\le n^2a^2.
\]

Consequently, for any positive integer \(T\),

\[
a^2<\frac{T}{2(T+1)}\quad\Longrightarrow\quad n\le T.
\]

There is also a dimension-specific bound in \(\mathbb R^4\): if
\(a^2<1/2\), then \(n\le5\). To prove it, subtract the ball center and
write the resulting vectors as \(z_i\). Distinct vectors satisfy

\[
z_i\cdot z_j
\le a^2-\frac12<0.
\]

Six vectors in \(\mathbb R^4\) have a nonzero affine dependence
\(\sum_i\alpha_i z_i=0\), \(\sum_i\alpha_i=0\). Separate its positive
and negative coefficients. Their two positive weighted sums agree, but
the squared norm of that common vector is a strictly negative weighted
sum of cross inner products, a contradiction.

Let \(T\) be any proved total ball-count bound. A probe ball centered at
\(d v\), where \(|v|=1\), contains the distinguished packing origin
exactly when \(d^2\le a^2\). The number of its **nonorigin** packing
centers is therefore at most

\[
M=T-\mathbf1_{d^2\le a^2}.
\]

The current accepted repair implementation uses \(d^2>a^2\) and the
first ball-count lemma, so its stored \(M\) is also its total ball-count
bound. The more general facts above do not by themselves certify an
unimplemented numerical probe family.

## 2. A nonnegative pointed energy

Fix a packing center at the origin. Let \(\chi(x)\) be the indicator of the
global kernel's active radial annulus, and let \(S\) consist of the
nonorigin packing neighbors in that annulus. Let \(w(x)\ge0\), supported
there. For each \(y\in S\), choose \(v\) uniformly among unit
vectors satisfying \(v\cdot\widehat y=\sigma\). At \(\sigma=\pm1\)
this means the corresponding point mass. Define

\[
P(x\mid y)=\Pr\{|x-dv|^2\le a^2\mid v\cdot\widehat y=\sigma\}.
\]

Every probe contains at most \(M\) selected nonorigin neighbors, hence

\[
E(S)=M\sum_{y\in S}w(y)
 -\sum_{x,y\in S}w(y)P(x\mid y)\ge0.
\]

Omitting packing points outside the compact radial support preserves the
inequality. Symmetrization gives the augmented kernel

\[
K_p(0,0)=0,\qquad K_p(0,x)=\frac M2w(x),\qquad
K_p(x,y)=-\frac{\chi(x)\chi(y)}2
 [w(y)P(x\mid y)+w(x)P(y\mid x)].
\]

In particular, the nonorigin kernel is zero if either argument lies
outside the active annulus. Support of w alone would not imply this
convention, since the other term in the symmetric formula could remain.

Its total energy is exactly \(E(S)\). Nonnegative combinations
\(\sum_j\beta_jK_{p,j}\) are therefore valid packing-local energies.
All their nonorigin values, including their diagonals, are nonpositive.

In the currently validated parameter range
\(d>0\), \(d^2>a^2\), and \(-1\le\sigma\le0\), the threshold

\[
c(r)=\frac{r^2+d^2-a^2}{2dr}
\]

is positive. Since \(v\cdot\widehat x=\sigma\) for the self term,
\(P(x\mid x)=0\). The pointed diagonal is then zero and its pair
contribution is exactly \(M w(x)\). For other valid probe families the
diagonal may be negative, which only improves the retained pair bound.

## 3. Origin shift and objective

Let \(b(x)\) be the radial feature vector, and let \(A_0\succeq0\) be
the degree-zero nonorigin Gram matrix of the base kernel. Write its
origin coupling as \(v_0^Tb(x)\), and its mass as
\(\kappa_0=v_0^Tu_0\), where \(A_0u_0=v_0\). Suppose
\(w_j(x)=w_j^Tb(x)\), and set

\[
h=\frac12\sum_j\beta_jM_jw_j,\qquad v_1=v_0-h.
\]

Require \(h\in\operatorname{range}(A_0)\); it suffices that every
\(w_j\) with \(\beta_jM_j\ne0\) lie in that range. Solve
\(A_0s=h\), and put

\[
u_1=u_0-s,\qquad
\kappa_1=v_1^Tu_1
=\kappa_0-2u_0^Th+h^Ts.
\]

The new degree-zero augmented block is PSD, since its energy is
\((u_1+\sum_x b(x))^TA_0(u_1+\sum_x b(x))\). This also covers
singular \(A_0\). The mass is independent of the chosen solution of the
linear system, and rational input admits rational solutions when the
range condition holds.

Now add the pointed energies. Their origin couplings cancel the shift
\(-h\). In the validated zero-diagonal case the original pair expression
is unchanged; in the general case it can only decrease. Every full
triangle cyclic expression can only decrease as well. Exterior triangles
may use this decrease to repair a positive base residual, but that
inequality must be established over the whole exterior domain.

The pointed energies have zero origin constant. If the other components
of the global certificate are fixed, its density objective changes by

\[
\frac{\kappa_1-\kappa_0}{8}
\]

under the normalization used by the supplied \(\mathbb R^4\) certificate.
Applying the global theorem in `LOCAL_ENERGY_THEOREMS.md` then gives the
same bound for arbitrary packings through its finite, periodic, and
boundary-limit argument.

## 4. Exact verification contract

`build_repair.py` solves the origin system in the 11 functional radial
coordinates and embeds it into the 15 piecewise coordinates. It checks
the embedding of every Gram ray. The full verifier independently checks
\(A_0u_1=v_1\), \(v_1^Tu_1=\kappa_1\ge0\), and
\(v_1+h=v_0\) in those actual expanded coordinates. Thus rank deficiency
in the expanded representation does not evade the range condition.

For \(-1<\sigma\le0\), putting \(t=\widehat x\cdot\widehat y\) gives

\[
P(x\mid y)=\operatorname{clip}_{[0,1]}
 \frac12\left(1-
 \frac{c(r)-\sigma t}{\sqrt{(1-\sigma^2)(1-t^2)}}\right),
\]

with the geometric indicator used when the denominator vanishes.
`interval_certificate.py` encloses this expression using rational bounds
and integer square-root enclosures. It combines these with outward
integer Bernstein bounds for the homogenized exterior polynomial. The
sigma=-1 point mass is handled directly. The complete derivation and
independent exact edge-case checks appear in
[POINTED_VERIFICATION.md](../audit/POINTED_VERIFICATION.md).

`verify_augmented_certificate.py` validates every probe, checks the origin
identities, reproves the scalar Fourier and tail conditions and all pair
intervals, checks every full cyclic branch conservatively, and checks
every exterior branch with the pointed correction. It recomputes the
whole objective exactly. UNKNOWN, incomplete subdivisions, failed
identities, numerical solver output, and positive witnesses are not
acceptance. Families beyond the currently validated parameters require
an extended verifier before they can establish a new bound.

## 5. Separately audited general helper

`interval_certificate_general.py` leaves the accepted v1 helper unchanged.
It allows every \(\sigma\in[-1,1]\), uses the active nonorigin count M,
and validates the total count \(T=M+\mathbf1_{d^2\le a^2}\) by either
ball-count lemma above. An optional `total_M` field must equal that exact T.
The origin-on-boundary case is included because probe balls are closed.

For \(b=d^2-a^2\le0\), the radial threshold c(r) increases, so its
endpoint values give its exact interval extrema. For b>0, its maximum is
at an endpoint and its minimum is at an endpoint or at \(\sqrt b\),
which is enclosed by the integer square-root routine. For positive sigma,
the endpoints bounding \(c(r)-\sigma t\) reverse; the helper uses the
minimum and maximum of the two products \(\sigma t\) directly. The
sigma=+1 point mass is the exact indicator of \(t\ge c(r)\).

The same sign-dependent radical bounds and exterior Bernstein comparison
then apply. The nonpositive self diagonal may be omitted conservatively
when checking the original pair expression after the origin shift.
`audit_pointed_general.py` passed 47,566 exact radical comparisons on
5,250 rational rectangles, including positive sigma, probe balls
containing the origin, origin-on-boundary cases, and deterministic closed
ball equalities. It accepted five valid count records and rejected four
incorrect records at the strict threshold or with an insufficient count.
This is an audit of the helper, not a global certificate for any future
jointly optimized candidate; all global conditions still require replay.
