# A stronger conditional count by affine latitude compression

For a packing whose center separation is at least one, write \(N_r\) for
the number of other centers at distance at most \(r\) from a selected center.
The following two local implications are verified by exact rational
arithmetic:

\[
\begin{aligned}
N_{1.005643}\ge24&\quad\Longrightarrow\quad N_{1.01534}=24,\\
N_{1.002946}\ge24&\quad\Longrightarrow\quad N_{1.02597}=24.
\end{aligned}
\]

The antecedent also implies that its count is exactly 24. More directly,
24 selected inner neighbors cannot coexist with any additional neighbor
inside the stated outer radius. The statements include all radial endpoints.
They depend on the retained exact spherical-code certificate
\(A(4,203/400)\le24\), not on uniqueness or stability of the 24-cell.
They are local results and do not themselves establish a global density bound.

## 1. Radial geometry and the missing pole

Suppose there are 24 neighbors of radii in \([1,r]\) and an additional
neighbor of radius at most \(R\), with \(1\le r<R<\sqrt2\).
Their directions are distinct. Let \(p\) be the direction of the additional
neighbor. The inner directions \(u_i\) satisfy

\[
u_i\cdot u_j\le s=1-\frac{1}{2r^2},\qquad
u_i\cdot p\le c=\frac{r^2+R^2-1}{2rR}.
\]

These inequalities follow from separation one and monotonicity of
\((a^2+b^2-1)/(2ab)\) in both radii on \([1,R]^2\).

Put \(S=203/400\). The task is to move the 24 inner directions into the
cap complement with pole coordinates at most \(S\), while keeping every
inner pair's dot product at most \(S\). The unchanged pole \(p\) would
then complete a forbidden 25-point spherical code.

## 2. An affine latitude map

For \(a=u\cdot p\le c<1\), set

\[
\beta=\max\left(0,\frac{c-S}{1+c}\right),\qquad
\alpha=1-\beta,\qquad
\phi(a)=\alpha a-\beta.
\]

Here \(0<\alpha\le1\), \(\phi(-1)=-1\), and \(\phi(c)\le S\).
Writing \(u=ap+u_\perp\), define

\[
T(u)=\phi(a)p+k(a)u_\perp,
\qquad
k(a)^2=
\alpha\frac{(1-a)+\beta(1+a)}{1-a}.
\tag{1}
\]

The denominator is positive on the entire cap complement. The identity
\(k(a)^2(1-a^2)=1-\phi(a)^2\) shows that \(T(u)\) has norm one.
Formula (1) is well-defined at \(a=-1\), where it returns \(-p\).
The latitude map is strictly increasing, and the perpendicular multiplier
is positive, so \(T\) is injective on the cap complement. Its image has
latitude at most \(S<1\), and therefore contains no point equal to \(p\).
The extra pole is kept fixed; (1) need not be defined at the northern pole.

For two inner directions with latitudes \(a,b\) and original dot product
\(t\), the transformed dot product is

\[
T(u)\cdot T(v)=\phi(a)\phi(b)+k(a)k(b)(t-ab).
\tag{2}
\]

The coefficient of \(t\) is positive.

## 3. Why it suffices to check original dot product exactly s

For fixed latitudes \(a,b\), the allowed dot products form the interval

\[
ab-\sqrt{(1-a^2)(1-b^2)}
\ \le t\le\
ab+\sqrt{(1-a^2)(1-b^2)}.
\]

If \(s\) belongs to this interval, increasing \(t\) to \(s\) increases
(2), so it suffices to check that endpoint. If the upper endpoint is
less than \(s\), first put the points on a common meridian, which
maximizes (2) for their latitudes. Write their northern polar angles as
\(\theta_a\ge\theta_b\), and let \(\theta_s=\arccos(s)\).
The upper-endpoint condition means
\(\theta_a-\theta_b>\theta_s\). Replace \(\theta_a\) by

\[
\theta_a^*=\theta_b+\theta_s.
\]

This lies between \(\theta_b\) and \(\theta_a\), and the corresponding
latitude obeys \(a\le a^*\le b\le c\). The original dot product becomes
exactly \(s\). Since \(\phi\) is increasing, the transformed meridian
angle decreases, so its dot product increases. Thus every original pair
with dot product at most \(s\) is dominated by a feasible pair with dot
product exactly \(s\). This reduction also covers the polar endpoints.

The complete remaining domain is the compact semialgebraic set

\[
-1\le a,b\le c,\qquad
D(a,b)=1+2abs-a^2-b^2-s^2\ge0.
\tag{3}
\]

## 4. Eliminating the radicals without changing the sign

Define rational polynomials

\[
P=S-\phi(a)\phi(b),\qquad V=s-ab,
\]

and

\[
\begin{aligned}
H={}&P^2(1-a)(1-b)\\
&-\alpha^2\bigl((1-a)+\beta(1+a)\bigr)
             \bigl((1-b)+\beta(1+b)\bigr)V^2.
\end{aligned}
\tag{4}
\]

The degrees in either variable are at most one for \(P,V\), two for
\(D\), and three for \(H\). Let \(L=k(a)k(b)>0\) and
\(U=(1-a)(1-b)>0\). Equation (4) says

\[
H=U(P^2-L^2V^2).
\]

The desired inequality (2) \(\le S\) is precisely \(P\ge LV\).
Each of the following sign cases is sufficient:

1. \(P\ge0\) and \(V\le0\).
2. \(P\ge0\) and \(H\ge0\).
3. \(V\le0\) and \(H\le0\).

In the second case, \(P\ge|LV|\). In the third case,
\(|P|\le|LV|=-LV\), and hence \(P\ge-|P|\ge LV\).
This treatment keeps the sign information that would be lost by blindly
squaring the inequality.

## 5. Exact verification over the continuous domain

[affine_latitude_counts.py](affine_latitude_counts.py) constructs all four
polynomials with rational coefficients and converts them to Bernstein
coefficients on \([-1,c]^2\). The convex-hull property of Bernstein
coefficients gives exact upper and lower bounds on each polynomial over
every box. The script subdivides boxes at their midpoints using integer
de Casteljau formulas. Denominators stay positive, so only integer signs
are needed during subdivision.

A box is discarded only when its upper bound for \(D\) is strictly
negative. Otherwise it is accepted only if one of the three sufficient
sign cases holds throughout the box. Every other box is subdivided.
An unresolved depth or node limit returns `UNKNOWN`, never success.
The numerical search and the conjectured location of a worst pair are
not inputs to this proof.

| Inner radius | Outer radius | Subdivision nodes | Greatest depth |
|---:|---:|---:|---:|
| 1.005643 | 1.01534 | 12,619 | 28 |
| 1.002946 | 1.02597 | 5,575 | 30 |

Every feasible box was accepted in both runs. The same command then
replays the exact count-24 spherical certificate on a byte copy. Its
result, input hash, rational parameters, and all domain-check statistics
are retained in [affine_latitude_certificates.json](affine_latitude_certificates.json).

The larger inner radii strengthen the corresponding sufficient Möbius
tests, which used 1.004727 and 1.000729 respectively. No claim is made that
the affine thresholds are the optimal geometric thresholds.

## 6. The quadratic local energy used in the global search

For either proved pair \((r,R)\), an occupied outer shell forces
\(N_r\le23\). Therefore

\[
Q=(23-N_r)(N_R-N_r)\ge0.
\]

Let \(A(x)=1_{1\le|x|\le r}\) and
\(B(x)=1_{r<|x|\le R}\). This local energy has origin constant zero,
origin coupling \(h(x)=23B(x)/2\), and symmetric nonorigin kernel

\[
K(x,y)=-\frac12\bigl(A(x)B(y)+B(x)A(y)\bigr).
\]

Its diagonal is zero, its pair contribution is \(23B(x)\), and every
nonorigin kernel term is nonpositive. These formulas make it compatible
with the retained global local-energy framework. Adding this energy
still requires a new global pair and cyclic verification; the local
lemma alone does not lower the accepted density bound.
