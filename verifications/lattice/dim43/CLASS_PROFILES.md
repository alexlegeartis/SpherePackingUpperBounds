# Standalone class, weighted-rank and AA foundations for the 43D certificate

All squared norms below use minimum squared lattice norm1. A nonzero class of
Lambda/2Lambda contains opposite pairs, since v and -v have the same class.
Representatives of distinct opposite pairs v,w satisfy
 |<v,w>| <= (|v|^2+|w|^2-4)/2,                         (1)
because (v+w)/2 and (v-w)/2 are nonzero lattice vectors. In particular there
is at most one opposite pair of squared norm<2 in a class. The zero class
has no nonzero vector of squared norm<4.

The certificate uses disjoint half-open squared-norm bins [E_j,E_{j+1}),
with E_0=1, a boundary at2, and last endpoint9/4. Exactly2 is a long vector;
exactly9/4 belongs to the nonpositive tail. Source prefixes mean [1,E_p),
not a second assignment of boundary vectors. Closed intervals appearing in
polynomial or angular checks are conservative containing intervals only.

## 1. Cumulative class capacities B and A

For N unit vectors representing distinct lines in R^d with all absolute
pairwise inner products at most c and dc^2<1, their Gram matrix satisfies
 N^2 <= d*tr(G^2) <= d[N+N(N-1)c^2].
Therefore N <= d(1-c^2)/(1-dc^2). Taking the integer floor is valid.

For opposite pairs in one class with squared norms<t<=9/4, (1) gives
normalized absolute inner products at most c=1-2/t. The function
(a+b-4)/(2sqrt(ab)) increases in either a,b in the relevant interval[1,9/4]
whenever a+b>=4; if a+b<4, (1) makes the pair impossible. Thus the cumulative
number of pairs is at most
 B_n(t)=floor[n(1-c^2)/(1-nc^2)],                       (2)
whenever its denominator is positive. The verifier checks this denominator.

Suppose the class also has its unique short pair with squared norm<a<2,
where a is a bin upper endpoint. If a+t<=4, there are no other pairs below t:
the actual short and long squared norms have sum strictly less than4. This
explains the non-strict endpoint test in the code; it would not be valid for
closed-bin counts at actual equality4.

Otherwise let beta=(a+t-4)^2/(4at). The normalized correlation of a long
representative with the normalized short vector is at most sqrt(beta), by
(1) and monotonicity. Project the normalized long representatives onto the
orthogonal complement of the short direction, then normalize. Their absolute
correlations are at most
 c'=(c+beta)/(1-beta).
Indeed the projected numerator is bounded by c+beta and the product of the
projected norms is at least1-beta. If (n-1)c'^2<1 the preceding relative
bound applies in dimension n-1. In all cases the short pair already occupies
one of the B_n(t) pairs. Thus the number of long pairs below t is bounded by
 A_n(a,t)=min(B_n(t)-1,
             floor[(n-1)(1-c'^2)/(1-(n-1)c'^2)])        (3)
when the second denominator is positive, and by B_n(t)-1 otherwise. Use zero
when a+t<=4. All denominators and case conditions are checked exactly.
Taking a backward minimum of cumulative capacities at later cutoffs remains
valid, because cumulative counts are nondecreasing.

## 2. Weighted-rank class constraints

These are finite-grid theorems for the exact cutoffs used in the certificate;
the verifier regenerates their exact proofs, not merely their numerical values.
For a cutoff2+d, 0<d<=1/4, set
 alpha=(n-1)d/[4(1+d)].
For all-long pairs with squared norms2+delta_i, delta_i in[0,d], the constraint
 sum_i(1-alpha*delta_i)<=n                            (4)
is verified as follows. Put N=#pairs, S=sum delta_i, Q=sum delta_i^2. Equation(1)
and the trace/Frobenius rank inequality imply
 (2N+S)^2 <= n[4N+4S+(N/2)Q+S^2/2].                 (5)
For clarity, the off-diagonal square sum is bounded by
 [(N-2)Q+S^2]/2; adding the diagonal sum4N+4S+Q gives (5).
Here Q<=dS. If N<=n, (4) is immediate. Otherwise normalized correlations are
at most d/(2+d), so (2) gives the finite cardinality upper bound
 N<=4n(1+d)/[4+4d-(n-1)d^2].                          (6)
Let S*=(N-n)/alpha. Substitution Q=dS into the right minus left of(5) gives
 F(S)=n[4N+4S+(N/2)dS+S^2/2]-(2N+S)^2.
Its quadratic coefficient is n/2-1>0. For each possible integer N, the exact
verifier checks F(0)<0 and F(S*)<=0. Convexity then gives F(S)<0 for every
0<=S<S*, contradicting(5). This proves(4) for each used cutoff.

If a short pair of squared norm2-z is included, and the other N-1 pairs are
long and at most2+d, equation(1) gives every long delta>=z. Consequently
 S>=(N-2)z, S<=(N-1)d-z,
 Q<=(d+z)S+2z^2-(N-2)dz.                             (7)
The last inequality follows from delta^2<=(d+z)delta-dz for delta in[z,d],
then adding the short contribution z^2. The same rank bound(5) applies.

For the narrow conditioned rows, restrict0<=z<=1/96 and use the same alpha.
If N>=2 and N<=n, S>=0 proves the total weighted bound. The singleton case
is also trivial:1+alpha*z<=n. For N>n, (6) still holds: short-long normalized
correlations do not exceed d/(2+d), and long-long correlations have that bound.
To exclude S<S*, substitute(7) into(5) and consider the polygon
 0<=z<=min(d,1/96), (N-2)z<=S<=(N-1)d-z, S<=S*.
The resulting upper polynomial is convex in S and, as checked exactly, along
each of the five possible polygon edge directions. Move any point horizontally
to a convex combination of boundary points, then along each edge to its vertices.
The verifier checks that every vertex has value<=0, and every vertex with
S<S* has value<0. Since such convex combinations preserve S, every point with
S<S* has strictly negative value. This contradicts(5).
Thus the total weighted bound includes the short contribution1+alpha*z,
and the long-only right-hand side is n-1-alpha*z.

For the all-short rows, replace alpha by alpha'=(21/20)alpha and allow any
short squared norm in[1,2), or0<z<=1. If z>d no long pair is possible, and the
claimed right-hand side n-1-alpha'*z is nonnegative. Otherwise use the preceding
polygon argument on0<=z<=d, with one additional necessary condition from(3).
Writing N-1 for the number of projected long lines, admissibility requires
 c'(z)^2 >= (N-n)/[(n-1)(N-2)] for N>n.                (8)
The projected correlation c'(z) decreases with z. Cases failing(8) already at
z=0 are impossible. Otherwise exact rational bisection finds an excluded upper
endpoint, used as a conservative additional polygon boundary. The verifier
again checks every vertex and every edge curvature. It does this for every
used cutoff and cardinality. The separate all-short proof document explains
the same finite argument in detail. This is not an assumed all-real-cutoff
claim. For N<=n the previous elementary reasoning applies unchanged.

Finally, in a radial bin replace a long delta by its bin upper endpoint. Since
alpha is nonnegative this lowers its weight, preserving the upper inequality.
Replace the short delta=-z by its upper endpoint minus2, which increases the
right-hand side and is likewise conservative. Thus the exact profile rows are
precisely those used in the LP dual replay.

## 3. Projected decomposition caps for AA cuts

Let A be the classes whose unique short representative has squared norm in
[1,U), with1<U<5/4, and N=|A|. For a fixed short target orientation w with
squared norm t in[a,b] subset[1,2], consider all decompositions w=u+v with
both source norms in[1,U]. Their first endpoints form a set invariant under
u -> w-u. Distinct points are lattice vectors at mutual distance at least1.
Centered at w/2, this is an antipodal set +/-y with
 1-t/4 <= |y|^2 <= U-t/4.
For two distinct centered lines, both distances between representatives and
their negatives are at least1. Normalizing gives absolute correlations at most
 s=1-1/[2(U-a/4)].
This follows from (r1+r2-1)/(2sqrt(r1*r2)), increasing in either squared radius
when both are below1. Here U-a/4<1, as required.

Also |y dot w|<= (U-1)/2. Therefore the squared correlation of y with w is at
most epsilon=(U-1)^2/[a(4-a)], using t<=2 and the increasing function t(4-t).
Projecting onto w-perpendicular and normalizing yields a projective code in
dimension d=n-1 with absolute correlations at most
 c=(s+epsilon)/(1-epsilon).                           (9)
The checks ensure epsilon<1 and a usable cap below. The projection is injective
on distinct lines because its resulting absolute correlation bound is<1.

Here is the exact projective cap family used, not an unstated numerical bound.
For an integer r>=0 define
 mu_r=(2r-1)!!/[d(d+2)...(d+2r-2)], with mu_0=1,
 theta_r=(2r+1)/(d+2r).
If c^2<theta_r, the even polynomial H(x)=x^(2r)(x^2-c^2) is nonpositive on
[-c,c], and all its normalized dimension-d Gegenbauer coefficients are
nonnegative, with constant coefficient mu_r(theta_r-c^2)>0. One way to verify
the coefficient assertion: the expansion of x^(2r) has nonnegative coefficients;
for a common degree2j the coefficient ratio from x^(2r) to x^(2r+2) is
 (r+1)(2r+1)/[(r+1-j)(d+2r+2j)].
Its minimum is at j=0, because the denominator's difference from j=0 is
j(2-d)-2j^2<=0. The new highest coefficient is positive. The constant coefficient
is the spherical moment displayed above.

The spherical addition formula then gives, for L distinct lines,
 L^2 H_0 <= sum_{i,j}H(<x_i,x_j>) <= L H(1).
Hence L <= (1-c^2)/[mu_r(theta_r-c^2)]. The ordered decomposition-point count is
at most twice the integer floor of this bound. Taking the smallest valid bound
for r=0,...,8 is allowed; all thresholds and moments are checked rationally.
This explains the stored C_j, including its factor2.

## 4. From a scalar angular polynomial to a global AA population cut

Take F(x)=sum_k b_k G_{2k}^{(n)}(x), b_0=1,b_k>=0, and H=F(1). Require F<=0
on[0,U-1]. If two distinct source representatives have both squared sum and
squared difference at least2, their normalized absolute correlation is at most
(U-1)/U<=U-1. Thus a positive F term forces the shorter sum/difference to be
a nonzero vector of squared norm<2, the unique short target pair in its class.
Fixing a target orientation gives unique source orientations with u+v=w.
There is no assertion that arbitrary additive class triples have this property;
only positive contributions are charged.

For a target bin[a,b], the absolute source correlation lies in the safe interval
 [max(0,(2-b)/(2U)), min(1-1/(2U),(2U-a)/2)].
The stored M_j>=0 majorizes F throughout this closed interval, checked by exact
polynomial sign tests. The preceding ordered decomposition cap C_j therefore
bounds the positive off-diagonal contribution for each target class by C_j M_j.
Summing the addition-formula inequality gives
 N^2 <= H N + sum_j C_j M_j m_j,
where m_j counts ALL short target classes in binj, including bins outside the
source prefix. Classes without a short representative incur no positive charge.
For any fixed N0>0, N^2>=2N0*N-N0^2 produces the linear inequality
 sum_j [ (2-H/N0)*1_{j in source prefix} - C_j M_j/N0 ] m_j <= N0.  (10)
This is the precise stored weight vector and constant. Different source-prefix
cuts may overlap: each inequality is independently valid and their combination
uses nonnegative multipliers. No disjointness between separate cuts is required.

## 5. Profile duals and final aggregation

For each short-bin branch (and the branch with no short pair), let x_b>=0 be
long opposite-pair counts. Rows(2),(3),(4) and the applicable conditioned rows
form a matrix inequality A*x<=r. A stored nonnegative rational dual vector y
with A^T*y>=g_long proves g_long dot x<=y dot r, without an integrality or
optimization assumption. Add the short height g_j if present and subtract the
nonnegative combination of AA weights from(10). The resulting branch upper
bound is checked exactly; kappa is their maximum, covering every class.

There are m=2^n-1 nonzero classes. Their opposite pairs contribute twice the
per-class heights, so summing(10) and the class budgets gives
 sum_{v!=0}G(|v|^2) <= 2m*kappa + 2 sum_l nu_l*N0_l.
The cyclic-collinear bridge consequently yields
 rho*Fourier(f)(0) <= f(0)+C+2m*kappa+2 sum_l nu_l*N0_l.
This is the exact objective used. It requires no prior density upper bound,
no quotient classification, and no low-density alternative branch. Every
continuous sign and mixed-exponential envelope is a separate verified gate.

For this 43D package, replace the final objective origin term C by r²C with r=7/10, and use the modified cyclic correction from ORIGIN_RATIO_PROOF.md. The class capacities and AA inequalities themselves are unchanged.
