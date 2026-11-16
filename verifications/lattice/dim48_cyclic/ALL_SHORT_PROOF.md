# All-short conditioned weighted class bound: exact finite audit

Scope verified by `relaxed_refined.py`: n=13,23,38,39,40,41,42,43,44 and
cutoffs d=k/384, 1<=k<=96. This is a class inequality, not a density result.
Set alpha=(21/20)(n-1)d/[4(1+d)]. If a nonzero class modulo2Lambda contains a
short opposite pair of squared norm2-z, 0<z<=1, then its long opposite pairs
with squared norms2+delta in[2,2+d] obey
 sum_long (1-alpha*delta) <= n-1-alpha*z.

Include the short pair and let N be total opposite-pair count, S=sum delta
including delta_short=-z, and Q=sum delta^2. Pairwise halfsums imply
 |<v_i,v_j>| <= (delta_i+delta_j)/2, so every long delta>=z.
If z>d there are no long pairs; the claimed RHS is nonnegative. For N>=2,
S>=(N-2)z>=0, hence N<=n is immediate. N=1 is also immediate because
alpha*z<n-1.

For N>n we have 0<=z<=d and
 (N-2)z<=S<=(N-1)d-z,
 Q<=(d+z)S+2z^2-(N-2)dz.
Rank and trace/Frobenius give the necessary inequality
 n[4N+4S+(N/2)Q+S^2/2]-(2N+S)^2 >=0.
The upper-bound polynomial after replacing Q is convex in S and along every
edge of the feasible polygon. To prove S>=(N-n)/alpha, it suffices to show
nonpositivity at all polygon vertices with S<=(N-n)/alpha, with strictness
at vertices having smaller S, and check all edge-direction curvatures.
These exact rational tests are implemented independently of any LP solver.

Additional legitimate restriction: N<=1+A(2-z,2+d), where A is the usual
conditional projected projective-code cardinality bound. Write t=2+d,
c=1-2/t, beta=(d-z)^2/[4(2-z)t], c'=(c+beta)/(1-beta). Then
 A<=min(B_n(t)-1, floor((n-1)(1-c'^2)/(1-(n-1)c'^2)))
when the second denominator is positive, otherwise use B_n(t)-1.
Here B_n(t)=floor(n(1-c^2)/(1-n*c^2)).
A is nonincreasing in z: beta decreases, as its derivative has sign
(d-z)(d+z-4)<0; the relative bound increases in c'^2. For each integer N,
exact bisection finds an infeasible upper z endpoint, which conservatively
clips the polygon. The upper endpoint remains included in the proof polygon,
so no strict-endpoint issue is hidden. Cases already infeasible at z=0 are
skipped. A separate ordinary relative bound gives
 N<=4n(1+d)/(4+4d-(n-1)d^2), making enumeration finite.

For half-open radial bins, replace each long delta by its bin upper endpoint
and short -z by its bin upper endpoint minus2; this only relaxes the inequality.
The new constraints may be added alongside the original unconditional rows
and stronger near2 conditioned rows. Do not claim a global density improvement
until exact radial certificates and all profile duals are replayed.
