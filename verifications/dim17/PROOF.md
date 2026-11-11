# Solid-harmonic cyclic certificate in dimension 17

The primary exact/interval verifier proves an unrestricted packing-density bound

Delta_17 < 0.018335552.

The density enclosure is centered at
0.01833555157029049139332536094236693008327772532524438023943130167025927979684,
with Arb radius below 3.20e-78 in the original replay. An independently implemented audit also passed; see ../dim17_cyclic_independent_20261111/REVIEW.md. This is not external peer review.
This improves the pinned published upper 0.0184640903350649 by 0.696153088% using
the outward decimal. No worldwide-priority or optimality claim is made. The
comparison source is Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini,
https://arxiv.org/abs/2006.02560v2, as recorded in the prior repository experiments.
The independent review checked recent indexed sources; its freshness caveat remains.

## Mathematical input and normalization

This certificate uses a cyclic generalization, not the individual-kernel-sign hypothesis of the printed Theorem 1.1. Use the cyclic three-point summation mechanism of Cohn, de Laat and Salmon,
*Three-point bounds for sphere packing*, Theorem 2.4 and equation (2.4),
https://arxiv.org/abs/2206.15373. The argument below also states the required
summation directly. This is a new computational witness in an existing framework,
not a claim of a new density theorem.

Set n=17, s=r^2=17/8 and T=R^2=51/20. Let P and H be the exact rational
polynomials in certificate.json and N=H(0)>0. Define

f(x)=P(pi |x|^2) exp(-pi |x|^2)/N.

The exact Gaussian-Laguerre Fourier identity gives
fhat(xi)=H(pi |xi|^2) exp(-pi |xi|^2)/N. The verifier checks this identity
coefficient by coefficient, checks H>0 on the nonnegative ray, and P<0 for
its argument at least pi*T using exact rational Sturm sequences and a rigorously
lower rational enclosure for pi. Thus fhat>=0, fhat(0)=1, and f<=0 beyond R.

## PSD kernel with solid spherical harmonics

Let w(u) be the five piecewise-linear cardinal hats at u_i=s+(T-s)i/4 on
[s,T], zero outside that interval. At zero define phi_0(0)=(1,0,...,0),
while for u=|x|^2 in [s,T] define phi_0(x)=(0,w(u)); otherwise use zero.
For nonzero annular x,y with u=|x|^2, v=|y|^2, t=xhat dot yhat, set

K(x,y) = [w(u)^T A_0,ann w(v)
 + sum_(ell=1)^16 w(u)^T A_ell w(v) (sqrt(u*v)/s)^ell P_ell^(17)(t)]/N.

The degree-zero term on all arguments is phi_0(x)^T A_0 phi_0(y)/N.
Positive-degree terms vanish if either argument is not in the annulus.
All 17 rational matrices are positive definite by exact Schur elimination.
The normalized Gegenbauer addition theorem and the radial factors
(u/s)^(ell/2) prove K is globally positive semidefinite. Neither the isolated
origin nor the support cutoff must be continuous for the finite kernel sums.

The solid radial factors are the representation change from the prior generic
annular engine. They turn the triangle constraints into rational polynomials in
squared side lengths without any square-root denominators.

## Polynomial triangle identities

For squared side lengths u,v,w define
Q_ell(u,v,w)=(sqrt(u*v)/s)^ell P_ell^(17)((u+v-w)/(2sqrt(u*v))).
These are exact polynomials, using

Q_0=1, Q_1=(u+v-w)/(2s),
Q_(ell+1)=[(2ell+n-2)(u+v-w)Q_ell/(2s)
            - ell*u*v*Q_(ell-1)/s^2]/(ell+n-2).

On each radial-hat cell the unnormalized kernel is therefore a rational
polynomial k(u,v,w) of total degree at most 18. For every compact triple of
squared sides u,v,w in [s,T], the verifier proves

k(u,v,w)+k(u,w,v)+k(v,w,u)<0.

There are 20 unordered triples of radial cells. Every initial tensor-Bernstein
coefficient in each cell is strictly negative, so no compact-cell subdivision
was necessary. Since R<2r, every such triple of side lengths is a triangle;
verifying the full cubes would be safe even without this fact.

If exactly one side lies outside the support, only the kernel at the opposite
vertex remains. The verifier proves k(u,v,w)<0 for u,v in [s,T], w>=T,
up to the maximum feasible w <= 4T. There are 10 unordered pairs of radial
cells. Exact rational Bernstein initialization and outward-rounded dyadic
de Casteljau subdivision certify all of them in 96 boxes. This witness needs
no triangle-domain exclusion: the entire larger rectangular region passes.
If two or three sides are outside support, every vertex term vanishes.
Compact and far certificates both include support boundary points, so the
closed-support convention leaves no uncovered case.

## Coupled radial inequality

On each radial cell, q(u)=N[2K(x,0)+K(x,x)] is a rational polynomial:
the degree-ell diagonal multiplier is (u/s)^ell. The verifier proves

P(pi*u) exp(-pi*u)+q(u)<0 for every u in [s,T].

A 256-bit Arb Taylor estimate bounds each interval by
F(m)+|F'(m)|h+sup|F''|h^2/2. Polynomial second derivatives are enclosed using
centered Taylor coefficients, not sampling. All 24 boxes pass, maximum depth 5.
Beyond R both kernel terms vanish and the already verified tail sign applies.

## Density conclusion

For a periodic packing with M centers per period cell of volume V, sum the PSD
kernel over the finite displacement set consisting of zero and all supported
neighbors around each center. Off-diagonal terms group into twice the cyclic
sums over unordered triangles and are nonpositive. Consequently the sum over
ordered distinct centers of 2K(x-y,0)+K(x-y,x-y) is at least -M K(0,0).
The pair inequality then bounds the ordered off-diagonal sum of f by M K(0,0).
Poisson summation and fhat>=0, fhat(0)=1 imply

M/V <= f(0)+K(0,0) = [P(0)+A_0[0,0]]/N.

Multiplying by pi^(17/2)/Gamma(19/2)*(s/4)^(17/2) yields the stated density.
Periodic approximation extends this to arbitrary packings. No lattice-only
assumption or external neighbor-count estimate is used.

## Reproduce and limits

Install requirements.txt. Run without Python -O:

python verify.py

Verification requires only NumPy and python-flint. The search additionally uses
SciPy, CVXPY and Clarabel. Exact witnesses, not optimizer statuses, define the
functions. The verifier trusts python-flint/Arb, exact rational arithmetic and
correctly rounded IEEE-754 NumPy operations for outward interval subdivision.
This is not a formal proof-assistant certificate or external peer review.

The initial candidate was rejected by a genuine far-domain positive value at
squared sides (2.125,2.125,8.16796875). The rejected witness and logs are retained.
Checking all real angular stationary points repaired discovery before the final
complete exact/interval verification. Numerical pairwise controls and 16D/22D
experiments are not certified results and must not be registered as bounds.
