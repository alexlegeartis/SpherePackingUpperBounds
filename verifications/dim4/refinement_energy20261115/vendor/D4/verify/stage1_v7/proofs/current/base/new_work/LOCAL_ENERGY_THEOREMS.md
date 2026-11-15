# Local energy inequalities for four-dimensional sphere packing

These are analytic derivations obtained during this investigation. They do not
by themselves prove an upper bound below 0.63. The accompanying exact data and
verifiers distinguish accepted numerical certificates from search proposals.

## 1. The global theorem needs positivity on actual neighborhoods

Normalize the minimum distance between packing centers to 1. Let f be a real,
even Schwartz function on R^4, with nonnegative Fourier transform and
I=fhat(0)>0. Use Fourier exponent -2 pi i x.xi. Let K be a bounded, symmetric
kernel, supported when each nonzero argument has radius between 1 and a fixed
finite R. Write K(0,0)=A.

For every finite set C of separation-one points in this annulus, assume

    E(C)=A+2 sum_(x in C) K(x,0)+sum_(x,y in C) K(x,y) >= 0.       (1)

Positive semidefiniteness of the full kernel implies (1), but the theorem only
uses (1). In addition require

    f(x)+2K(x,0)+K(x,x) <= 0                     (|x|>=1),       (2)
    K(x,y)+K(-x,y-x)+K(-y,x-y) <= 0              (all sides>=1). (3)

For a finite packing, sum (1) over its centers, using all neighbors in the
kernel support. A pair contributes 2K(x,0)+K(x,x), and the terms indexed by
three distinct centers group into cyclic sums (3). Therefore

    sum_(i,j) f(x_i-x_j) <= N[f(0)+A].                         (4)

For periodic packings, the same argument is applied per period and Poisson
summation gives N^2 I / V as a lower bound for the scalar sum per cell. Hence
their center density is at most [f(0)+A]/I. Arbitrary packings have the same
upper bound: repeat a large finite patch with a separating corridor, and let
the boundary-to-volume ratio tend to zero. Bounded finite support of K and
Schwartz decay of f make the boundary errors vanish. The ordinary packing
density is consequently at most

    (pi^2/32) [f(0)+A]/I.                                    (5)

The finite grouping and periodic-to-arbitrary-packing argument are also
written out in the retained d4_research/UPPER_BOUND.md. No lattice assumption
is used for the upper bound.

## 2. Integer neighbor counts provide additional local energies

Let N_j count neighbors in a specified radial interval, and let
q(N_1,...,N_k) be a quadratic polynomial nonnegative on every possible integer
count vector. Its constant, linear, and quadratic terms can be represented by
an additional kernel. Its local energy is exactly q(N), so it may be added to
(1), whether or not that additional kernel is positive semidefinite on all
finite sets of vectors.

For a single closed ball cutoff R0, write chi(r)=1_[1,R0](r).
The polynomial q(N)=a+bN+cN^2 has the representation

    K_q(0,0)=a,
    K_q(x,0)=b chi(|x|)/2,
    K_q(x,y)=c chi(|x|)chi(|y|).

Its pair contribution is (b+c)chi, its cyclic contribution is the sum of the
three corresponding products multiplied by c, and its objective contribution
is a. When N<=M is known, exact conditions q(n)>=0 for n=0,...,M suffice.

The accepted global count certificate uses the simpler energies

    lambda_j(M_j-N_j) >= 0, lambda_j>=0.

They subtract lambda_j from the pair expression inside the corresponding
radius and add M_j lambda_j to the objective numerator. They do not change
the cyclic expression. In the normalization f=e^(-2r^2)F(2r^2), G(0)=1,
I=pi^2/4, the resulting density bound is the exact rational number

    [F(0)+kappa+sum_j M_j lambda_j]/8.

In particular, exact verification of the supplied scalar, a positively
scaled supplied nonorigin kernel, a new rational origin row, and the two
neighbor bounds 24 at 201/200 and 25 at 100557/100000 gives

    Delta_4 <= 0.6345865969160774... < 0.634587.

The displayed long decimal is an approximation to the exact rational value
stored in counts/exact_fixed_scalar_count/data.json. The strict rounded
upper bound is verified by rational comparison.

### Count-times-moment inequalities

If N<=M and w is a nonnegative radial function, then

    (M-N) sum_(x in C) w(|x|) >= 0.

This is represented by

    K_w(0,0)=0,
    K_w(x,0)=M w(|x|)/2,
    K_w(x,y)=-(chi(|x|)w(|y|)+chi(|y|)w(|x|))/2.

Its pair contribution is (M-chi)w. It supplies negative cyclic terms even
when only one of the triangle's edges lies in the counted interval. Each
nonnegative Bernstein basis function is an exact polynomial choice of w.

## 3. Transfer an entire spherical certificate into a local energy

Let C be a code of N unit vectors on S^3 with distinct inner products at most
s<1. Suppose a spherical scalar g and an anchored kernel J satisfy the
spherical pair and cyclic conditions stated in spherical_sdp/THEOREM.md.
Let kappa_s be the origin entry of J, and put B=g(1)+kappa_s. Summing the
anchored nonnegative energies gives

    sum_(u,v in C) g(u.v) <= B N.                             (6)

This implication uses the positivity of the anchored local energies and
their pair/cyclic signs. It does not use nonnegativity of the Gegenbauer
coefficients of g. That latter condition is needed only if (6) is subsequently
used to infer the usual numerical code-size upper bound.

Consequently

    R_g(C)=B N-sum_(u,v in C)g(u.v) >= 0                      (7)

is a valid local quadratic energy retaining angular information.

Now take a separation-one packing neighborhood whose radii are in [1,R0].
For 1<=R0<2, its projected directions have inner products at most

    s(R0)=max(R0/2, 1-1/(2R0^2)).

Indeed the maximum of (a^2+b^2-1)/(2ab) over [1,R0]^2 occurs at a corner:
holding either radius fixed gives only an interior minimum. For R0 below
the golden ratio, the second term in the displayed maximum is the larger.

Whenever s(R0)<=s, use (7) with these projected neighbors. With chi as above,
the resulting four-dimensional local kernel is

    K_g(0,0)=0,
    K_g(x,0)=B chi(|x|)/2,
    K_g(x,y)=-chi(|x|)chi(|y|)g(x.y/(|x||y|)).

Its pair contribution is exactly kappa_s chi, its objective contribution is
zero, and its cyclic terms are the corresponding negative angular values of
g. Thus it may be combined with a global packing certificate via (1)-(5).

For a repaired spherical certificate, g here is the repaired scalar: its
constant coefficient is decreased by the exact repair amount. Kappa_s is
unchanged. Using the unrepaired scalar would not be justified.

The anchored proof of (6) involves three projected neighbors in addition to
the original packing center. This is the source of the information about four
centers. Co-optimizing g and J together with the global certificate is a
possible larger search space. No claim that this co-optimization reaches
0.63 is made. The accepted spherical certificate at s=.5075 happens to have
very small nonconstant scalar coefficients, so its particular transferred
ray largely reduces to a count inequality; that diagnostic is recorded to
avoid overstating its computational contribution.

## 4. Ball counts and cross-cap energies

If n separation-one points lie in a Euclidean ball of radius a, then

    n(n-1)/2 <= sum_(i<j)|x_i-x_j|^2
              =n sum_i |x_i-xbar|^2 <= n^2 a^2.

Thus a^2<M/[2(M+1)] implies n<=M. In R^4 one also has n<=5 whenever
a^2<1/2: the vectors relative to the ball center have strictly negative
mutual inner products. Six such vectors would have an affine dependence;
separating its positive and negative coefficients makes the squared norm of
their common weighted sum a strictly negative cross inner-product sum.

Choose two probe balls B(dv,a) and B(dw,a), with v.w fixed, and average over
rotations. If each first ball contains at most M eligible neighbors, then

    (M-N_v)N_w >= 0.

Here eligible neighbors are restricted to the global kernel's active
annulus. Write chi(x) for its indicator. Let alpha(r) be the single-probe
membership probability multiplied by chi(r), and let J(x,y) be the averaged
joint membership in the two different balls, symmetrized in x,y and
multiplied by chi(x)chi(y). Thus the nonorigin kernel is zero if either
argument is outside the active annulus. The average energy has kernel

    K(0,0)=0, K(x,0)=M alpha(|x|)/2, K(x,y)=-J(x,y).

The pair cost is M alpha(r)-J(x,x); its cyclic correction is nonpositive.
Opposite and orthogonal probe directions can supply negative terms at
exterior triangles. These inequalities are geometrically exact; any
numerically optimized combination must still have its resulting global pair
and cyclic inequalities independently verified.

For opposite probes, the additional nonnegative square (N_v-N_(-v))^2
gives a positive semidefinite kernel consisting entirely of odd spherical
harmonics. It has zero energy on a centrally symmetric neighborhood such as
every centered finite radial truncation of D4. This supplies a useful
comparison direction beyond a fixed finite angular truncation.
