# A rigorous general-packing upper bound in dimension 20

**Two independently implemented verifiers pass. No worldwide-priority claim is made.** This is an unrestricted sphere-packing density bound, not a lattice-only or kissing-number bound.

The bound is

Delta_20 < 0.006979502.

The exact interval evaluation of the constructed certificate is
0.00697950190417318581318448850239932436503950767230416847431912053614315849425,
with Arb radius below 4.45e-78. The current general-packing table https://cohn.mit.edu/sphere-packing/ lists upper bound .0071270536033763 and lower density .0040241650630681… . Against that reference, this certificate improves the upper bound by about2.07%, closing about4.7553% of the current gap, **not the requested20%**.

## The theorem and normalization

Use Cohn–de Laat–Salmon, Three-point bounds for sphere packing, Theorem1.1, https://arxiv.org/pdf/2206.15373. Set n=20, r²=5/2, R²=3. The theorem requires f2 continuous/integrable, nonnegative Fourier transform with value1 at zero, f3 positive definite as a kernel, f2<=0 outside R, the coupled annular inequality, and f3(x,y)<=0 on the admissible annular pair domain. The theorem does not require f3 to be continuous.

The JSON certificate supplies exact rational polynomials P,H and17 exact rational symmetric matrices A_0,...,A_16. Let N=H(0)>0. Define
f2(x)=P(pi|x|²) exp(-pi|x|²)/N.
Its Fourier transform is H(pi|xi|²) exp(-pi|xi|²)/N. The Fourier identity is verified through the Laguerre basis L_k^9(2X), with eigenvalues(-1)^k.

## Positive-definite kernel

There are five radial nodes u_i=5/2+i/8, i=0,...,4. Let w(u) be their usual piecewise-linear nodal interpolation vector for u in[5/2,3]. Thus w(u) has at most two nonzero entries and is linear within each radial interval.

Define the degree-zero feature vector phi_0(x):
- at the origin: (1,0,0,0,0,0);
- in the annulus: (0,w(|x|²));
- elsewhere: zero.

For x,y in the annulus, the unnormalized kernel is
phi_0(x)^T A_0 phi_0(y) + sum_(l=1)^16 [w(|x|²)^T A_l w(|y|²)] P_l^(20)(xhat.yhat).
Outside the annulus the positive-degree terms are zero; the degree-zero expression remains as defined. Divide the whole kernel by N to obtain f3.

Every matrix is positive definite by exact rational LDL^T decomposition. The normalized Gegenbauer polynomials P_l^(20), with P_l(1)=1, are positive-definite spherical kernels by the addition theorem. Products and sums of positive-definite kernels are positive definite. Hence this construction gives a valid global f3, including at the isolated origin coordinate.

## Global Fourier and exterior signs

Exact Sturm sequences show H has no positive real roots and H(0)>0. For a rational number pi_lower<pi, exact Sturm sequences show -P(X+3*pi_lower) has no positive roots and has positive constant term. Thus the Fourier transform is positive everywhere and f2 is negative for |x|>=R. pi_lower is independently enclosed below pi with Arb.

## Entire off-diagonal annular domain

Subdivide every coarse radial interpolation interval into8 parts, giving33 nodes from5/2 to3 with step1/64. The kernel's angular polynomial is bilinear in u,v within every fine radial rectangle.

For admissible radii and separation>=r,
t=cos(theta)<=(u+v-r²)/(2sqrt(uv)).
The right side is increasing in both u,v on this domain, since u-v+r²>=2r²-R²=2>0. For each fine-grid node pair(i,j), choose upper radial neighbors i+1,j+1 (clipped at the final node), and a rational c_ij upper bound for the resulting maximum cosine. Its validity is proved exactly by
4*c_ij²*u_upper*v_upper >= (u_upper+v_upper-r²)²,
with all terms positive.

At each of561 unordered node pairs, form the exact rational angular polynomial Q_ij(t). Exact Sturm sequences show it has no roots on[-1,c_ij], and its endpoint values are negative. Therefore Q_ij<0 throughout that interval.

At any actual radii in a fine rectangle, the angular polynomial is a convex combination of the four corner polynomials. Their verified angular domains contain every allowable angle for that rectangle. This proves f3<=0 on the full continuum domain, not merely the sampled nodes.

## Coupled annular inequality

On each coarse radial interval, the kernel part
2 f3(x,0)+f3(x,x)
is a rational quadratic polynomial q_i(u)/N, because all P_l(1)=1. Hence the coupled inequality is equivalent to
P(pi u) exp(-pi u)+q_i(u)<=0.

A separate verifier uses256-bit Arb arithmetic and adaptive rational interval subdivision. At midpoint m and half-width h it bounds the function by
F(m)+|F'(m)|h+sup|F''|h²/2.
All quantities are interval-enclosed. The bound on the polynomial part of F'' uses its Taylor coefficients around pi*m, with absolute-value bounds on the remainder interval. The exponential is bounded above by its value at the interval's left endpoint. The verifier examines 100 boxes and certifies 52 terminal intervals; maximum subdivision depth 7. No floating sign decision is used.

## Objective

The exact unnormalized objective divided by N is
(P(0)+A_0[0,0])/N = 1134354088304391233619/38146972656261249493.

Multiply by vol(B_(r/2)^20)=pi^10/10!*(5/8)^10. Arb gives the interval quoted above and verifies it is strictly below .006979502.

## Verification package

certificate.json contains the exact rational witness. Run verify.py for the primary verifier and independent_verify.py for the independently implemented audit. Both verify every continuum hypothesis and the objective. Only Python and python-flint are required; the discovery optimizer is unnecessary.

The primary verifier uses exact rational LDL and Sturm sequences plus256-bit Arb interval arithmetic. The independent verifier uses separate matrix elimination, a direct monomial Fourier transform, explicit Gegenbauer coefficients, integer-square-root angle bounds, and384-bit Arb arithmetic. Both certify52 terminal annular intervals. Saved logs document successful reproduction.

The optimization solver output is not part of the proof. The exact witness defines the functions. The20%-gap target has not been reached, and worldwide literature priority has not been established.
