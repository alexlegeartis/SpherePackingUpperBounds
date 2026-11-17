# Independent audit of the 43D compact harmonic kernel

**Status: COMPLETE PASS. Main geometry, exact algebra, caps, fresh scalar closure, tail and normalization have passed the described independent audit.**

Only frozen input SHA256 `0619be767fbb2d7b15b157af8bd282126425442d1339ce24ec33f250aed61615` is covered.

## Kernel and plane decomposition
In dimension43, use normalized Gegenbauer polynomials of parameter41/2, even degrees0 through18. Let R=19/13, h(t)=max((19-13t)/6,0), H=(h^6,h^9,h^12), c0=1, ck>=0, q=sum ck. Each PSD matrix Bk gives the positive-definite harmonic feature kernel K(x,y)=q^-2 sum ck G_(2k)(cos(x,y)) H(|x|²)^T Bk H(|y|²). The origin extension K(0,x)=z^T H/q and K(0,0)=s is PSD because the constant harmonic block [[s,z^T],[z,B0]] is PSD.

A Gauss-reduced saturated lattice-plane basis has squared norms1<=a<=b and inner product0<=d<=a/2. Apart from the three root pairs ±v, ±w, ±(w-v), every lattice vector has squared norm>=2a>=2>R. For three supported lines, c=a+b-2d satisfies b<=c<=R. The normalized cube substitution a=1+(6/13)x, b=a+(6/13)(1-x)y, c=b+(6/13)(1-x)(1-y)z gives feature arguments (1-x)^3, (1-x)^3(1-y)^3, and (1-x)^3(1-y)^3(1-z)^3. The other pair weights are (1-z)^6 and ((1-y)(1-z))^6. These are the formulas checked by the interval gate. At profile zeros the original pair terms vanish; the normalized formula extends continuously with positive a,b,c denominators.

Two supported lines have 0<=cos<=1/2. Independently constructed quadratic Bernstein controls (1,a,a²), (1,(a+b)/2,ab), (1,b,b²) on eight equal subintervals produce17 distinct vectors. All153 unordered control-pair angular polynomials were independently checked by direct cosine Sturm root counts on[0,1/2], with strict margin1e-7. This differs from the producer's cosine-squared representation.

At most one opposite pair on each primitive line is supported because a second multiple has squared norm>=4>R. Its four ordered collinear terms plus origin terms give G(t)=2/q (H^T W H+z^T H), W=sum ckBk/q. With nonpositive independent-plane sums, positivity of the full kernel implies sum_(v!=0) G(|v|²)>=-s. The complete main trace has now independently established the three-line signs:212423 nodes,106212 accepted leaves and an empty final stack at384bits. All incremental prefixes were SHA-guarded against the final trace, whose SHA256 is `44ee0c9f5c8bbc31782bd2e1ba12a1557e7aa707dc690d495a979754ff0ad9d4`.

## Exact input and repairs
The numerical discovery used quadratic Bernstein radial coordinates. The independently checked rational transformation is T=[[1,-2,1],[0,2,-2],[0,0,1]], B_monomial=T^T B_Bernstein T, z_monomial=T^T z_Bernstein. Chosen rational repair subtracts1/500 only from B0[0,0] and1/1000 only from z0. All final matrices pass exact Sylvester tests. Exact origin slack is1e-12. Discovery Gaussian/profile pieces are absent from the final proof; no fact about their numerical amplitudes is assumed.

## Scalar closure and normalization
Contact coordinate is x=17|v|², so h=(19/6)-(13/102)x and compact cutoff is323/13. Degree81 P and Q are related by the exact43D Fourier involution on polynomial times Gaussian: Laguerre parameter41/2, L1=43/2-x, recurrence coefficient39/2. Independent checker verifies both transform directions, Q(0)=1, and exact factorization into positive negative-root linear factors, nonnegative squares and strictly positive conjugate-pair quadratics.

Fresh isolated replay checked653 contiguous annuli from17 to180, including323/13 as a boundary, using761 interval nodes and1024-bit Arb. The replayed budget is exactly0. Therefore P(x)e^-x+G(x/17)<=0 on the entire finite interval, not merely after a weighted counting estimate. All81 physical roots were freshly enclosed; every possibly real root is below180, P(180)<0 and the leading coefficient is negative. The infinite tail is negative. Consequently no mod-2 counting theorem is actually needed for this zero-budget result; nevertheless all stored n43 cap counts and the2^44 aggregation factor were independently recomputed.

With the complete main geometry now verified, Poisson summation yields density <=(17/4)^(43/2)/Gamma(45/2) (P(0)+s). The independent odd-dimensional comparison squares this expression, uses Gamma(45/2)=44! sqrt(pi)/(4^22 22!), and proves the strict bound2.5363822051e-6 against a rational Machin lower bound for pi. No floating density approximation establishes the final rounding.

## Independence and completed coverage
Exact algebra and cap checks are separately implemented. Fresh scalar replay uses an isolated snapshot but shares the producer's Arb/Taylor code; main replay uses an independent preorder coverage parser at384bits and shares the reviewed Jet/Arb interval backend. No independent arithmetic implementation is claimed. Full main replay and empty-stack coverage have passed. The stored maximum-upper field is a conservative ceiling after prefix resume, not a measured optimized slack.
