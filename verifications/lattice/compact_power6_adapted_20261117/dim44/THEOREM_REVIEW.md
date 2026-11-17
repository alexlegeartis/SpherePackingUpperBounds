# Independent audit of the 44D compact harmonic kernel

**Status: COMPLETE PASS. Main geometry, exact algebra, caps, fresh scalar replay, tail and normalization have passed the described independent audit.**

Only frozen input SHA256 `4d67f38f02a67dfba251f088e471865d2919f4b95e864aa8786e293aa657ca01` is covered. This is a dimension-specific optimized kernel, not a dimension restriction of another input.

## Geometric and PSD theorem
Use normalized even Gegenbauer harmonics in dimension44 (parameter21), powers6,9,12, support R=19/13 and H=(h^6,h^9,h^12), h(t)=max((19-13t)/6,0). With c0=1, ck>=0, q=sum ck and PSD matrices Bk, the nonzero-vector kernel is K(x,y)=q^-2 sum ck G_(2k)(cos(x,y))H(|x|²)^T Bk H(|y|²). The origin extension K(0,x)=z^T H/q, K(0,0)=s is PSD by its exact constant-harmonic block.

A Gauss-reduced lattice plane has squared norms1<=a<=b and dot0<=d<=a/2. Only the three signed pairs ±v, ±w, ±(w-v) can be supported because all other vectors have squared norm>=2a>=2>R. The normalized three-line cube and outer powers6/feature stride3 are identical to the mathematical formulas reviewed for the43D package, with the44D angular polynomial substituted. Every independent pair belongs to one saturated plane. The complete trace has independently proved the normalized three-pair sum nonpositive on this cube:175733 nodes,87867 accepted leaves and an empty final stack at384bits. Every prefix was hash-guarded against final trace SHA256 `aa51e623cbc20cb5c81eaa4c9a68b9468860bb4c4f73d23b342686fa5bff5c77`. Two-supported-line geometry has |cos|<=1/2, where the separately implemented Bernstein-control and direct-cosine Sturm audit passed153 polynomials with margin1e-7.

At most one opposite pair per primitive line is supported. Four ordered collinear terms and origin cross terms give G(t)=2/q(H^T W H+z^T H), W=sum ckBk/q. With the plane inequalities now verified, positivity of the full finite kernel sum implies sum G>=-s.

## Exact chosen input
The numerical source is in quadratic Bernstein radial coordinates. Independent exact transformation T=[[1,-2,1],[0,2,-2],[0,0,1]] verifies B_monomial=T^T B_Bernstein T and z_monomial=T^T z_Bernstein, followed by subtraction1/10000 from B0[0,0] and1/200000 from z0. Exact Sylvester tests pass for every final block. Origin Schur slack is exactly1e-12. Numerical optimization accuracy is not a premise.

## Scalar Fourier closure
Here n/2-1=21 and contact x=22|v|². The exact Fourier-polynomial recurrence is L0=1, L1=22-x and Lk=((2k+20-x)L(k-1)-(k+20)L(k-2))/k. These constants were independently derived from dimension44; neither the discarded faulty generator nor its output is used. Both Fourier-transform directions are checked exactly for degree71, Q(0)=1, and the stored positive-factor product reconstructs Q exactly.

The compact radial polynomial is h=19/6-(13/132)x, cutoff418/13. The scalar checker covers633 contiguous annuli from22 to180 with cutoff as a boundary. It computes rigorous upper heights for P(x)e^-x+G(x/22), including the full exponential Taylor remainder, and records nonnegative heights even if a minimum-width cell is accepted. Fresh isolated replay passed633 annuli and2535 interval nodes with the same exact budget. All71 roots were freshly enclosed and lie away from the real tail beyond180, with P(180)<0.

For vectors in one mod-2 class, both half-sum and half-difference are lattice vectors. In a squared-radius annulus of width<=4, this implies |cos|<=1-2/T, T the upper squared radius. Below T=2 only one unoriented pair can occur. Otherwise the even polynomial t^(2k)(t²-s²), s=1-2/T, has positive constant harmonic mu_k((2k+1)/(44+2k)-s²) for the tested admissible k, and nonnegative other harmonics. To see nonconstant-harmonic positivity directly, the ratio of the degree-ell normalized harmonic coefficient in t^(2k+2) to that in t^(2k), for even ell<=2k, is (2k+2)(2k+1)/((2k+2-ell)(2k+44+ell)). This is at least the constant ratio(2k+1)/(2k+44), since the denominator changes by ell(2-44-ell)<=0. The new highest-degree coefficient is positive. Thus its diagonal/constant ratio bounds the number of unoriented directions; flooring and minimizing over k=0,...,500 is safe. Both signs and all2^44 classes give the2^45 factor. The exact checker independently recomputes every stored cap and the complete budget.

Consequently, with geometry and fresh scalar replay complete, density <=(22/4)^22/22! (P(0)+s+budget). This rational quantity is independently checked to be strictly below1.7676809433e-6. Unlike43D, this44D certificate has a nonzero residual budget and therefore retains the mod-2 counting argument.

## Audit scope
Exact repair/PSD/Fourier/factor/normalization checks and the cap-polynomial audit are independently implemented. Fresh scalar replay uses an isolated snapshot but shares the producer's reviewed Arb/Taylor implementation. Main replay uses an independent complete-tree coverage parser at384bits and the reviewed shared Jet/Arb interval backend. No independent arithmetic backend is claimed. Complete main geometry has passed. The stored maximum-upper field is a conservative ceiling after prefix resume, not a measured optimized slack.
