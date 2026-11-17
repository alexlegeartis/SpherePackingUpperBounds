# Independent audit of the 42D compact harmonic kernel

**Status: COMPLETE PASS. Main geometry, exact algebra, caps, fresh scalar replay, tail and normalization have passed the described independent audit.**

Only frozen input SHA256 `bd730f1b4d529af26ce3e33c9cb470ac050b16c02b06db67a2ed8838c8235d82` is covered. This is a dimension-specific optimized kernel, not a dimension restriction of another input.

## Geometric and PSD theorem
Use normalized even Gegenbauer harmonics in dimension42 (parameter20), powers6,9,12, support R=19/13 and H=(h^6,h^9,h^12), h(t)=max((19-13t)/6,0). With c0=1, ck>=0, q=sum ck and PSD matrices Bk, the nonzero-vector kernel is K(x,y)=q^-2 sum ck G_(2k)(cos(x,y))H(|x|²)^T Bk H(|y|²). The origin extension K(0,x)=z^T H/q, K(0,0)=s is PSD by its exact constant-harmonic block.

A Gauss-reduced lattice plane has squared norms1<=a<=b and dot0<=d<=a/2. Only the three signed pairs ±v, ±w, ±(w-v) can be supported because all other vectors have squared norm>=2a>=2>R. The normalized three-line cube and outer powers6/feature stride3 are identical to the mathematical formulas reviewed for the43D package, with the42D angular polynomial substituted. Every independent pair belongs to one saturated plane. The complete trace has independently proved the normalized three-pair sum nonpositive on this cube:203647 nodes,101824 accepted leaves and an empty final stack at384bits. Every prefix was hash-guarded against final trace SHA256 `e3a784900a2042c3c0694b592dfeb9579b55e989ea389a7bf2231cdc340dab03`. Two-supported-line geometry has |cos|<=1/2, where the separately implemented Bernstein-control and direct-cosine Sturm audit passed153 polynomials with margin1e-7.

At most one opposite pair per primitive line is supported. Four ordered collinear terms and origin cross terms give G(t)=2/q(H^T W H+z^T H), W=sum ckBk/q. With the plane inequalities now verified, positivity of the full finite kernel sum implies sum G>=-s.

## Exact chosen input
The numerical source is in quadratic Bernstein radial coordinates. Independent exact transformation T=[[1,-2,1],[0,2,-2],[0,0,1]] verifies B_monomial=T^T B_Bernstein T and z_monomial=T^T z_Bernstein, followed by subtraction1/500 from B0[0,0] and1/1000 from z0. Exact Sylvester tests pass for every final block. Origin Schur slack is exactly1e-12. Numerical optimization accuracy is not a premise.

## Scalar Fourier closure
Here n/2-1=20 and contact x=17|v|². The exact Fourier-polynomial recurrence is L0=1, L1=21-x and Lk=((2k+19-x)L(k-1)-(k+19)L(k-2))/k. These constants were independently derived from dimension42; neither the discarded faulty generator nor its output is used. Both Fourier-transform directions are checked exactly for degree81, Q(0)=1, and the stored positive-factor product reconstructs Q exactly.

The compact radial polynomial is h=19/6-(13/102)x, cutoff323/13. The scalar checker covers653 contiguous annuli from17 to180 with cutoff as a boundary. It computes rigorous upper heights for P(x)e^-x+G(x/17), including the full exponential Taylor remainder, and records nonnegative heights even if a minimum-width cell is accepted. Fresh isolated replay passed653 annuli and6471 interval nodes with the same exact budget (approximately7.17543e-6). All81 roots were freshly enclosed and lie away from the real tail beyond180, with P(180)<0.

For vectors in one mod-2 class, both half-sum and half-difference are lattice vectors. In a squared-radius annulus of width<=4, this implies |cos|<=1-2/T, T the upper squared radius. Below T=2 only one unoriented pair can occur. Otherwise the even polynomial t^(2k)(t²-s²), s=1-2/T, has positive constant harmonic mu_k((2k+1)/(42+2k)-s²) for the tested admissible k, and nonnegative other harmonics. To see nonconstant-harmonic positivity directly, the ratio of the degree-ell normalized harmonic coefficient in t^(2k+2) to that in t^(2k), for even ell<=2k, is (2k+2)(2k+1)/((2k+2-ell)(2k+42+ell)). This is at least the constant ratio(2k+1)/(2k+42), since the denominator changes by ell(2-42-ell)<=0. The new highest-degree coefficient is positive. Thus its diagonal/constant ratio bounds the number of unoriented directions; flooring and minimizing over k=0,...,500 is safe. Both signs and all2^42 classes give the2^43 factor. The exact checker independently recomputes every stored cap and the complete budget.

Consequently, with geometry and fresh scalar replay complete, density <=(17/4)^21/21! (P(0)+s+budget). This rational quantity is independently checked to be strictly below3.473037643e-6. Unlike43D, this42D certificate has a nonzero residual budget and therefore retains the mod-2 counting argument.

## Audit scope
Exact repair/PSD/Fourier/factor/normalization checks and the cap-polynomial audit are independently implemented. Fresh scalar replay uses an isolated snapshot but shares the producer's reviewed Arb/Taylor implementation. Main replay uses an independent complete-tree coverage parser at384bits and the reviewed shared Jet/Arb interval backend. No independent arithmetic backend is claimed. Complete main geometry has passed. The stored maximum-upper field is a conservative ceiling after prefix resume, not a measured optimized slack.

Only the final refit with the stated frozen hash is audited. The earlier rejected42D candidate is explicitly excluded.

## Additional exact tail proof

`check_tail_exact.py` independently computes every rational coefficient of P(180+t) by Horner composition and verifies that all coefficients are strictly negative. Thus P(x)<0 for every x>=180 directly, without any real-root isolation. This supplements the already-passed81-root enclosure check; it does not change the input polynomial or bound.
