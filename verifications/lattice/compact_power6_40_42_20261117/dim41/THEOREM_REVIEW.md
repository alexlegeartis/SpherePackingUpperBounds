# Independent audit of the 41D compact harmonic kernel

**Status: COMPLETE INDEPENDENT AUDIT PASS. Lattice-only density upper bound 4.3927514316e-6.**

Only frozen input SHA256 `312bb8e9c7bb637505048a7b28516b48dbc15add44494159552e7bc7e9535680` is covered. This is a dimension-specific optimized kernel, not a dimension restriction of another input.

## Geometric and PSD theorem
Use normalized even Gegenbauer harmonics in dimension41 (parameter39/2), powers6,9,12, support R=19/13 and H=(h^6,h^9,h^12), h(t)=max((19-13t)/6,0). With c0=1, ck>=0, q=sum ck and PSD matrices Bk, the nonzero-vector kernel is K(x,y)=q^-2 sum ck G_(2k)(cos(x,y))H(|x|²)^T Bk H(|y|²). The origin extension K(0,x)=z^T H/q, K(0,0)=s is PSD by its exact constant-harmonic block.

A Gauss-reduced lattice plane has squared norms1<=a<=b and dot0<=d<=a/2. Only the three signed pairs ±v, ±w, ±(w-v) can be supported because all other vectors have squared norm>=2a>=2>R. The normalized three-line cube and outer powers6/feature stride3 are identical to the mathematical formulas reviewed for the43D package, with the41D angular polynomial substituted. Every independent pair belongs to one saturated plane. A complete trace must prove the normalized three-pair sum nonpositive on this cube. Two-supported-line geometry has |cos|<=1/2, where the separately implemented Bernstein-control and direct-cosine Sturm audit passed153 polynomials with margin1e-7.

At most one opposite pair per primitive line is supported. Four ordered collinear terms and origin cross terms give G(t)=2/q(H^T W H+z^T H), W=sum ckBk/q. Since all plane inequalities have passed complete replay, positivity of the full finite kernel sum implies sum G>=-s.

## Exact chosen input
The numerical source is in quadratic Bernstein radial coordinates. Independent exact transformation T=[[1,-2,1],[0,2,-2],[0,0,1]] verifies B_monomial=T^T B_Bernstein T and z_monomial=T^T z_Bernstein, followed by subtraction1/100 from B0[0,0] and1/1000 from z0. Exact Sylvester tests pass for every final block. Origin Schur slack is exactly1e-12. Numerical optimization accuracy is not a premise.

## Scalar Fourier closure
Here n/2-1=39/2 and contact x=17|v|². The exact Fourier-polynomial recurrence is L0=1, L1=41/2-x and Lk=((2k+37/2-x)L(k-1)-(k+37/2)L(k-2))/k. These constants were independently derived from dimension41; neither the discarded faulty generator nor its output is used. Both Fourier-transform directions are checked exactly for degree81, Q(0)=1, and the stored positive-factor product reconstructs Q exactly.

The compact radial polynomial is h=19/6-(13/102)x, cutoff323/13. The scalar checker covers653 contiguous annuli from17 to180 with cutoff as a boundary. It computes rigorous upper heights for P(x)e^-x+G(x/17), including the full exponential Taylor remainder, and records nonnegative heights even if a minimum-width cell is accepted. Fresh isolated replay passed653 annuli and5669 interval nodes with the same exact budget (approximately1.44670074e-05). All81 roots were freshly enclosed and lie away from the real tail beyond180, with P(180)<0.

For vectors in one mod-2 class, both half-sum and half-difference are lattice vectors. In a squared-radius annulus of width<=4, this implies |cos|<=1-2/T, T the upper squared radius. Below T=2 only one unoriented pair can occur. Otherwise the even polynomial t^(2k)(t²-s²), s=1-2/T, has positive constant harmonic mu_k((2k+1)/(41+2*k)-s²) for the tested admissible k, and nonnegative other harmonics. To see nonconstant-harmonic positivity directly, the ratio of the degree-ell normalized harmonic coefficient in t^(2k+2) to that in t^(2k), for even ell<=2k, is (2k+2)(2k+1)/((2k+2-ell)(2k+41+ell)). This is at least the constant ratio(2k+1)/(2k+41), since the denominator changes by ell(2-41-ell)<=0. The new highest-degree coefficient is positive. Thus its diagonal/constant ratio bounds the number of unoriented directions; flooring and minimizing over k=0,...,500 is safe. Both signs and all2^41 classes give the2^42 factor. The exact checker independently recomputes every stored cap and the complete budget.

Consequently, after geometry and fresh scalar replay, density <=(17/4)^(41/2)/Gamma(43/2) (P(0)+s+budget). The independent checker squares this expression, uses Gamma(43/2)=42! sqrt(pi)/(4^21 21!), and proves the strict bound4.3927514316e-6 against a rational Machin lower bound for pi. Unlike43D, this41D certificate has a nonzero residual budget and therefore retains the mod-2 counting argument.

## Audit scope
Exact repair/PSD/Fourier/factor/normalization checks and the cap-polynomial audit are independently implemented. Fresh scalar replay uses an isolated snapshot but shares the producer's reviewed Arb/Taylor implementation. Main replay uses an independent complete-tree coverage parser at384bits and the reviewed shared Jet/Arb interval backend. No independent arithmetic backend is claimed. Complete main geometry has passed the full independent replay.

The numerical origin PSD constraint anticipated repairs1/100 and1/1000, as explained in PROOF_AWARE_OPTIMIZATION.md. This is only a search strategy: the final rational matrices, origin Schur slack and full function inequalities are checked independently. Only the stated immutable input is covered.

## Additional exact tail proof

`check_tail_exact.py` independently computes every rational coefficient of P(180+t) by Horner composition and verifies that all coefficients are strictly negative. Thus P(x)<0 for every x>=180 directly, without any real-root isolation. This supplements the already-passed81-root enclosure check; it does not change the input polynomial or bound.

## Complete main geometry replay

Independent 384-bit replay passed all 184501 nodes and 92251 accepted leaves with empty coverage stack. Trace SHA256 `86ec39b414fcd0ab1866329cf732eeccd6a706a85169ff31bf6b71f33e2d1d97`. The resumed prefix ceiling -1e-7 is conservative, not a measurement of the worst leaf.

## Complete main geometry replay

Independent 384-bit replay passed all 184501 nodes and 92251 accepted leaves with empty coverage stack. Trace SHA256 `86ec39b414fcd0ab1866329cf732eeccd6a706a85169ff31bf6b71f33e2d1d97`. The resumed prefix ceiling -1e-7 is conservative, not a measurement of the worst leaf.
