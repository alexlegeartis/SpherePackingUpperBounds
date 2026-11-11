# Independent review of the 22D and 23D Gaussian cyclic certificates

**PASS.** The reviewed exact witnesses prove unrestricted sphere-packing density

- Delta_22 < **0.003653383161**.
- Delta_23 < **0.002642540089**.

No mathematical defect was found. One documentation correction was requested and made: each angular polynomial has actual degree **14**, padded to 20 coefficients. The producer's degree-19 Bernstein calculation is valid degree elevation, not evidence that Q has degree 19. No mathematical input changed.

## Independent computations

This audit imports no producer verifier. It includes exact copies of the two radial certificates, two angular polynomials, and two comparator records. `source_hashes.json` binds these copies to their original paths. `verify_all.py` checks every copy's SHA256 and compares the originals byte-for-byte when those original paths exist.

1. **Fourier identity:** reconstructed all 112 coefficients in each direction by the three-term recurrence for k! L_k^(n/2-1), independently of the producer's coefficientwise rising-factorial formula. The generating-function identity FT[e^(-(1-z)X)] = (1-z)^(-n/2)e^(-Y/(1-z)) proves the normalization for exp(-2*pi*i*x.y).
2. **Global radial signs:** complete certified FLINT complex-root isolation covers all 111 roots of P and H in each dimension. Each disk either avoids the real axis or lies wholly below the relevant threshold. Endpoint signs then give H>0 on [0,infinity) and P<0 on [7*pi_lower,infinity). This independently reruns the same trusted root-isolation backend used by the producer; it is not an independent implementation of that backend. Exact-rational Sturm attempts on the degree-111 polynomials were stopped for cost before completing and are not counted as a passing check. Their logs are retained as incomplete attempts.
3. **Angular PSD:** explicit normalized Gegenbauer coefficient sums reconstruct each Q in its actual physical dimension. All harmonic coefficients are nonnegative. An independent even-moment calculation gives b0=1. With B=Q(1), the origin degree-zero Gram entries are MB, -M, M/B: positive diagonal and determinant zero.
4. **Pair inequality:** a second-derivative Taylor bound for the entire divided pair function G(X)=P(X)e^((kappa-1)X-C)-2M+Me^(-kappa X+C), C=51*pi/4, proves G<0 on the full interval [s*pi_lower,7*pi_upper]. The independent covers use **95 intervals in 22D** and **93 in 23D**, maximum depth 11. This differs from the producer's bound on P followed by monotone exponential factors. The global physical sign and 0<h<=1 finish the unbounded radial ray.
5. **Angular signs:** exact rational Sturm sequences prove Q<-1/100000 throughout [-1,1/2], Q<0 through 5001/10000, and exactly one subsequent root before 1. The sign variations are [7,7,6] in both dimensions.
6. **Compact cyclic domain:** a new parameterization v=2ut covers the complete rectangle 1<=u<=5, 0<=v<=1. It requires **no excluded rectangles**. With t=v/(2u), w=sqrt(1+u^2-v), the three terms are Q(t), exp[-40(1-v)]Q((1-v/2)/w), and exp[-40(u^2-v)]Q((u-t)/w). Exact degree-14 Bernstein bounds and Arb derivative enclosures certify **1,368 rectangles in 22D** and **1,415 in 23D**, maximum depth 15. This is a different cover and implementation from the producer's sorted (t,u) rectangles.
7. **Normalization and rounding:** s=n/8, kappa=102/n, R^2=7 and kappa*s=51/4 are checked exactly. The objective (P(0)+MB)/H(0) matches the certificate exactly. Rational alternating-series Machin bounds prove 51*pi/4>40. Rational outward rounding proves both strict decimals, using a squared exact inequality in odd dimension rather than a floating-point square-root approximation.

The resulting density enclosures begin 0.0036533831609501341267 (22D) and 0.0026425400886950499942 (23D), agreeing with the producer. Full intervals and partition hashes are in the JSON results.

## Mathematical review of the complete argument

The addition theorem supplies a finite PSD feature representation. A smooth radial cutoff, zero near the origin and one outside the excluded ball, extends the specified origin and nonorigin features to globally smooth Schwartz features. Gaussian h is nonzero at arbitrarily large radius; R^2=7 is only a physical-polynomial sign threshold. No compact-support or finite-neighbor-count assumption is used.

For every packing triangle, choose adjacent shortest and second-shortest sides, scale the shortest to 1, and write the second as u>=1. The remaining longest side imposes 2ut<=1 and t<=1/2. Dividing its cyclic kernel sum by the positive factor corresponding to the two shortest sides gives exactly the expression checked above, at A=kappa*pi*(actual shortest side)^2>=51*pi/4>40. The common M/B and Gaussian normalization cancel with the correct signs.

To extend A=40 to all larger A, put p=1-2ut and q=u^2-2ut, so 0<=p<=q and z1<=z2. The verified single sign change of Q excludes a positive first coefficient followed by a negative second coefficient. Two nonnegative coefficients decrease with A; two nonpositive ones have sum at most zero. In the mixed case a<=0<=b, the derivative of a e^(-pA)+b e^(-qA) changes sign at most once, from negative to positive. Its maximum is at A=40 or its nonpositive infinite limit. The boundary cases p=0 and p=q are also valid. Since Q(t)<0, compact negativity persists.

The omitted unbounded domains are rigorously controlled: for t<=0 both exponent gaps are at least 1 and |Q|<=B by positive harmonic coefficients. For t>=0 and u>=5, w>=u and z1<=1/u<=1/5, so the first secondary term is negative, while the other gap is at least u^2-1>=24. The exact check 40^30/30!>2B*100000 proves that all possibly positive tail terms are smaller than the cap margin. This covers degenerate and unbounded triangles without sampling.

For a periodic packing, summing PSD inequalities over representative centers produces NK(0,0), the ordered nonidentical-pair sum of 2K(x,0)+K(x,x), and twice each unordered triangle's cyclic sum, modulo translations. The triple term is nonpositive, so the pair sum is at least -NK(0,0). Combining f+2K(x,0)+K(x,x)<=0 gives the upper bound N(f(0)+K(0,0)) for the complete pair sum of f. Poisson summation with H>=0 gives its lower bound N^2 H(0)/V. Thus N/V<=(P(0)+MB)/H(0). Absolute convergence justifies regrouping, and finite-window periodization with a boundary gap extends the density bound to arbitrary packings.

The packing-ball volume is pi^(n/2)*(s/4)^(n/2)/Gamma(n/2+1). In odd dimension n=2m+1 it simplifies to 2*pi^m*s^(m+1/2)*(m+1)!/(2m+2)!, exactly as used in the independent squared rounding check. This is a direct cyclic periodic proof; the individual-sign hypothesis of printed Theorem 1.1 is not assumed.

## Comparison scope and replay

The pinned older comparisons are 0.0037259419689206 in 22D and 0.0026842798864291 in 23D. Their records identify https://cohn.mit.edu/sphere-packing/ and arXiv:2006.02560v2, with an indexed retrieval dated 2026-11-11, four-month-old crawl, and a 502 from the live page. This audit verifies the exact comparison arithmetic, not a fresh retrieval or exhaustive worldwide priority search. It establishes improvements of more than 1.947% and 1.554%, respectively, against those specific pinned numbers.

From the workspace root, run:

    research/density_solid_cyclic_20261111/.venv/bin/python research/dim22_23_gaussian_independent_20261111/verify_all.py

Alternatively use any Python environment with python-flint==0.9.0. Assertions must remain enabled. The copied inputs and these three scripts are sufficient; the original producer directories are optional. The final result is `audit.json`, with detailed per-dimension angular/radial JSON files. This is independent computer-assisted review, not external peer review or proof-assistant formalization.
