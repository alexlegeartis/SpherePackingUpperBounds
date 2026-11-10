# Independent audit of the fixed 21-dimensional certificate

Verdict: PASS. The rational witness snapshot21_final.json.exact.json proves an unrestricted sphere-packing density upper bound strictly below 0.005149554, under Cohn–de Laat–Salmon, Three-point bounds for sphere packing, Theorem 1.1 (https://arxiv.org/abs/2206.15373). This is not merely a lattice bound. No optimality or priority claim is made.

The independent executable is independent21_audit.py, with exact results in independent21_audit.json. It reads the rational witness directly and does not use the discovery solver or producer verification routines. Run from the workspace root with Python and python-flint 0.9.0 available. Solver status AlmostSolved has no role in the proof.

## Exact algebra and global signs

Here n=21, r²=21/8, R²=63/20, the radial degree is 39, and angular degrees are 0 through 16. All 17 symmetric rational matrices are strictly positive definite: the independent Schur calculation found 86 positive rational pivots. The degree-zero block has an additional origin coordinate.

Writing f2(x)=exp(-π|x|²)P(π|x|²)/N, the Fourier polynomial H was reconstructed directly from monomials, not from the stored Laguerre coefficients. For j=0,...,39 its coefficient is

H_j = sum_(k>=j) p_k k! (-1)^j binom(k+19/2,k-j)/j!.

Every coefficient matches the witness exactly; N=H(0)>0. This checks the half-integer dimension parameter independently. Exact rational Sturm chains prove H>0 on [0,infinity) and P<0 beyond R² times a rational lower bound for π. Forty-term alternating Machin bounds supply that lower bound, so the physical sign holds on the required ray |x|²>=R².

## Entire angular domain, not just samples

The angular polynomials were independently built using the explicit coefficient formula for normalized Gegenbauer polynomials with lambda=19/2, rather than the producer recurrence. All 561 unordered pairs of 33 refined squared-radius nodes pass exact Sturm negativity on [-1,cap]. Each cap is a rational upper bound obtained by integer square-root arithmetic for the allowed cosine at both upper neighboring nodes.

Within every original radial cell the weights are affine in squared radius, and the kernel is bilinear in the two weight vectors. Thus its value is a convex combination of four refined-corner kernels at the same angle. The allowed cosine (u+v-r²)/(2 sqrt(uv)) is increasing in each coordinate on the whole box: its derivative numerator is bounded below by 2r²-R²>0. Consequently the verified upper-neighbor caps cover every angle required at all four corners of every refined cell. Symmetry covers the remaining ordered pairs. Endpoint values and coincident cell boundaries are included.

## Coupled annular inequality and theorem

The coupled kernel term is reconstructed as a rational quadratic in squared radius on each of the four original cells, including both origin cross terms and all angular diagonals P_l(1)=1. Independent 384-bit Arb evaluation with a second-order Taylor remainder proves physical function plus this quadratic strictly negative throughout the annulus. The exhaustive binary cover has 50 certified terminal intervals, maximum depth 7.

The radial function is Schwartz. The kernel is positive definite globally by the spherical addition theorem and the matrix positivity, with features zero off the closed annulus and the separate degree-zero feature at the origin. Theorem 1.1 requires kernel positive definiteness, not continuity of these features. The Fourier normalization, exterior sign, full off-diagonal annular sign and coupled inequality therefore establish every hypothesis. Both f2 and f3 are divided by the same N.

## Objective and odd-dimensional volume

The exact objective is

U = 2346843009126750540349/76293945311687615555.

The density is U vol(B_(r/2)^21). Using Gamma(23/2)=22! sqrt(pi)/(4^11 11!), this is

U * 2*11!/22! * pi^10 * (21/8)^(21/2).

The independent script squares this positive expression, substitutes a rational Machin upper bound for pi, and proves the resulting rational strictly less than (5149554/10^9)^2. Thus the reported decimal upper bound does not depend on floating-point rounding of a half-integer power or Gamma function. A numerical evaluation is 0.0051495531437087255.

The saved source small_density/baseline.html lists upper comparator 0.0051596603948176 and lower density 0.0029298291067157 for dimension 21. Against those figures, the conservative rounded new upper improves the upper bound by about 0.195873% and closes about 0.453236% of the additive gap. This is not a 20% improvement. Current priority and a fresh literature search have not been verified here.

## Packaging cautions

The producer verifier prints several check results without asserting them before writing a certificate. This is harmless for this fixed witness because the independent verifier asserts every condition, but should be hardened before general reuse. Its JSON scope text still says annular verification pending; the independent annular proof now resolves that issue. Preserve the fixed witness rather than substituting later solver output under the same name.
