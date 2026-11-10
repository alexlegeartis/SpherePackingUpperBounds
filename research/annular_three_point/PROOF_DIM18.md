# Dimension 18: locally verified continuum certificate

## Claim and status

The exact witness `witnesses_global_fine31/dim18.json` yields the unrestricted sphere-packing density bound

    Delta_18 < 0.013408931.

All local exact/full-domain computational checks pass. Independent review passed; see ../../verifications/dim18/REVIEW.md. A worldwide-priority claim is not made. No worldwide-priority claim is made.

The density enclosure is

    0.01340893006264654494029079322832899292797624195812881555218302593199193264424

with Arb radius below 3.74e-78. The new bound is about 0.5666% below the pinned Cohn table comparator 0.0134853404450862; this is a pinned-source comparison, not a best-known assertion.

## Framework and exact witness

Apply Cohn–de Laat–Salmon, *Three-point bounds for sphere packing*, Theorem 1.1, https://arxiv.org/abs/2206.15373, as in the existing independently audited dimension-20 repository proof (`reference/PROOF.md`).

Set n=18, r²=9/4, and R²=26999999999999997/10000000000000000, exactly the finite decimal retained in the witness. The slight decimal deviation from 27/10 is intentional provenance preservation; verification uses the recorded rational value, not 27/10. Use five equally spaced radial hats, angular degrees 0 through 16, and radial polynomial degree 31.

The witness defines rational polynomials P,H and rational symmetric matrices A_0,...,A_16. Put N=H(0)>0 and

    f2(x) = P(pi |x|²) exp(-pi |x|²) / N.

The exact Fourier transform is H(pi |xi|²) exp(-pi |xi|²)/N. The normalized Laguerre eigenbasis is L_k^8(2X)/binomial(k+8,k), with Fourier eigenvalue (-1)^k. Rational reconstruction checks the identity exactly, rather than trusting optimizer coefficients.

At the origin the degree-zero feature is (1,0,...,0); in the annulus it is (0,w(|x|²)); elsewhere it is zero. The degree-zero kernel uses A_0. The positive angular-degree terms on the annulus are

    sum_l=1^16 [w(|x|²)^T A_l w(|y|²)] P_l^(18)(xhat·yhat),

with P_l^(18)(1)=1. Divide the full kernel by N. Exact positive LDL pivots prove all 17 blocks positive definite. The standard Gegenbauer addition theorem then proves global kernel positivity. The continuum theorem permits this isolated-origin, discontinuous kernel.

## Full-domain checks

1. Exact Sturm calculations prove H>0 on [0,infinity).
2. A rational lower bound for pi, itself enclosed by Arb, and exact Sturm calculations prove P(pi u)<0 for u>=R².
3. Each coarse radial cell is subdivided into eight pieces, making 33 fine nodes. For each of the 561 unordered fine-node pairs, a rational upper cosine bound is computed at the upper radial neighbors. Its squared inequality is checked exactly. The monotonicity hypothesis R²<2r² is checked.
4. Every corresponding rational angular polynomial is strictly negative throughout [-1,c], using exact root counts and endpoint signs. Bilinear interpolation extends this to the complete admissible annular pair domain, not just nodes or sampled angles.
5. On every coarse radial cell, 2f3(x,0)+f3(x,x) is an exact quadratic in |x|² divided by N. The verifier encloses the coupled annular expression with 256-bit Arb Taylor bounds, covering 64 boxes with maximum subdivision depth 6. Every terminal interval is strictly negative.

These checks cover the computational hypotheses of the cited theorem. The theorem and Gegenbauer addition theorem are mathematical inputs, not formalized in this package.

## Objective and reproduction

The exact normalized objective is

    777222638133492104364138274745 / 26843545600001902746210926592.

Multiply by pi^9/9! * (9/16)^9. The interval result is strictly below 0.013408931.

Run, without Python -O:

    python run_certificate.py witnesses_global_fine31/dim18.json

This invokes both `verify_global.py` and `verify_annular.py` and fails if either fails. `audit_global_fine31_18.log` records the complete replay. Numerical discovery used `search_global_fine.py --degree 31 --dimensions 17 18`; solver status was optimal_inaccurate for dimension 18 and is not proof evidence. The exact witness, not that numerical status, defines the claimed bound.
