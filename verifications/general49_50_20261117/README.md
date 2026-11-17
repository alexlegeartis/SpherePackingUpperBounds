# Certified unrestricted sphere-packing upper bounds in49D and50D

-49D: density **<2.85041752e-7**.
-50D: density **<2.04993038e-7**.

Method: exact positive-factor Fourier witnesses, dimension-adapted cyclic positive-semidefinite angular kernels, and universal thin-annulus neighbor budgets tightened by exact Levenshtein-type spherical-code polynomials. These are unrestricted packing bounds, not lattice-only results. No mod2 class or coset assumption enters the proof.

## Reproduction

Install numpy, mpmath and python-flint, then run `python -B run_all.py` without optimization. It replays Fourier positivity/transform/tail, the exact angular structure, full compact interval gate, every radial annulus and the independent rational cap/density checks for both dimensions. `--quick` checks exact algebra and caps against the saved interval receipts only; it does not freshly establish interval signs. The complete independent interval runs use1024-bit Arb and are recorded here. This review and high-precision replay share FLINT/Arb with the producer, not an independent arithmetic backend.

The dim49/dim50 polynomial, angular and radial inputs are unchanged from the independently audited preceding unrestricted proof packages. The saved caps replace only the old larger neighbor-count allowances. Independent check_caps.py reconstructs normalized Gegenbauer expansions by recurrence and descending elimination rather than the producer's closed coefficient formula. check_exact49.py is the prior independent49/51 checker scoped explicitly to49; no51 certificate is included. The50 independent Gaussian transform checker is preserved unchanged.

## Complete proof chain

Each dim*/THEOREM.md gives the cyclic PSD plus per-center error-budget theorem. Each radial interval bounds the corrected real-space expression by a nonnegative height. All632 bins partition x in[22,180], where x=22r²; the global tail handles the remainder. When the interval target is not met, the actual outward upper bound is charged, never an unproved target.

For a squared-radius bin[L,R] of width at most1, all neighboring centers have angular products at mosts=1-1/(2R). Every new cap is exactly f(t)=(t-s)K(t)^2 or(1+t)(t-s)K(t)^2, is nonpositive on[-1,s], and has nonnegative normalized Gegenbauer coefficients with f0>0. The spherical addition theorem gives cardinality at mostfloor(f(1)/f0). Cap times the independently replayed nonnegative height is charged per center, with no lattice factor. The exact rational source polynomials and cap coefficients are in caps49/ and caps50/.

The origin Schur complement includes the actual1e-12 guard. Density normalization uses exact factorial arithmetic, with a rational Machin lower pi enclosure for odd49D. Both safe decimals are verified by strict rational inequalities.

## Comparison scope

The attributable external comparison is Afkhami-Jeddi,Cohn,Hartman,de Laat,Tajdini(2020), numerical spinless bootstrap, using the highest supplied truncation in their MIT dataset. These are numerical benchmarks, not certified packing bounds; the general same-radius Cohn–Elkies lift is conjectural and functional positivity is not supplied by the spectra. See [the complete49–96 benchmark table](../../NUMERICAL_BENCHMARKS49_96.md).

Our safe decimals are9.7387573447% and6.8379867092% below those numerical benchmarks in49D and50D. This is not an exhaustive priority claim, nor a claim that the external estimates have a rigorous error certificate. No51D or later new bound is included.
