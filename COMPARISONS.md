# Comparison sources and conventions

Reviewed 10 November 2026. Except for the updated 4D and 5D rows below, the numerical baseline is [Henry Cohn's sphere-packing table](https://cohn.mit.edu/sphere-packing/), retrieved through the search index (reported crawl age: four months). Direct live retrieval failed with HTTP 502. Thus these are pinned, source-attributed comparisons, not a claim that no newer results exist.

## Lower bounds

- **4D:** the D4 lattice has density pi^2/16, as tabulated by [Cohn](https://cohn.mit.edu/sphere-packing/). This is the known construction, not an optimality claim for unrestricted packings.

- **21D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Center density 2^(-41)*3^(-23/2)*43^(21/2), as listed in Cohn's general-packing table. Multiply by 4^11*11!*pi^10/22! for density L.

- **19D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Center density 3^(-9) * 5^(-21/2) * 13^(19/2), as listed in Cohn's general-packing table. Multiply by the unit-ball volume 4^10*10!*pi^9/20! for the density L. This is a cited construction, not a new lower bound.

- **20D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Non-lattice periodic packing, five spheres per underlying lattice cell. Exact density L = (pi^10/10!) * 3^20/(2^10 * 5^(21/2)).
- **48D:** Leech and Sloane, [Sphere packings and error-correcting codes](https://doi.org/10.4153/CJM-1971-081-3), Canadian Journal of Mathematics 23 (1971), 718–745. Extremal lattice packing density L = (3*pi/2)^24/24!. The attribution and center density are independently tabulated by Cohn (reference 21). Lattice constructions also supply lower bounds for unrestricted packings; they do not establish optimality among all packings.

## Previous upper bounds

**4D current comparison:** [arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), version dated 8 December 2025, reports 0.63610277 for unrestricted packing density. Checked directly on 10 November 2026. This is the paper’s printed value at eight decimal places; no additional digits, directed-rounding guarantee or independent certificate replay are asserted. Percentages therefore compare against that reported decimal, not a newly certified exact baseline. The repository bound <0.635766243 remains smaller. This paper reports dimensions 4–7 and 9–16 in Table 1, with a separate 8D case study; it provides no replacement comparison for our 19D, 20D, 21D or 48D rows.

**4D historical comparison:** Cohn–de Laat–Salmon, [Three-point bounds for sphere packing](https://arxiv.org/abs/2206.15373), supplies the earlier upper 0.6361073321551329 from Cohn's table. The new package also assumes four published near-neighbor average-degree bounds from [de Laat–de Oliveira Filho–Vallentin, Table 2](https://arxiv.org/pdf/1206.2608). See the 4D proof for exact thresholds and the limits of its independent computational audit.

For dimensions 19, 20, 21 and 48, Cohn's table attributes its upper bounds to Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini, [High-dimensional sphere packing and the modular bootstrap](https://arxiv.org/abs/2006.02560v2), JHEP 12 (2020), 066, [DOI](https://doi.org/10.1007/JHEP12(2020)066). The table's upward-rounded entries are 0.0098179551395438 (19D), 0.0071270536033763 (20D), 0.0051596603948176 (21D), and 0.0000004529067791 (48D). These concern unrestricted packing density, not a lattice-only bound. We use these published-table decimals, not the unknown exact optimum of the Cohn–Elkies optimization problem.

## Arithmetic and interpretation

Use the same lower bound L before and after the update, the referenced older upper U, and the repository's displayed certified strict upper N:

- Upper-bound decrease (%) = 100 * (U - N) / U.
- Additive-gap decrease (%) = 100 * ((U - L) - (N - L)) / (U - L) = 100 * (U - N) / (U - L).

The second percentage is the fraction of the earlier uncertainty interval removed. It is not a claim that N lies that percentage above the optimum. Both comparisons use densities, not center densities. N is rounded outward, so we do not take credit for unreported extra digits of the certificate.

`python comparison.py` checks all displayed values, exact-rational upper-percentage arithmetic, gap identities and agreement at 70 and 100 decimal digits. `python comparison.py --write` regenerates README/results tables and the registry. Lower formulas are numerical evaluations of cited exact expressions, not newly proved lower-bound constructions. The readable tables show densities to five significant digits with an approximation mark and percentages to two decimal places. These display values are not substitute rigorous bounds. The machine-readable registry retains its previous precision, including six decimal places for percentages; calculations use the unshortened inputs. `RESULTS.json` also records both additive gaps.

## Dimension5 addition

The lower density is pi²/(15sqrt(2)), from the standard D5 lattice, tabulated by [Cohn](https://cohn.mit.edu/sphere-packing/). The current comparison uses the paper-reported0.5122645 from [arXiv:2512.04829v2 Table1](https://arxiv.org/html/2512.04829v2#S3.T1), not the older0.5126451306253027 three-point benchmark. We have not independently replayed the newer paper's own numerical certificates. The decimal is a pinned reported comparison, not an exact optimum claim.

The new5D certificate uses the rigorous neighbor-count results in [de Laat–de Oliveira Filho–Vallentin, Table2](https://arxiv.org/html/1206.2608#S6.T2) as external mathematical inputs. Their use is audited; those original spherical SDP certificates are not reconstructed here. The new witness is verified over complete continuous domains. Its audit is source review plus independent reruns and objective reconstruction, not a clean-room implementation or formal proof.

## Four-dimensional refinement

The previous p4 repository witness has outward upper 0.634526023. The earlier audited repository certificate 0.635766243 is preserved in the parent package. Relative to the updated paper-reported 0.63610277 baseline and D4 lower pi²/16, the conservative rounded comparison removes 8.189832% of the additive gap and lowers the upper by 0.247876%. The same four externally published neighbor-count inputs are required; audit consists of source review, isolated replay and independent objective reconstruction, not a clean-room implementation or formal verification.

The earlier 5D refinement uses the same external neighbor-count inputs and yields the outward-rounded bound0.510584101, closing3.5747932% of the gap to the D5 lower density relative to the pinned paper-reported0.5122645. Audit scope and benchmark caveats above remain unchanged.

The earlier radial-degree-six 4D witness strengthens the outward bound to **0.633430327**. Using this rounded value and the same updated paper-reported baseline, the additive-gap decrease is 13.881022%. Previous p3 and p4 certificates are preserved. The external count assumptions and audit scope are unchanged; the new far cosine cap 1/2 is proved by a radial-corner convexity argument in the refinement proof.

The earlier radial-degree-four 5D witness strengthens the outward bound to **0.510522525**, removing 3.705787% of the additive gap against the same paper-reported 0.5122645 baseline. Both previous 5D certificates are preserved. Neither the current 4D nor 5D result reaches a 20% gap decrease.

## Origin-coordinate refinements

The current registered bounds are 4D <0.633428734445, 5D <0.510501094478, 19D <0.009618343139, 20D <0.006979497272, 21D <0.005149549312. They result from exact Schur-complement slack removal in the existing PSD kernels, not new dimensions or a new theorem. Full source-domain replays and independent exact factorization and Arb objective reconstruction passed. Comparators and external mathematical inputs are unchanged; the generated tables recompute percentages from the strict outward decimals. See [proof and replay](verifications/origin_slack/).
