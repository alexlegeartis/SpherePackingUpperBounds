# Comparison sources and conventions

Reviewed 10 November 2026. The numerical baseline is [Henry Cohn's sphere-packing table](https://cohn.mit.edu/sphere-packing/), retrieved through the search index (reported crawl age: four months). Direct live retrieval failed with HTTP 502. Thus these are pinned, source-attributed comparisons, not a claim that no newer results exist.

## Lower bounds

- **20D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Non-lattice periodic packing, five spheres per underlying lattice cell. Exact density L = (pi^10/10!) * 3^20/(2^10 * 5^(21/2)).
- **48D:** Leech and Sloane, [Sphere packings and error-correcting codes](https://doi.org/10.4153/CJM-1971-081-3), Canadian Journal of Mathematics 23 (1971), 718–745. Extremal lattice packing density L = (3*pi/2)^24/24!. The attribution and center density are independently tabulated by Cohn (reference 21). Lattice constructions also supply lower bounds for unrestricted packings; they do not establish optimality among all packings.

## Previous upper bounds

For both rows, Cohn's table attributes its upper bounds to Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini, [High-dimensional sphere packing and the modular bootstrap](https://arxiv.org/abs/2006.02560v2), JHEP 12 (2020), 066, [DOI](https://doi.org/10.1007/JHEP12(2020)066). The table's upward-rounded entries are 0.0071270536033763 (20D) and 0.0000004529067791 (48D). These concern unrestricted packing density, not a lattice-only bound. We use these published-table decimals, not the unknown exact optimum of the Cohn–Elkies optimization problem.

## Arithmetic and interpretation

Use the same lower bound L before and after the update, the referenced older upper U, and the repository's displayed certified strict upper N:

- Upper-bound decrease (%) = 100 * (U - N) / U.
- Additive-gap decrease (%) = 100 * ((U - L) - (N - L)) / (U - L) = 100 * (U - N) / (U - L).

The second percentage is the fraction of the earlier uncertainty interval removed. It is not a claim that N lies that percentage above the optimum. Both comparisons use densities, not center densities. N is rounded outward, so we do not take credit for unreported extra digits of the certificate.

`python comparison.py` checks all displayed values, exact-rational upper-percentage arithmetic, gap identities and agreement at 70 and 100 decimal digits. `python comparison.py --write` regenerates README/results tables and the registry. Lower formulas are numerical evaluations of cited exact expressions, not newly proved lower-bound constructions. Percentages are rounded to six decimal places. `RESULTS.json` also records both additive gaps.
