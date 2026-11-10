# Comparison sources and conventions

Reviewed 10 November 2026. The numerical baseline is [Henry Cohn's sphere-packing table](https://cohn.mit.edu/sphere-packing/), retrieved through the search index (reported crawl age: four months). Direct live retrieval failed with HTTP 502. Thus these are pinned, source-attributed comparisons, not a claim that no newer results exist.

## Lower bounds

- **21D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Center density 2^(-41)*3^(-23/2)*43^(21/2), as listed in Cohn's general-packing table. Multiply by 4^11*11!*pi^10/22! for density L.

- **19D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Center density 3^(-9) * 5^(-21/2) * 13^(19/2), as listed in Cohn's general-packing table. Multiply by the unit-ball volume 4^10*10!*pi^9/20! for the density L. This is a cited construction, not a new lower bound.

- **20D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Non-lattice periodic packing, five spheres per underlying lattice cell. Exact density L = (pi^10/10!) * 3^20/(2^10 * 5^(21/2)).
- **48D:** Leech and Sloane, [Sphere packings and error-correcting codes](https://doi.org/10.4153/CJM-1971-081-3), Canadian Journal of Mathematics 23 (1971), 718–745. Extremal lattice packing density L = (3*pi/2)^24/24!. The attribution and center density are independently tabulated by Cohn (reference 21). Lattice constructions also supply lower bounds for unrestricted packings; they do not establish optimality among all packings.

## Previous upper bounds

For all four rows, Cohn's table attributes its upper bounds to Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini, [High-dimensional sphere packing and the modular bootstrap](https://arxiv.org/abs/2006.02560v2), JHEP 12 (2020), 066, [DOI](https://doi.org/10.1007/JHEP12(2020)066). The table's upward-rounded entries are 0.0098179551395438 (19D), 0.0071270536033763 (20D), 0.0051596603948176 (21D), and 0.0000004529067791 (48D). These concern unrestricted packing density, not a lattice-only bound. We use these published-table decimals, not the unknown exact optimum of the Cohn–Elkies optimization problem.

## Arithmetic and interpretation

Use the same lower bound L before and after the update, the referenced older upper U, and the repository's displayed certified strict upper N:

- Upper-bound decrease (%) = 100 * (U - N) / U.
- Additive-gap decrease (%) = 100 * ((U - L) - (N - L)) / (U - L) = 100 * (U - N) / (U - L).

The second percentage is the fraction of the earlier uncertainty interval removed. It is not a claim that N lies that percentage above the optimum. Both comparisons use densities, not center densities. N is rounded outward, so we do not take credit for unreported extra digits of the certificate.

`python comparison.py` checks all displayed values, exact-rational upper-percentage arithmetic, gap identities and agreement at 70 and 100 decimal digits. `python comparison.py --write` regenerates README/results tables and the registry. Lower formulas are numerical evaluations of cited exact expressions, not newly proved lower-bound constructions. The readable tables show densities to five significant digits with an approximation mark and percentages to two decimal places. These display values are not substitute rigorous bounds. The machine-readable registry retains its previous precision, including six decimal places for percentages; calculations use the unshortened inputs. `RESULTS.json` also records both additive gaps.
