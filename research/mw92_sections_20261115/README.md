# New best-known lower bounds in six previously missing dimensions

**86D–91D now have explicit improved constructions.** These add six dimensions to the dated literature review, bringing the total to **39 of the 49 dimensions from 48D through 96D**. All six objects are explicit periodic packings; their period bases and rational shifts can be exported from the certificate.

All densities below are **center densities**. “Previous” means the exact repository checkpoint `a30feca`. External comparisons use the named published work or its stated classical consequences.

| Dimension | New center density | Gain over preceding bound | Gain over named external comparison | Translates |
| --- | ---: | ---: | ---: | ---: |
| 86 | 3.532468471752e+12 | **39.116648%** | **39.116648%** | 11 |
| 87 | 1.220042544453e+13 | **96.155351%** | **96.155351%** | 6 |
| 88 | 5.045199642804e+13 | **186.786397%** | **186.786397%** | 5 |
| 89 | 2.039717945339e+14 | **185.488837%** | **>185.413282%** | 4 |
| 90 | 7.844458559109e+14 | **117.127769%** | **>117.089650%** | 3 |
| 91 | 2.478889293267e+15 | **33.154794%** | **>33.143162%** | 2 |

The external comparisons are **BCS1982 sections in 86D–87D**, **BCS1982 in 88D** (also the density of Sakamoto’s 2026 minimum-8 unimodular lattices), and the **full Rogers-simplex section consequence of Shioda’s MW92 lattice in 89D–91D**. The last three gains use rigorous upper endpoints for the complete external theorem evaluation, rather than the weaker preceding second-moment approximation. Derived section numbers are not attributed as printed records to the original authors.

## What unlocked the new dimensions

We reconstructed the full 92D Mordell-Weil Gram from 180 polynomial sections supplied by Shioda’s published construction. The certificate verifies every curve equation, all 16,110 pairwise heights by polynomial GCD, the complete integer span, and the published determinant `47^14`. The source minimum 16 is inherited.

Actual short dual normals then give explicit sections. Unitary isometries supply further normals, and finite clusters of nearby projected layers improve those sections by the classical antipode construction. The resulting period, shifts, separation guarantees and density squares are checked exactly. Credit for the parent lattice and antipode principle remains with the cited authors; the new results are these explicit packings and certified density bounds.

[Complete proof](PROOF.md) · [Sources and hashes](sources.json) · [All exact results](results.json) · [Normal and cluster labels](sections.json) · [Parent Gram and point coordinates](parent_certificate.npz)

## Reproduce

```sh
python research/mw92_sections_20261115/verify.py
python research/mw92_sections_20261115/verify.py --export /tmp/new_mw92_packings
```

NumPy and python-flint are required. Run without `-O`. The second command requires a new output directory and emits the complete period Gram and rational shifts in period coordinates. No search, optimizer, IEEE-754 assumption or minimum enumeration is needed by this auditor.

## Dimensions still missing an improvement

**48, 54, 55, 56, 58, 63, 64, 72, 80 and 92.** No new SOTA is claimed in those ten dimensions. A 58D half-coset experiment was capped without a complete positive certificate. [Search budgets and outcomes](discovery_summary.json).
