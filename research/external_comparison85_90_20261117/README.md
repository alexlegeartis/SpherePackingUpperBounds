# External comparisons for the 85–90D packing improvements

The principal packing table compares against external work preceding our constructions, not earlier results by this repository. This corrects the17 November release's choice of comparison column. Our new packings and their original proof packages are unchanged.

## Sources and exact interpretation

* **85–87D:** derived sections of the Bos–Conway–Sloane construction, perpendicular to primitive dual normals sqrt(2)A3, sqrt(2)A2 and sqrt(2)A1. [Bos, Conway and Sloane (1982), Theorems1–2 and Eq.(15)](https://neilsloane.com/doc/Me103.pdf); [publisher DOI](https://doi.org/10.1112/S0025579300012262). These dimension-specific numbers are derived here, not asserted to be printed records in the paper. The exact center-density squares are2^80,2^84/3 and2^85. [Construction, code parameters, dual-section derivation and exact matrix audit](../priority_resolution_20261115/CLASSICAL_SECTIONS.md).
* **88D:** minimum8 unimodular lattices have center-density square2^88. [Sakamoto (2026), Theorem1.1 and the extremal-minimum formula](https://arxiv.org/html/2609.15633v2) explicitly give such rank88 lattices. This is the same packing density as the older BCS construction, which need not be integral. No claim of global optimality is intended.
* **89–90D:** derived MW92 sections using Levenshtein spherical linear programming and [Cohn–Zhao (2014), Proposition2.1](https://arxiv.org/html/1212.5966). The parent is still Schütt–Shioda's MW92 lattice. [Full exact derivation and independent certificate](KL_PROOF.md). This new evaluation of pre-existing methods strengthens the earlier Rogers-derived comparison. Neither these dimension-specific numerical values nor global optimization are attributed to the original papers.

The 17 November source check reconfirmed Sakamoto's theorem, and a separate direct BCS PDF check confirmed Theorem2 and Eq.(15), printed pp.174–175 (SHA256 `8ce46a5a6af6e58bd2a2f7ab2ceefcd629a106b9c61b0ddae99c59940054d448`). The archived Schütt–Shioda PDF was not retrievable through the fresh web check; its prior precise locators, PDF hash and exact derivations remain in the linked source registry. The exact Rogers guarantee enclosures were freshly replayed, then superseded in89D/90D by the stronger spherical-LP consequence. The table uses the strongest external-derived guarantees certified in this follow-up. Published optimized Euclidean LP calculations could give stronger section guarantees; this is not an exhaustive global record survey or an optimization claim.

## Reproducible comparison arithmetic

Run `python research/external_comparison85_90_20261117/verify.py` from the repository root. It reads the unchanged frozen candidates and exact external center-density squares, and proves each displayed two-decimal percentage by rational squaring. The common unit-ball volume cancels in each same-dimensional density ratio. For89D/90D it uses the exact integer center-density-square floors certified by check_kl.py and the separately implemented independent_kl.py. No floating-point result is trusted for the printed strict gain.

| Dimension | Increase over named external comparison | Separate increase over our preceding repository construction |
| ---: | ---: | ---: |
|85|>118.72%|50.43%|
|86|>183.82%|71.27%|
|87|>292.21%|47.73%|
|88|>388.59%|38.68%|
|89|>265.02%|16.86%|
|90|>124.41%|4.69%|

The original release's RESULTS.json and replay outputs preserve the internal-step comparisons as historical data. They are not the source of the corrected external-comparison column. This separation preserves all original verification-package hashes.
