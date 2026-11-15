# Worldwide-status check: 15 November 2026

**84D, 85D and 96D are best-known lower bounds in our literature review
dated 15 November 2026.** No stronger published bound was located.
The [current review](../covering_and_deep_cosets_20261115/PRIORITY.md)
registers 32 construction dimensions and gives each source and exact gain.
[All remaining dimensions and the 55D literature correction](../covering_and_deep_cosets_20261115/REMAINING_DIMENSIONS.md).

* **84D/85D:** classical BCS88 sections give stronger comparisons than
  Chen's guarantees. Our gains over those sections are **156%** and
  **32.197832%**.
* **86D/87D:** classical BCS88 sections exceed our candidates, which
  remain rejected.
* **88D:** the recent source computation has been reproduced and bound
  to our parent Gram. BCS's classical theorem already gives density `2^44`.
* **93–95D:** new periodic existence proofs improve the derived MW92
  products named in the complete table.
* **96D:** our current bound exceeds the explicitly reported BCS1982
  construction `3^48/2^24` by **118.554980%**.

We compare **proved lower bounds**, including theorem guarantees. An older
packing's unknown actual minimum and a previous-record table are not
prerequisites for an improvement claim. The numerical percentages below
remain valid historical comparisons with Chen; the stronger BCS sections
are used for the headline gains. The research ceiling remains 96D.

## Direct published comparison in 84D and 85D

Hao Chen, [*On a generalization of Craig lattices*](https://doi.org/10.5802/jtnb.825),
Journal de Théorie des Nombres de Bordeaux **25** (2013), 59–70,
Theorem 2.4 and §3, gives generalized Craig lattices with binary codes.
The parameters specified for these dimensions are `m=10`, `ell=89`,
`k=1`, using repetition-type codes of distance 84. The stated theorem
therefore guarantees

```
delta_n >= 2 * 5^(n/2) / (89^9 * sqrt(n+1)).
```

| Dimension | Chen's theorem guarantee, center density | Our center density | Ratio | Increase over that guarantee |
| --- | ---: | ---: | ---: | ---: |
| 84 | 140,783,142,502.490… | 1,407,374,883,553.280 | 9.99675713× | 899.675713% |
| 85 | 312,965,088,614.593… | 1,453,530,529,565.988… | 4.64438553× | 364.438553% |

The paper's Table 3.1 displays logarithms `37.006` and `38.1616`.
Evaluating its theorem with its specified parameters instead gives
`37.0346836385…` and `38.1872107767…`. This is a numerical discrepancy,
not ordinary rounding. We use the larger theorem-derived guarantees.
The exact comparison concerns these **proved lower bounds**; Chen's actual
lattice minima are not computed here.

Our bounds are the existence constructions already proved in
[COUNTING_PROOF.md, §4](COUNTING_PROOF.md#4-the-80d-parent-gives-deeper-nonsingular-half-cosets).
The gains are certified by rational-square inequalities in
[verify_external_comparisons.py](verify_external_comparisons.py), with
[machine-readable results](external_comparisons.json). No explicit
horizontal representatives have been supplied by the counting proof.

## What was checked

* **Current reference tables.** [Henry Cohn's table](https://cohn.mit.edu/sphere-packing/)
  covers 1–48, 56, 64 and 72; its [MIT archived record](https://dspace.mit.edu/entities/product/d8a88a76-3fa7-4b75-89ef-ea71a0e05b77)
  confirms that scope. [Nebe–Sloane's density table](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html)
  is dated February 2012 and is sparse above 48. NumberDB's
  [T148 generation source](https://numberdb.org/files/T148/generate.py?rev=998)
  explicitly stops at 48. **Omission from these tables is not evidence of
  a record.** The density and unimodular catalogue pages were downloaded.
* **Shimada and Chen.** The full Chen journal paper and
  [arXiv:1110.1950v5](https://arxiv.org/abs/1110.1950v5) were checked.
  The preprint includes the same 84D/85D guarantee and some weaker code
  examples. [Shimada's Frobenius-incidence construction](https://arxiv.org/abs/1007.4614)
  gives the earlier 84D bound `2^34/sqrt(85)`, below Chen's. Searches by
  the paper's title, the displayed density values, the dimensions and the
  constructions did not locate a subsequent stronger 84D/85D guarantee.
* **Current code tables.** Selected Grassl entries relevant to applying
  Chen's theorem at `m=4` were inspected. In particular `[84,21]` has
  reported distance bounds 27–31, and `[85,21]` has 28–32. An upper bound
  of 32 is **not** an existence assertion for a distance-32 code. These
  checks do not optimize every possible code in Chen's construction.
* **New 88D claims.** [Sakamoto, arXiv:2609.15633v2](https://arxiv.org/abs/2609.15633v2)
  claims two extremal even unimodular 88D lattices. The catalogue's
  [unimodular table](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/unimodular.html)
  labels these minimum claims as claims by the authors. The preprint's
  88D result is relevant to possible sections but is not a same-dimensional
  84D/85D density comparison. Its full computational proof has now been freshly
  [reproduced and bound to our parent Gram](../priority_resolution_20261115/REPLAY88.md).
* **Mordell–Weil families.** The [Schütt–Shioda 2017 author draft](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf),
  §13.4, Theorem 13.35, gives the even-dimensional family and the 80D
  guarantee `7^40/(2^40*41^6)`. These are stronger parent data than L80
  in 80D, but their products alone do not establish a denser 84D/85D
  packing. Higher-minimum sections and other extensions have not been
  classified. These are directions for further research, not stronger
  published numerical bounds found in this review.
* **Citation searches.** The OpenAlex DOI record and its indexed citing
  works were inspected. It returned one citing bibliography; Crossref
  reported zero citations. This sparse indexing cannot certify that no
  relevant later work exists. No messages were sent to authors.

Download hashes, locators, URLs and a record of the search queries are in
[literature_sources.json](literature_sources.json). Search coverage has
limits; the best-known claim uses the dated comparison of located papers,
not a missing table entry.

## Status by dimension

| Dimensions | Evidence available | Literature status |
| --- | --- | --- |
| 73–77 | Exact existence proofs; named derived-product comparisons | **Best-known; literature checked 15 November 2026** |
| 78–79 | Explicit L80 antipode certificates; derived-product and section comparisons | **Best-known; literature checked 15 November 2026** |
| **81** | Stronger MW80 covering proof | **Best-known; literature checked 15 November 2026** |
| 82–83 | Exact existence proofs; named derived-product comparisons | **Best-known; literature checked 15 November 2026** |
| **84–85** | Exact gains of 156% and 32.197832% over classical BCS88 sections | **Best-known; literature checked 15 November 2026** |
| **86–87** | Classical BCS88 sections exceed our results | **Our candidates rejected** |
| 93–95 | Stronger periodic existence proofs exceed derived MW92 products | **Best-known; literature checked 15 November 2026** |
| **96** | Exact gain of 118.554980% over the reported BCS96 construction | **Best-known; literature checked 15 November 2026** |

The [59D source investigation](../../constructions/dim59/#literature-derived-comparison)
likewise locates no stronger published lower bound; our result exceeds both
the classical κ13 section and the reported Elkies guarantee. Computing the
actual minimum of the historical E59 lattice is unnecessary for this claim.
For 36D, the July 2026 [Jumagulov preprint](https://arxiv.org/abs/2607.11319)
still names KP36 as the best-known packing; our certificate improves KP36
by 18.97%.

## Remaining research

The computational parent-minimum check in 88D and the BCS96 source comparison
are complete. The 86D/87D candidates are weaker than classical sections.
Structured Mordell–Weil sections and further compatible-coset constructions
remain possible ways to improve the current bounds. A stronger prior paper,
if identified, would update the dated best-known claims and comparisons.

Replay:

```sh
python research/cross_dimension_lower_20261115/verify_counting.py
python research/cross_dimension_lower_20261115/verify80_counting.py
python research/cross_dimension_lower_20261115/verify_external_comparisons.py
```
