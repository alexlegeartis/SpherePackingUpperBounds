# Worldwide-status check: 15 November 2026

The [continued primary-source investigation and stronger covering proof](../covering_and_deep_cosets_20261115/)
now gives specific verdicts instead of one blanket priority label:

* **84D/85D:** the stronger comparison is a section of the classical
  BCS88 construction, derived from its published theorem. Our gains are
  **156%** and **32.197832%**. These are record candidates against the
  strongest comparisons located; the earlier Chen comparison is weaker.
* **86D/87D:** the same classical construction gives larger densities
  than our results. **Our candidates are rejected.**
* **88D:** the full recent source computation has now been reproduced
  and bound to our actual parent Gram. Its density `2^44` was already
  obtained by BCS's classical theorem.
* **93–95D:** new periodic existence proofs exceed the derived products;
  no reported dimension-specific record has been identified.
* **96D:** the full BCS1982 paper explicitly reports the construction
  of exact center density `3^48/2^24`. Our new proof exceeds it by
  **118.554980%**, a record candidate against the strongest published
  numerical benchmark located.

See the [current table through 96D](TO96.md). The research ceiling is 96D.
The historical Chen arithmetic below remains valid; it no longer defines
our leading comparison. A record-candidate label records the verified
gain and the stated literature coverage, rather than an exhaustive
classification of all packings.

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
  classified. They remain a material gap in a worldwide claim.
* **Citation searches.** The OpenAlex DOI record and its indexed citing
  works were inspected. It returned one citing bibliography; Crossref
  reported zero citations. This sparse indexing cannot certify that no
  relevant later work exists. No messages were sent to authors.

Download hashes, locators, URLs and a record of the search queries are in
[literature_sources.json](literature_sources.json). Search coverage has
limits; neither a search result count nor a missing table entry is a proof
of global priority.

## Status by dimension

| Dimensions | Evidence available | Worldwide status |
| --- | --- | --- |
| 73–77 | Exact existence proofs; named derived-product comparisons | Previous worldwide record not identified |
| 78–79 | Explicit L80 antipode certificates; derived-product and section-method comparisons | Previous worldwide record not identified |
| 81–83 | Exact existence proofs; named derived-product comparisons | Previous worldwide record not identified |
| **84–85** | Exact gains of 156% and 32.197832% over stronger classical BCS88 sections | **Record candidates against the strongest comparisons located** |
| **86–87** | Classical BCS88 sections exceed our results | **Our record candidates rejected** |
| **93–95** | New periodic existence bounds exceed derived MW92 products | Reported dimension-specific record not identified |
| **96** | Exact gain of 118.554980% over the explicitly reported BCS96 construction | **Record candidate against the strongest published numerical benchmark located** |

The earlier 59D construction retains its separate
[source investigation](../../constructions/dim59/#literature-derived-comparison):
its classical comparison is verified, but the actual density of Elkies's
historical E59 and a current worldwide record are not identified.
For 36D, the July 2026 [Jumagulov preprint](https://arxiv.org/abs/2607.11319)
still names KP36 as the best known packing. That supports the relevance of
our 18.97% KP36 comparison, but is not an independent validation of our
packing certificate or an exhaustive November 2026 priority check.

## Remaining scope

The computational parent-minimum gap in 88D and the missing numerical
BCS96 source comparison have been closed. The 86D/87D record question
has a negative answer for our candidates. The positive 84D/85D/96D
claims are exact improvements over the specified strongest comparisons
located, with published versus derived values identified explicitly.
Further structured Mordell–Weil sections and other constructions remain
relevant checks before assigning an unqualified worldwide SOTA label.
No complete classification is asserted by this literature search.

Replay:

```sh
python research/cross_dimension_lower_20261115/verify_counting.py
python research/cross_dimension_lower_20261115/verify80_counting.py
python research/cross_dimension_lower_20261115/verify_external_comparisons.py
```
