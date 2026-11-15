# Best-known lower bounds: literature checked 15 November 2026

**Our certified bounds in 49D–52D, 59D, 81D, 84D, 85D and 96D are the
best-known lower bounds found in the documented literature review.** No
stronger published lower bound was located in those dimensions. The
sources and numerical comparisons are given below; the review date and
claim are also retained in [literature_review.json](literature_review.json).

The comparison criterion is the **proved lower bound**. If an earlier
paper proves density at least A and our certificate proves B>A, we have
improved that paper's lower bound. Computing the earlier packing's
unknown actual density is unnecessary. A previous-record table is also
unnecessary: we compare the papers, their theorem guarantees and named
classical consequences. We revise a best-known claim if a stronger prior
published bound is identified.

## Numerical comparisons and attribution

All values below are **center densities**, not occupied fractions.
The percentage is the gain over the named comparator in that row.
**Derived** means that we calculate the comparison from the cited work;
it does not mean that its authors printed that same-dimensional number.

| Dimension | Our certified lower bound | Named external comparison | Comparison center density | Increase |
| --- | ---: | --- | ---: | ---: |
| 49 | 14,578.768811918… | **Derived:** P48q × A1, from [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf) | 8,417.056098014… | **73.205081%** |
| 50 | 12,625.584147021… | **Derived:** P48q × A2, same source | 4,859.589603973… | **159.807621%** |
| 51 | 11,794.804713454… | **Derived:** P48q × A3, same source | 2,975.878722267… | **296.346955%** |
| 52 | 12,370.495545910… | **Derived:** two Blichfeldt sections of [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0) | 5,237.552386373… | **136.188483%** |
| 59 | 274,522.589741150… | **Derived:** Γ72 perpendicular to √2κ13, from [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862) and [κ13 catalogue data](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html) | 269,064.727051… | **2.028457%** |
| 81 | 1,066,248,855,605.296… | **Derived:** MW80 × A1, from [Schütt–Shioda (2017 author draft), Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 609,520,950,801.816… | **74.932273%** |
| 84 | 1,407,374,883,553.280 | **Derived:** BCS88 perpendicular to √2D4, from [Bos–Conway–Sloane (1982), Theorems 1–2, Eq. (15)](https://neilsloane.com/doc/Me103.pdf) | 549,755,813,888 | **156%** |
| 85 | 1,453,530,529,565.988… | **Derived:** BCS88 perpendicular to √2A3, same source | 1,099,511,627,776 | **32.197832%** |
| 96 | 10,391,088,344,672,478.939… | **Reported:** Bos–Conway–Sloane (1982), Table 1, p.177 | `3^48/2^24`, from Eq. (15) | **118.554980%** |

The exact density squares, comparison provenance and rational gain
intervals are checked by [verify_literature_comparisons.py](verify_literature_comparisons.py)
and retained in [literature_comparisons.json](literature_comparisons.json).
The 49D–52D and 59D packings are explicit; the 81D, 84D, 85D and 96D
bounds are existence proofs. Both establish sphere-packing lower bounds.

For 59D, Vardy's 1999 paper reports a weaker Elkies guarantee with
`log2(delta)` beginning `17.75…`. Our proof improves that published
guarantee too. Its exact actual lattice minimum is not needed for this
conclusion. [Source transcription and exact comparison](../../constructions/dim59/#verified-external-comparison).

## Literature reviewed

* **Reference tables:** Cohn, Nebe–Sloane, the SPLAG third-edition
  preface and NumberDB. These tables have sparse coverage above 48D;
  the review therefore also checked construction papers directly.
* **Construction sources:** Bos–Conway–Sloane (1982), Vardy (1999),
  Nebe's Ne64 and Γ72 constructions, Schütt–Shioda's Mordell–Weil
  chapter, Chen (2013), Shimada (2012), Leterrier (2022), Nigam (2022)
  and the recent papers listed in the linked source registries.
  The classical BCS88 sections give stronger 84D/85D comparisons than
  Chen's guarantees. Leterrier54 supplies the 52D comparison above.
* **Methods:** BCS supplies the P48q parent and its norm-doubler;
  SPLAG's preface describes Elkies's half-lattice-hole extensions.
  We do not claim to originate half-coset lamination. The new results
  are the stated density bounds and their proofs.
* **Covering input:** Coxeter–Few–Rogers (1959), with the precise
  regular-simplex theorem checked in Fejes Tóth's 2022 survey, §2.2.
  Our finite rational estimate and torus-layer argument are given in
  [COVERING_PROOF.md](COVERING_PROOF.md). The original publisher's
  abstract was checked; the full 1959 article was not downloaded.

Dimension-specific queries included written-out and numerical dimensions,
P49–P52, P48/P48q holes and covering radii, MW54/MW80 and their sections,
laminated extensions, citing papers and current tables. Source URLs,
locators, download hashes, queries and access details are recorded in
[sources.json](sources.json), the [48D–96D registry](../range48_96_20261115/sources.json),
the [classical-section registry](../priority_resolution_20261115/sources.json)
and the [earlier literature review](../cross_dimension_lower_20261115/literature_sources.json).

## Earlier values and rejected candidates

The five values announced at checkpoint `e9d03c1` in 49D–52D and 81D
are superseded by the current bounds above. Their exact improvements
are checked against [previous_bounds.json](previous_bounds.json).

**86D and 87D are not improvements:** classical BCS88 sections exceed
our candidates. **63D is not an improvement:** the reported Elkies
guarantee exceeds our Ne64 antipode. These negative decisions remain
in the [complete 48D–96D table](../range48_96_20261115/TABLE.md).
