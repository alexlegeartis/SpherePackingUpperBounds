# Best-known lower bounds: literature checked 16 November 2026

**Our certified bounds in 49D–53D, 57D–62D, 65D–71D, 73D–79D, 81D–91D and 93D–96D are the
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

The 16 November update carries forward the previous source review and adds
the new 58D construction and its dimension-specific literature searches.

## Numerical comparisons and attribution

All values below are **center densities**, not occupied fractions.
The percentage is the gain over the named comparator in that row.
**Derived** means that we calculate the comparison from the cited work;
it does not mean that its authors printed that same-dimensional number.

| Dimension | Our certified lower bound | Named external comparison | Comparison center density | Increase |
| --- | ---: | --- | ---: | ---: |
| 49 | 14,578.768811918… | **Derived:** P48q × A1, from [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf) | 8,417.056098014… | **73.205081%** |
| 50 | 14578.768811917893… | **Derived:** P48q × A2, same source | 4,859.589603973… | **200.000000%** |
| 51 | 13186.992560309189… | **Derived:** P48q × A3, same source | 2,975.878722267… | **343.129368%** |
| 52 | 13773.364524023099… | **Derived:** two Blichfeldt sections of [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0) | 5,237.552386373… | **162.973303%** |
| 53 | 20890.804898853429… | **Derived:** one Blichfeldt section of [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0) | 17532.82564220611 | **19.152528%** |
| 57 | 163093.23475965115… | **Derived:** MW56 × A1 — [Schütt–Shioda, Theorem 13.35, p=29](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 98106.6069059854 | **66.240827%** |
| 58 | 161833.931161289967 | **Derived:** Γ72 perpendicular to a rank-14 Leech block — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section audit](../novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 151348.908966445256 | **6.927716%** |
| 59 | 274,522.589741150… | **Derived:** Γ72 perpendicular to √2κ13, from [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862) and [κ13 catalogue data](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html) | 269,064.727051… | **2.028457%** |
| 60 | 773597.1467710571 | **Reported:** Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 19.04… | `log₂ δ=19.04…` | **>42.525785%** |
| 61 | 2043183.683386255 | **Reported:** Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 20.54… | `log₂ δ=20.54…` | **>33.088644%** |
| 62 | 5188971.717319936 | **Reported:** Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 22.13… | `log₂ δ=22.13…` | **>12.273483%** |
| 65 | 40837705.24854031 | **Derived:** Γ72/Leech-block section, normal rank 7 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 33554432.0 | **21.705846%** |
| 66 | 82897968.18530319 | **Derived:** Γ72/Leech-block section, normal rank 6 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 77490641.39082 | **6.978039%** |
| 67 | 247629167.9764929 | **Derived:** Γ72/Leech-block section, normal rank 5 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 189812531.248503 | **30.459863%** |
| 68 | 812141667.4695614 | **Derived:** Γ72/Leech-block section, normal rank 4 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 536870912.0 | **51.273174%** |
| 69 | 2342574740.001862 | **Derived:** Γ72/Leech-block section, normal rank 3 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 1518500249.98802 | **54.268973%** |
| 70 | 7120856007.29983 | **Derived:** Γ72/Leech-block section, normal rank 2 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 4959401049.01248 | **43.582984%** |
| 71 | 19663555319.93719 | **Derived:** Γ72/Leech-block section, normal rank 1 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](../../research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | 17179869184.0 | **14.456956%** |
| 73 | 56109217798.08641 | **Derived:** Γ72 × A1 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | 34359738368.0 | **63.299316%** |
| 74 | 52900277856.1331 | **Derived:** Γ72 × A2 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | 19837604196.0499 | **166.666667%** |
| 75 | 52900277856.1331 | **Derived:** Γ72 × A3 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | 12148001999.9042 | **335.464843%** |
| 76 | 61083979320.88889 | **Derived:** Γ72 × D4 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | 8589934592.0 | **611.111111%** |
| 77 | 70533703808.177467… | **Derived:** Bare L80 section, normal rank 3 — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf); [normal matrix audit](../../research/cross_dimension_lower_20261115/l80_results.json) | 22443687119.2346 | **214.269681%** |
| 78 | 104731824371.1325 | **Derived:** Bare L80 section, normal rank 2 — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf); [normal matrix audit](../../research/cross_dimension_lower_20261115/l80_results.json) | 79350416784.1997 | **31.986483%** |
| 79 | 299178543929.5057 | **Derived:** L80 perpendicular to a minimum vector — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) | 274877906944.0 | **8.840520%** |
| 81 | 1066394783124.8914… | **Derived:** MW80 × A1, from [Schütt–Shioda (2017 author draft), Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 609,520,950,801.816… | **74.956215%** |
| 82 | 1040341684972.145… | **Derived:** MW80 × A2 — [Schütt–Shioda, Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 351907085022.145 | **195.629650%** |
| 83 | 1112627538435.653 | **Derived:** MW80 × A3 — [Schütt–Shioda, Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 215498198793.618 | **416.304797%** |
| 84 | 1,407,374,883,553.280 | **Derived:** BCS88 perpendicular to √2D4, from [Bos–Conway–Sloane (1982), Theorems 1–2, Eq. (15)](https://neilsloane.com/doc/Me103.pdf) | 549,755,813,888 | **156%** |
| 85 | 1598670440672.3403… | **Derived:** BCS88 perpendicular to √2A3, same source | 1,099,511,627,776 | **45.398230%** |
| 86 | 3.532468471752e+12 | **Derived:** BCS88 perpendicular to √2A2 — [Bos–Conway–Sloane (1982), Theorems 1–2](https://neilsloane.com/doc/Me103.pdf); [derived section](../../research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | 2.539213337094e+12 | **39.116648%** |
| 87 | 1.220042544453e+13 | **Derived:** BCS88 perpendicular to √2A1 — [Bos–Conway–Sloane (1982), Theorems 1–2](https://neilsloane.com/doc/Me103.pdf); [derived section](../../research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | 6.219777023951e+12 | **96.155351%** |
| 88 | 5.045199642804e+13 | **Reported:** Minimum-8 unimodular lattices — [Sakamoto (2026)](https://arxiv.org/abs/2609.15633v2); same density as [BCS1982](https://neilsloane.com/doc/Me103.pdf) | 1.759218604442e+13 | **186.786397%** |
| 89 | 2.039717945339e+14 | **Derived:** MW92 / full Rogers-simplex section — [Schütt–Shioda, Theorem 13.35](<https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf>) + [Rogers (1958), via Fejes Tóth §2.2](https://arxiv.org/abs/2202.11358); [derived interval](../../research/mw92_sections_20261115/PROOF.md#external-comparisons-including-the-full-rogers-consequence) | ≈[7.144650443536e+13, 7.146541762714e+13] | **>185.413282%** |
| 90 | 7.844458559109e+14 | **Derived:** MW92 / full Rogers-simplex section — [Schütt–Shioda, Theorem 13.35](<https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf>) + [Rogers (1958), via Fejes Tóth §2.2](https://arxiv.org/abs/2202.11358); [derived interval](../../research/mw92_sections_20261115/PROOF.md#external-comparisons-including-the-full-rogers-consequence) | ≈[3.612830642882e+14, 3.613465011337e+14] | **>117.089650%** |
| 91 | 2.478889293267e+15 | **Derived:** MW92 / full Rogers-simplex section — [Schütt–Shioda, Theorem 13.35](<https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf>) + [Rogers (1958), via Fejes Tóth §2.2](https://arxiv.org/abs/2202.11358); [derived interval](../../research/mw92_sections_20261115/PROOF.md#external-comparisons-including-the-full-rogers-consequence) | ≈[1.861659817306e+15, 1.861822446835e+15] | **>33.143162%** |
| 93 | 8743033982260848.9… | **Derived:** MW92 × A1 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 4.88702544074196e+15 | **78.902977%** |
| 94 | 8584803057662693.3… | **Derived:** MW92 × A2 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 2.82152545374892e+15 | **204.261053%** |
| 95 | 8868241996705554.8… | **Derived:** MW92 × A3 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | 1.72782441448991e+15 | **413.260602%** |
| 96 | 10395177038683104.0… | **Reported:** Bos–Conway–Sloane (1982), Table 1, p.177 | `3^48/2^24`, from Eq. (15) | **118.640977%** |

The exact density squares, comparison provenance and rational gain
intervals are checked by [verify_literature_comparisons.py](verify_literature_comparisons.py)
and retained in [literature_comparisons.json](literature_comparisons.json).
The 49D–53D, 58D and 65D lattices and the 59D–62D, 66D–71D, 78D–79D and
86D–91D periodic packings are explicit. The other registered bounds are existence
proofs. Both establish sphere-packing lower bounds. For 60D–62D, the stated
strict gains hold even at the upper ends of the source's decimal-prefix
intervals; the full outward gain intervals are retained in the exact report.

The [remaining-dimensions review](REMAINING_DIMENSIONS.md) accounts for the
whole range. The subsequent [adaptive-layer continuation](../adaptive_layers_20261115/)
strengthens 13 bounds and adds an explicit 53D lattice to the dated claims,
bringing that checkpoint to 33. The [MW92 continuation](../mw92_sections_20261115/)
adds six explicit construction dimensions, **86D–91D**, for a checkpoint total of **39**. The [new CM58 construction](../cm58_20261116/)
adds 58D, bringing the current total to **40**.
The 89D–91D gains enclose and exceed the full Rogers-simplex consequence;
they do not rely on weakening its numerical evaluation. The former 55D
literature correction is superseded by the short-period correction below.

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
* **Leterrier's 2024 thesis:** Proposition 3.2.22, pp.155–156, prints
  `delta²=5^55/(2^56*3^24)` in 55D, but its displayed lattice has a
  norm 6 period vector and does not establish that claim. Our prior
  minimum 10/local-optimality claim is withdrawn. A
  [cross-coset tilted lattice](../cross_coset55_20261116/PROOF.md) gives
  center density 91,333.867143… with **no new 55D record claim**.
  The full printed 116,783.149646… remains the comparison target in the
  complete range table and the [MW56 projection search](../mw56_projection55_20261116/PROOF.md).
  Its table on p.204 also reports
  weaker 68D and 74D bounds. See [review details](remaining_dimensions_review.json).
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

The **earlier 86D and 87D candidates were not improvements** over the
classical BCS88 sections. The new MW92 constructions now exceed those
comparisons. **63D still has no improvement:** the reported Elkies
guarantee exceeds our Ne64 antipode. The [complete 48D–96D table](../range48_96_20261115/TABLE.md)
records the current replacements and the remaining negative decisions.
