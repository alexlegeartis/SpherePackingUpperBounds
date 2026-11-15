# Published references and repository comparisons, dimensions 57–71

**The strongest 59D comparison calculated here is now the classical kappa13 section of Gamma72, with center density 269064.727051…; our increase is 2.028457%. This comparison is derived from published inputs, not a 59D result quoted from Nebe's paper. The current literature review registers our 59D bound as best-known, with no stronger published bound located.**

[Dated review and current comparisons](../../covering_and_deep_cosets_20261115/PRIORITY.md). Comparison is between proved lower bounds; a previous-record table and an older lattice's unknown actual minimum are not prerequisites.

The same distinction applies to the earlier 65–71D percentage gains: their comparators were our derived constructions or an inherited-minimum section bound. They were not identified published same-dimensional records.

All densities are **center densities** `delta`, not occupied-volume fractions. The latter are `Delta_n = pi^(n/2)/Gamma(n/2+1) * delta`. A comparison of exact rational squared densities proves the percentages below.

## Which previous work is being used?

| Reference | Actual role in this report |
|---|---|
| [V99: Alexander Vardy, *Density doubling, double-circulants, and new sphere packings*](https://doi.org/10.1090/S0002-9947-99-02169-8) | Historical same-dimensional numerical comparison: Elkies's sections E_n, reported by Vardy. These are not being attributed to Vardy as their constructor. |
| [N98: Gabriele Nebe, *Some cyclo-quaternionic lattices*](https://doi.org/10.1006/jabr.1997.7163) | Construction of the 64D parent Ne64. Used with N99 for its minimum theorem, not as a published 59D or 65D comparison result. |
| [N99: Gabriele Nebe, *Construction and investigation of lattices with matrix groups*](https://doi.org/10.1090/conm/249/03759) | Published minimum-12 input for Ne64. The resulting center density in dimension 64 is 3^16. No same-dimensional 59D or 65D result is attributed to this paper. |
| [N12: Gabriele Nebe, *An even unimodular 72-dimensional lattice of minimum 8*](https://doi.org/10.1515/crelle.2011.175) | Published 72D parent and structural input for the classical kappa13 comparison. The derived 59D density 2^22/(9*sqrt(3)) is calculated here, not a 59D numerical result reported by Nebe. Other selected sections and antipode clusters are also calculations in this repository. |
| [K13: Gabriele Nebe and N. J. A. Sloane; classical construction discussed by J. H. Conway and N. J. A. Sloane, *Catalogue of Lattices: KAPPA13 and the Leech lattice; Sphere Packings, Lattices and Groups, Chapter 6*](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html) | Classical normal for a stronger derived 59D comparison. Exact catalogue coordinates certify a primitive kappa13 embedding in Leech. Together with N12 this yields 59D center density 2^22/(9*sqrt(3)) = 269064.727051..., not a quoted 59D record. |
| [S19: Matthias Schütt and Tetsuji Shioda, *Mordell–Weil Lattices, chapter Ranks and Sphere Packings*](https://doi.org/10.1007/978-981-32-9301-4_13) | Published 56D density input for our derived torus-layer and product bounds. The 57D value 151706.083320 is calculated here and must not be attributed to this chapter. |
| [CS96: J. H. Conway and N. J. A. Sloane, *The Antipode Construction for Sphere Packings*](https://neilsloane.com/doc/anti.pdf) | The construction method. This paper's examples do not supply the 59D or 65–71D numerical comparison values used here. |

Full publication details and precise locators are listed under **Source details** below. Parent-lattice papers and construction methods are credited as inputs; using a new section of a published parent does not mean we have improved that paper's density in its original dimension.

## Published numerical comparison: Elkies's sections, reported by Vardy

[Vardy (1999), Section 4, p. 278](https://doi.org/10.1090/S0002-9947-99-02169-8) reports lower bounds for Elkies's sections `E_n`, in a table headed `log2(delta(E_n))`. Here both numerical columns use **log2(center density)**, so their units agree. The ellipses are retained: the printed entries are not exact density formulas.

| Dimension | Published `log2(delta(E_n))`, Vardy §4 | Strongest bound calculated here, `log2(delta)` | What the comparison establishes |
|---|---:|---:|---|
| 57 | 15.37… | 17.210919 | Improves the printed published lower-bound guarantee. |
| 58 | 16.46… | 17.207519 | Improves the printed published lower-bound guarantee. |
| 59 | 17.75… | 18.066565 | 23.67–24.54% improvement over the published guarantee (verified interval). |
| 60 | 19.04… | 19.561223 | Improves the printed published lower-bound guarantee. |
| 61 | 20.54… | 20.962387 | Improves the printed published lower-bound guarantee. |
| 62 | 22.13… | 22.307017 | Improves the printed published lower-bound guarantee. |
| 63 | 23.92… | 23.643718 | Below even the printed published lower bound; no record gain. |

There is **no percentage claim against the actual Elkies lattices** in this table. The original PDF was checked on 15 November 2026. Reference 10, p. 283 attributes the section results to **Elkies, personal communication, February 1996**. No exact E59 section parameters are supplied there.

For 59D, let `B` be the guarantee reported by the printed prefix. Then `2^17.75 <= B < 2^17.76`, and exact integer-power comparisons prove **23.675981% < 100*(delta_ours/B-1) < 24.536216%**. The upper bound here is on the **reported guarantee**, not on the actual Elkies density. This establishes an improvement over the published lower bound; the unknown actual minimum is not required for that conclusion. [Source verification and derivation](../../../constructions/dim59/README.md#verified-external-comparison) · [Exact arithmetic report](../../../constructions/dim59/external_comparison.json).

For **65–71D**, this audit has not identified a published same-dimensional numerical record to use as the comparator. Nebe's 64D and 72D parent results do not fill that gap. Their exact parent densities are respectively `3^16` and `2^36`; no improvement in either parent dimension is claimed here.

## Exact comparisons with calculated reference constructions

The 59D comparator is a consequence of **published Gamma72 structure and classical kappa13**. Its numerical 59D density is derived here. Other rows use the reference constructions calculated in this repository. The links identify each denominator; these percentages must not be described as gains over identified published records.

| Dimension | Reference construction / source | Reference center density | Our certified construction | Center density | Increase |
|---|---|---:|---|---:|---:|
| 59 | [Gamma72 perpendicular to sqrt(2)*kappa13 (derived from N12 + K13)](../../../constructions/dim59/literature_section_comparison.json) | 269,064.727051 | [Six-translate Ne64 antipode](../../../constructions/dim59/result.json) | 274,522.589741 | 2.028457% |
| 60 | [Gamma72 laminated normal rank 12](gamma72_leech_section_rank12.npz) | 524,288.000000 | [Ne64 antipode, 5 translates](../dual64/ne64_simplex_rank4.npz) | 773,597.146771 | 47.551946% |
| 61 | [Gamma72 laminated normal rank 11](gamma72_leech_section_rank11.npz) | 1,048,576.000000 | [Ne64 antipode, 4 translates](../dual64/ne64_simplex_rank3.npz) | 2,043,183.683386 | 94.853180% |
| 62 | [Gamma72 laminated normal rank 10](gamma72_leech_section_rank10.npz) | 2,421,582.543463 | [Ne64 antipode, 3 translates](../dual64/ne64_simplex_rank2.npz) | 5,188,971.717320 | 114.280192% |
| 65 | [Ne64 discriminant-count extension (coset minimum≥8)](../code_routes/NE64_THETA_COUNT.md) | 37,279,553.935621 | [Explicit Ne65, certified depth 26/3](../higher_parent/ne65_verified.json) | 40,837,705.248540 | 9.544512% |
| 66 | [Gamma72 laminated normal rank 6](gamma72_leech_section_rank6.npz) | 77,490,641.390820 | [Gamma72 antipode, 9 translates](../higher_parent/gamma72_root_cluster_rank6.npz) | 82,897,968.185303 | 6.978039% |
| 67 | [Gamma72 laminated normal rank 5](gamma72_leech_section_rank5.npz) | 189,812,531.248503 | [Gamma72 antipode, 6 translates](../higher_parent/gamma72_simplex_rank5.npz) | 247,629,167.976493 | 30.459863% |
| 68 | [Gamma72 laminated normal rank 4](gamma72_leech_section_rank4.npz) | 536,870,912.000000 | [Gamma72 antipode, 5 translates](../higher_parent/gamma72_simplex_rank4.npz) | 812,141,667.469561 | 51.273174% |
| 69 | [Gamma72 laminated normal rank 3](gamma72_leech_section_rank3.npz) | 1,518,500,249.988025 | [Gamma72 antipode, 4 translates](../higher_parent/gamma72_simplex_rank3.npz) | 2,342,574,740.001863 | 54.268973% |
| 70 | [Gamma72 laminated normal rank 2](gamma72_leech_section_rank2.npz) | 4,959,401,049.012478 | [Gamma72 antipode, 3 translates](../higher_parent/gamma72_simplex_rank2.npz) | 7,120,856,007.299829 | 43.582984% |
| 71 | [Gamma72 laminated normal rank 1](gamma72_leech_section_rank1.npz) | 17,179,869,184.000000 | [Gamma72 antipode, 2 translates](../higher_parent/gamma72_simplex_rank1.npz) | 19,663,555,319.937195 | 14.456956% |

In particular, the 59D row is now the six-translate construction `274522.589741…`, replacing the obsolete bare-section value `244079.864699…` that the previous version called the current candidate. Its exact period, translate data and attained minimum are checked by [two implementations](../../../constructions/dim59/README.md#replay). They also verify the earlier internal reference `delta=2^18`. The stronger classical reference now used in this table is checked separately by [verify_literature_section.py](../../../constructions/dim59/verify_literature_section.py).

The Gamma72 section references in 66–71D also attain a ceiling for the **bare-section density formula using inherited minimum 8**, via the known small-dimensional Hermite constants. That ceiling does not bound arbitrary packings or sections with a higher actual minimum, and it is not a literature record.

## Controls that do not give record improvements

| Dimension | Construction tested here | Center density | Why it is not an improvement |
|---|---|---:|---|
| 57 | [Ne64 bare E7 section](ne64_section_E7.npz) | 57,530.175827 | Below our 56D-derived torus bound 151706.083320. |
| 58 | [Ne64 bare E6 section](ne64_section_E6.npz) | 115,060.351653 | Below our Gamma72 section 151348.908966. |
| 63 | [Ne64 antipode, 2 translates](../dual64/ne64_rank1.npz) | 13,105,940.958282 | Below the Elkies lower bound printed in Vardy §4 (`log2(delta)=23.92…`). |
| 64 | [Gamma72 antipode, 9 translates](../higher_parent/gamma72_root_cluster_rank8.npz) | 19,144,465.698823 | Below Nebe's published Ne64 density 43046721; this is a failed control, not our best 64D bound. |

## Where our reference numbers come from

* **57D torus bound:** we derive it from the published Schütt 56D input `5^28/(2^28*29^4)` using separated offsets and colored layers. The 57D number is our calculation, not a result attributed to Schütt. [Derivation](../code_routes/TORUS_LAYER_BENCHMARKS.md) and [exact radius arithmetic](propagated_benchmarks.py).
* **Gamma72 sections:** we embed primitive laminated normals into the scaled Leech block of Nebe's parent and verify them in [verify_gamma72_leech_sections.py](verify_gamma72_leech_sections.py). For rank `k`, `delta^2=2^(72-k)/det(K)`. In 59D, `det(K)=2^23`, hence `delta=2^18=262144`. This specific 59D value was calculated here; it is not cited as a numerical result of Nebe's paper.
* **Stronger 59D classical section:** Nebe's Section 4.2 supplies a primitive `sqrt(2)*Leech` section. The catalogue's classical `kappa13` embeds primitively in Leech with Gram determinant 972, as verified from both published coordinate matrices. The scaled normal and its 59D perpendicular have determinant `2^13*972`. The perpendicular has minimum exactly 8 because it contains Nebe's minimum-8 48D section. Hence `delta=2^22/(9*sqrt(3))=269064.727051…`, and our gain is **2.028457%**. This is a derivation from N12 and K13, not a density quoted from either source. [Proof and source survey](../../../constructions/dim59/README.md#literature-derived-comparison).
* **65D depth-8 extension:** the existence proof using the Ne64 theta series was developed here: [NE64_THETA_COUNT.md](../code_routes/NE64_THETA_COUNT.md). Its center density is `3^(33/2)/2`. The explicit deeper coset gives `3^17/sqrt(10)`; the internal ratio is `sqrt(6/5)`.
* **Additional propagated bounds:** products, torus colorings and averaged ball-window projections are retained in [propagated_benchmarks.json](propagated_benchmarks.json). Their numerical consequences are marked as our computations, even when their parent is published.

## Source details

**V99. Alexander Vardy.** [Density doubling, double-circulants, and new sphere packings](https://doi.org/10.1090/S0002-9947-99-02169-8). Transactions of the American Mathematical Society 351 (1999), 271–283. **Locator:** Section 4, p. 278, the table for n=57,...,63; reference 10 on p. 283 identifies the Elkies section results as personal communication, February 1996.

The table gives lower bounds with ellipses, not exact section determinants or a current record survey. The verified 59D interval 23.675981%–24.536216% compares with the reported guarantee only; it does not certify an improvement over the actual Elkies lattice.

**N98. Gabriele Nebe.** [Some cyclo-quaternionic lattices](https://doi.org/10.1006/jabr.1997.7163). Journal of Algebra 199 (1998), 472–498. **Locator:** The Ne64 lattice; see the catalogue's REFERENCES and GRAM entries.

[Archived source](../dual64/CQ64b.html).

**N99. Gabriele Nebe.** [Construction and investigation of lattices with matrix groups](https://doi.org/10.1090/conm/249/03759). Contemporary Mathematics 249 (1999), 205–219. **Locator:** The minimum-norm reference identified in the Ne64 catalogue's REFERENCES entry.

[Archived source](../dual64/CQ64b.html).

**N12. Gabriele Nebe.** [An even unimodular 72-dimensional lattice of minimum 8](https://doi.org/10.1515/crelle.2011.175). Journal für die reine und angewandte Mathematik 673 (2012), 237–247; preprint 2010. **Locator:** Theorem 2.3 and Section 4 prove minimum 8. Section 4.2, pp. 8–9 of arXiv:1008.2862, gives the primitive sqrt(2)*Leech section K1 and its perpendicular K2 of minimum 8.

[Archived source](../higher_parent/neb72.html).

**K13. Gabriele Nebe and N. J. A. Sloane; classical construction discussed by J. H. Conway and N. J. A. Sloane.** [Catalogue of Lattices: KAPPA13 and the Leech lattice; Sphere Packings, Lattices and Groups, Chapter 6](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html). Catalogue entries last modified 18 July 2014; Conway–Sloane book, third edition, Springer (1999). **Locator:** KAPPA13: dimension 13, minimum 4, Gram determinant 972, BASIS and GRAM; REFERENCES cites SPLAG Chapter 6.

[Archived source](../../../constructions/dim59/comparison_sources/KAPPA13.html).

**S19. Matthias Schütt and Tetsuji Shioda.** [Mordell–Weil Lattices, chapter Ranks and Sphere Packings](https://doi.org/10.1007/978-981-32-9301-4_13). Springer (2019), chapter 13; 56D construction cited at p. 400 by the Nebe–Sloane catalogue. **Locator:** 56D input, delta=5^28/(2^28*29^4); catalogue row 56 cites M. Schütt, MW19 p. 400.

[Archived source](../code_routes/density_catalogue.html).

The 2019 final edition was not available in the inspected publisher preview. On 15 November 2026 the archived author draft dated 16 March 2017 was checked in full: Theorem 13.35 and equation (13.24), p. 357, independently confirm the 56D input, but do not give exact 59D section parameters.

**CS96. J. H. Conway and N. J. A. Sloane.** [The Antipode Construction for Sphere Packings](https://neilsloane.com/doc/anti.pdf). Inventiones Mathematicae 123 (1996), 309–313. **Locator:** Section 2, Theorem 1 and density formula (1).

**TABLE. Gabriele Nebe and N. J. A. Sloane.** [Table of Densest Packings Presently Known](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html). Catalogue of Lattices; archived copy retained by the earlier audit. **Locator:** Rows 56, 64, 72 and the small-dimensional inputs used in direct products.

[Archived source](../code_routes/density_catalogue.html).

## Replay and claim status

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/propagated_benchmarks.py
python constructions/dim59/verify.py
python constructions/dim59/verify_independent.py
python constructions/dim59/verify_external_comparison.py
python constructions/dim59/verify_literature_section.py
```

The first command regenerates this page, the exact comparison arithmetic and source attribution from [benchmark_sources.json](benchmark_sources.json). It consumes previously verified construction data; it does not replace their minimum proofs. The next two commands replay the 59D construction and its internal reference independently. The fourth checks the historical reported-guarantee interval using exact rational powers; its source transcription is an input. The last verifies the classical kappa13 embedding and the stronger derived 59D comparison.

**Status:** certified constructions with named published and derived comparisons. The dated literature review registers 59D as best-known. Other rows here establish the specified numerical comparisons. Sources and proved guarantees determine the comparison; no previous-record registry is required.
