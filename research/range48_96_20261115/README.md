# Lower bounds in every dimension 48D–96D

**Latest:** [explicit MW92 antipodes](../mw92_sections_20261115/) add **86D–91D**. The [new CM58 lattice](../cm58_20261116/) adds 58D, bringing the dated total to **40**; the nine remaining dimensions are 48, 54–56, 63–64, 72, 80 and 92.

This continuation checks **all 49 dimensions**, rather than only the
84D/85D/96D candidates. [The complete table](TABLE.md) names the source
or construction of every comparison. All values are **center densities**;
multiply by `pi^(n/2)/Gamma(1+n/2)` for occupied fractions.

The table now includes the [subsequent deeper-coset and covering pass](../covering_and_deep_cosets_20261115/).
Its exact verifiers are separate from this aggregation script. The [adaptive continuation](../adaptive_layers_20261115/) now strengthens13 entries and adds an explicit 53D lattice.

New positive work in the other dimensions:

* **49–52D:** [deeper compatible half-cosets above P48](../adaptive_layers_20261115/P48_PROOF.md)
  now give **explicit minimum-6 lattices** with densities
  **14,578.768812; 14,578.768812; 13,186.992560; 13,773.364524**.
  These supersede the earlier explicit root-metric constructions
  and are checked by full coset enumerations and fitted vertical metrics.
  The same counting calculation in 53D/54D is weaker than known
  constructions, so those variants are not improvements.
* **81D:** the covering-based two-color extension above **Shioda's MW80**
  now gives **1,066,394,783,124.8914…**, superseding the earlier
  996,397,626,643.8245… torus checkpoint. The same covering argument
  improves 57D, 82D and 93–96D. These are new calculations from named
  classical inputs, not published numerical records attributed to them.
* **53D:** an [explicit primitive section](../adaptive_layers_20261115/SECTION53_PROOF.md)
  of Leterrier's published54D Gram now gives **20,890.8048988534…**.
* **54D literature correction:** Leterrier's 2022 Corollary 1.4 supplies
  the proved guarantee `sqrt(3)*5^27/(2^27*3^13)` =
  **60,306.5591597049…**, matching Elkies's older MW54 guarantee.
  Its derived 53D section is stronger than the new 53D counting variant.

The 49–52D shifts are explicit and every required coset has been fully
enumerated. The 81D torus shifts remain existential.
**49D–53D, 57D–62D, 65D–71D, 73D–79D, 81D–91D and 93D–96D are best-known lower bounds in
the literature review dated 16 November 2026:** no stronger published
bound was located. The [current review and comparison table](../covering_and_deep_cosets_20261115/PRIORITY.md)
name each source and distinguish reported numbers from derived bounds.
Comparison is between proved guarantees; a prior lattice's unknown
actual minimum and a previous-record table are not prerequisites.

## Decisions across the full range

| Dimensions | Best route currently checked / decision |
| --- | --- |
| **48** | Classical P48 minimum-6 lattices; no improvement found |
| **49–52** | New explicit P48 compatible-coset lattices, with full coset enumerations; seek deeper compatible shifts and further external comparisons |
| **53–56** | New explicit 53D section; retain Leterrier54, Leterrier2024's55D lamination and MW56 |
| **57–58** | Covering-based MW56 extension in 57D; explicit CM plane above MW56 now improves 58D |
| **59–62** | Retain existing Ne64 antipodes; their comparisons have the scope recorded in the source audit |
| **63** | Our Ne64 antipode is below Elkies's historical guarantee; do not claim an improvement |
| **64** | Nebe's Ne64, center density `3^16`; no improvement found |
| **65** | Existing explicit Ne65 remains the best repository construction checked here |
| **66–71** | Existing Gamma72 antipodes remain strongest in this scan; no new gain found |
| **72** | Nebe's Gamma72, center density `2^36`; no improvement found |
| **73–77** | Retain73D–76D counting bounds; the adaptive exceptional-coset argument improves77D |
| **78–79** | Existing L80 antipodes remain strongest; comparisons now include actual bare sections and MW80 section controls |
| **80** | Use the stronger published **MW80**, not L80, as the density input |
| **81** | New MW80 two-color torus bound exceeds our previous counting bound |
| **82** | New covering-based three-color MW80 extension exceeds the L80 counting checkpoint |
| **83–85** | Retain83D/84D counting bounds; the fitted adaptive85D metric gives a stronger bound |
| **86–88** | Explicit MW92 antipodes now exceed BCS88 and its sections |
| **89–91** | Explicit MW92 antipodes exceed the full Rogers-based section consequences |
| **92** | Classical MW92, center density `2^92/47^7`; no improvement found |
| **93–96** | Stronger covering-based torus proofs; 96D exceeds the explicitly reported BCS96 construction by 118.640977% |

## What the complete table actually compares

The external-input family contains the classical P48, Leterrier54,
Shioda MW56/MW80/MW92, Nebe Ne64/Gamma72, BCS88 and BCS96 inputs.
For each lattice input, recursive sections use Blichfeldt's Hermite
upper bound, checked by exact rational-power inequalities. We also
include explicit Gamma72/Leech-block sections, the classical kappa13
normal, BCS88 root sections, and bare sections using the audited L80
normal embeddings. These **derived** numbers are not attributed as
dimension-specific records printed by the source authors.

Vardy's published prefixes for Elkies's 57–63D guarantees are retained
as lower bounds on those reported lattices. Exact integer-power checks
give conservative rational squared-density floors. Neither those floors
nor the printed prefixes are upper bounds on the actual Elkies densities.
The 63D comparison already suffices to reject our weaker candidate.

The family is closed under products with the named classical lattices
in dimensions 1–24. The same product closure is applied before and
after adding the repository constructions, so inherited product gains
do not masquerade as independent constructions. **This is a reproducible
comparison family, not a complete catalogue of worldwide records.**

## Broad torus scan

The checker tests up to twelve added dimensions above every even-dimensional
parent, stopping at 96D. It uses period `2R` (index `2^r`, same-color
minimum ratio 4), the exact small periods from the preceding 93–96D
proof, and the standard `(r+1)`-coloring of A_r. For r≥3 the latter has
same-color minimum ratio at least 2: a norm-2 root changes color, while
every positive root-lattice norm is an even integer.

Normalize the horizontal lattice to covolume `1/delta`, minimum at
least 4. If only a center-density guarantee is known, first normalize
its actual minimum to 4 and then enlarge the lattice until its covolume
is `1/delta`; its minimum remains at least 4. Thus a density guarantee
alone is enough, without assuming equality in Leterrier's estimate.

For N colors, pick rational d satisfying

```
(N-1)^2 * delta^2 * V_n^2 * d^n < 1.
```

This is the square of the greedy forbidden-volume condition. With
`t=max(4-d,4/lambda)` and vertical root Gram Q, the resulting periodic
packing has separation at least 2 and center-density square

```
delta_new^2 = delta^2 * (2/t)^r / det(Q).
```

The [preceding torus proof](../priority_resolution_20261115/TORUS_PROOF.md)
gives the full argument. Every volume inequality, same-color condition
and density square in this scan is checked exactly. Horizontal offsets
are existential. For 81D, the normalized depth is `625793/250000` and
the vertical squared minimum is `374207/250000`; the two-color period
has minimum four times the vertical minimum.

## Sources and reproduction

* **P48 and BCS:** Bos–Conway–Sloane,
  [*Further lattice packings in high dimensions* (1982)](https://neilsloane.com/doc/Me103.pdf),
  Theorems 1–2, equation (15), Table 1 and §5.
* **54D:** Gauthier Leterrier,
  [*On the Mordell–Weil lattice of y²=x³+bx+t^(3^n+1) in characteristic 3*](https://doi.org/10.1007/s40993-022-00321-0),
  Research in Number Theory **8**, 23 (2022), **Corollary 1.4 and
  the table on p.3**, at n=3. The source states a lower guarantee,
  not equality for the actual lattice density.
* **MW56/MW80/MW92:** Schütt–Shioda,
  [2017 author draft](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf),
  Theorem 13.35 and equation (13.24), p.357, at p=29,41,47.
* **Ne64:** Nebe, [*Some cyclo-quaternionic lattices* (1998)](https://doi.org/10.1006/jabr.1997.7163),
  and [*Construction and investigation of lattices with matrix groups* (1999)](https://doi.org/10.1090/conm/249/03759);
  [catalogue and minimum references](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html).
* **Gamma72:** Nebe,
  [*An even unimodular 72-dimensional lattice of minimum 8*](https://arxiv.org/abs/1008.2862),
  Theorem 2.3 and §4.2.
* **L80:** Bachoc–Nebe,
  [*Extremal lattices of minimum 8 related to the Mathieu group M22* (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf).
* **Historical 57–63D guarantees:** Vardy,
  [*Density doubling, double-circulants, and new sphere packings* (1999)](https://doi.org/10.1090/S0002-9947-99-02169-8),
  §4, p.278. Attribution and limitations are retained in the
  [existing source audit](../novel_routes_20261113/audits/PROPAGATED_BENCHMARKS.md).
* **Small lattice inputs:** the
  [Nebe–Sloane density catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html)
  and Conway–Sloane's classical lattice table. We use the lattice values,
  with no rounded nonlattice values substituted as exact numbers.

The additional Nigam [2022 preprint](https://arxiv.org/abs/2206.06525v2)
was checked in full. Its displayed computations do not supply stronger
48–96D numerical guarantees; no extra density input is inferred from
an uncomputed family member. Cohn's current table covers 1–48 and
56/64/72, leaving most of this interval untabulated; that omission is
not evidence of a record. The wider searches and remaining source
limitations are recorded in the linked earlier literature reviews.
The [Klartag 2026 theorem](https://doi.org/10.1007/s00222-026-01412-w)
gives a general bound `c*n^2*2^-n`; an unevaluated universal constant is
not inserted as a numerical 48–96D benchmark. Download hashes and query
records are in [sources.json](sources.json).

Replay from the repository root, using Python's standard library:

```sh
python research/range48_96_20261115/verify_range.py
```

It regenerates the calculations in memory and rejects stale
[results.json](results.json) or [TABLE.md](TABLE.md).
`--write-result` deliberately updates those two generated files.
Existing certificate reports are consumed as audited inputs, not
re-proved by this arithmetic aggregation. The parent minimum theorems,
theta modularity and finite quadratic-space transitivity remain stated
mathematical dependencies. Every input and output dimension is at most 96.

For the new explicit lattices, also run the
[exact setup audit and full interval replay](P48_PROOF.md).
The full replay uses NumPy, SymPy, python-flint and GNU C++; the
standard-library range script only aggregates the completed certificate.
