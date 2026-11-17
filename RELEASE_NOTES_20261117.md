# Historical release notices

The following notices were moved out of the README on 17 November 2026. They preserve original numerical values, proof links and caveats; later certificates may supersede them. Consult the current tables and linked packages for the latest results.


> **17 November stronger unrestricted 51D bound:** density **<1.39812853e-7**, an **8.6090%** improvement over the preceding repository certificate. [Complete proof and independent 512-bit replay](verifications/general51_pc29_geomfix_20261117/). The reduction relative to the attributed external numerical benchmark is **8.7465%**, not a certification of that benchmark.

> **17 November stronger unrestricted 52D bound:** density **<9.7440293e-8**, an **8.5210%** improvement over the preceding repository certificate. [Complete proof and independent 512-bit replay](verifications/general52_pc29_geomfix_20261117/). The comparison with the attributed external numerical benchmark is **8.6032%**, not a certification of that benchmark.

> **17 November stronger 45–48D lattice bounds:** **<1.1728551765e-6**, **<7.721596228e-7**, **<5.130033534e-7**, and **<3.4559190288e-7**, respectively. [Three-feature compact kernels, exact dimensional transfers, and complete independent audit](verifications/lattice/compact_power6_45_48_20261117/). The 48D result improves the preceding repository bound by 1.23137%; table comparisons remain against attributed external results.

> **17 November 94/95D packing update:** verified density floors **7.84434540405174e-21** and **2.06770138493108e-21** from anisotropic torus-offset flags. [Exact proofs and independent source replay](research/mw92_flags9495_verified_20261117/). These are periodic existence constructions, not explicit offsets or lattice claims.

> **Earlier 17 November 44D/45D lattice update:** density **<1.776548487e-6** in 44D and **<1.1735752753e-6** in 45D, with [complete proof and independent replay](verifications/lattice/compact_matrix44_45_20261117/). Comparisons below retain Zubrilina’s external bounds.

> **17 November stronger unrestricted 50D bound:** density **<1.99548924e-7**, with [complete exact proof and independent interval replay](verifications/general50_pc29_20261117/). This is 9.3121% below the attributed external numerical benchmark, not a certification of that benchmark.

> **Earlier 17 November dimension-specific update:** lattice-only 46D **<7.766101459e-7**, with [complete independent proof and replay](verifications/lattice/compact_matrix46_new_20261117/). The earlier transferred result remains preserved below and in its proof package.

> **Earlier 17 November lattice transfer:** 46D **<7.790880151e-7** and 47D **<5.186745468e-7**, with [complete independent audit and replay](verifications/lattice/compact_matrix46_47_20261117/). These retain the certified compact 48D plane geometry and use separate target-dimensional scalar proofs.

> **17 November update:** unrestricted53D density **<7.40755520e-8** ([proof](verifications/general53_20261117/)); lattice-only48D density **<3.499004606e-7** ([proof](verifications/lattice/compact_matrix48_20261117/)). These are different problems and are compared with their respective external sources.

> **17 November 51D/52D unrestricted update:** density **<1.52983244e-7** and **<1.06516626e-7**, with [complete independently audited proofs](verifications/general51_52_20261117/). The comparison remains explicitly numerical, not a certified prior bound.

> **Earlier 17 November unrestricted49D/50D update:** certified density **<2.85041752e-7** in 49D and **<2.04993038e-7** in 50D. [Complete proof and audit](verifications/general49_50_20261117/) · [Numerical benchmark comparison](NUMERICAL_BENCHMARKS49_96.md#certified-results-compared-with-these-numerical-benchmarks). These do not extend the lattice-only table.

[Numerical comparison benchmarks for dimensions 49–96](NUMERICAL_BENCHMARKS49_96.md): Afkhami-Jeddi, Cohn et al. (2020), highest-order spinless-bootstrap spectra. These are numerical estimates, **not certified upper bounds**; they are kept separate from the proof tables.

Reproducible certificates for unrestricted sphere packings and lattice-only upper bounds, organized by dimension on the model of KissingNumbers. Upper-bound certificates and packing constructions are separated below. Neither establishes optimality or a kissing-number result. Comparisons are source-attributed; literature searches are not exhaustive priority guarantees.

> **17 November general upper update:** 45D **<1.2613300683e-6** and 48D **<4.0691468600e-7**, with [complete proof and audit scope](verifications/general_angular_20261117/).

> **17 November verified release:** six periodic packing improvements in 85D–90D and lattice-only refinements in 36D–48D. [Exact before/after values, methods and reproduction](CHANGES_20261117.md).


<!-- lower-return-20261115:start -->
> [!IMPORTANT]
> **59D packing certificate — 15 November 2026.** The six-translate Ne64 antipode has packing fraction **2.637602590720003e-12** and center density **274,522.589741…**. A stronger literature-based comparison than the old Elkies guarantee comes from [Nebe's Γ72, §4.2](https://arxiv.org/pdf/1008.2862) and the [classical κ13 section](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html): these imply a 59D packing of center density **269,064.727051…**. Our gain is **2.028457%**. The comparison is derived here from published inputs; it is not a 59D record quoted from those sources. [Proof, source survey and exact replay](constructions/dim59/#literature-derived-comparison).
<!-- lower-return-20261115:end -->

<!-- larger-decreases-20261115:start -->
> [!IMPORTANT]
> **Six further upper bounds — 15 November 2026.** New signed-profile certificates lower the preceding main by **2.79% in 20D**, **1.63% in 19D**, **1.52% in 21D**, **1.20% in 48D**, **1.11% in 18D** and **0.35% in 17D**. All six complete primary and independent proof chains pass. [Exact bounds, comparison and replay](research/larger_decreases_20261115/). The broader **4D** searches did **not** yield a stronger registered result; an exact cyclic counterexample rejects the promising numerical proposal.
<!-- larger-decreases-20261115:end -->

<!-- d4-energy-20261115:start -->
> [!IMPORTANT]
> **4D upper-bound improvement — 15 November 2026.** Every packing in R4 has density **<0.629520189662**, below [Felpix-Studios’ pinned <0.629521681653](https://github.com/Felpix-Studios/sphere-packing/tree/bb97550ba47ad92bacfee567843888be16895a23/preprints/D4). The new degree-16 energy polynomial uses five previously discarded angular Gram credits. The complete base proof is replayed, followed by **308 new interval leaves and 245 exact Gram equations**. [Proof and reproduction](verifications/dim4/refinement_energy20261115/) · [Attribution](verifications/dim4/refinement_energy20261115/NOTICE.md). This is an extension of their certificate; comparison does not establish exhaustive worldwide priority.
<!-- d4-energy-20261115:end -->

<!-- combined-20261115:start -->
> [!IMPORTANT]
> **Nine further upper-bound refinements — 15 November 2026.** The merged profile tools improve **23D, 26D–31D and 33D–34D**, with further reductions up to **2.62% in 29D** and **2.60% in 30D**. All nine registered packages pass complete primary and independent proof programs. [Exact bounds and before/after table](CHANGES_20261115.md) · [Merged research and replay scope](research/combined_branches_20261115/README.md).
<!-- combined-20261115:end -->

> [!IMPORTANT]
> **Sixteen stronger upper bounds — 14 November 2026.** New signed logarithmic profiles and degree-95 radial enrichment improve **22D, 25D, 32D and every dimension 35D–47D**. The largest further reductions from the preceding main are **4.03% in 39D**, **3.97% in 38D** and **3.54% in 37D**. All sixteen complete primary and independent proof chains pass. [Exact bounds, before/after table and what changed](CHANGES_20261114.md) · [Proof and verification scope](verifications/profile_enrichment_20261114/PROOF.md). The earlier same-day 45D–47D results are superseded and retained as historical witnesses.

> [!IMPORTANT]
> **New advances — 12 November 2026.** This repository now includes both density upper bounds and an explicit packing-density lower bound.
>
> - **36D lower bound:** **≈7.3303e-7**, **18.97% above KP36**. A [contact-preserving off-diagonal deformation of the Coxeter–Todd/Turyn seed](constructions/dim36/) has squared minimum **8** and center density **16777216/3176523 ≈5.2816**. Exact reconstruction, full minimum enumeration and an independent audit/rerun are included. This is a construction, not an upper bound or an optimality claim.
> - **9D upper bound:** **<0.190926892**. [Signed-hinge cyclic profile and average-neighbor count](verifications/dim9/), with the external average-neighbor bound **373** explicitly retained as an input.
> - **10D upper bound:** **<0.143330282**. [Log-hinge cyclic profile and root-cut radial repair](verifications/dim10/).
> - **16D upper bound:** **<0.024900975**. [Widened shell count and coupled log-hinge radial kernel](verifications/dim16/widened_reconstruction/), including the independently replayed shell-count certificate.
>
> Upper comparisons use the cited paper's displayed values, not an independently replayed baseline. The 36D minimum proof uses outward-rounded binary64 enumeration under its documented arithmetic/compiler assumptions; the independent rerun uses the same audited implementation. See each package for proof dependencies and exact reproduction commands.

## Packing comparison and research notes

Archived from the README on 17 November 2026. These full notes preserve the original source qualifications, exact values and research history; some descriptions may be superseded by newer linked certificates.

The 84D/85D rows now use stronger classical comparisons than Chen's 2013
guarantees. The 96D comparator is an explicitly reported published
construction, with its exact density recovered from the full paper.
The 49D–53D, 58D and 65D lattices and the 59D–62D, 66D–71D, 78D–79D and
86D–91D periodic packings are explicit. The 57D, 73D–77D, 81D–85D and 93D–96D
bounds are existence proofs. The three Vardy comparisons retain intervals
for the published decimal prefixes; their gains use the upper endpoints.
[Sources, review scope and exact gains](research/covering_and_deep_cosets_20261115/PRIORITY.md).

The new **86D–91D** MW92 constructions supersede the earlier 86D/87D
candidates, which were below the classical BCS88 sections.
The 89D–91D gains use upper endpoints enclosing the full classical Rogers
consequence. Those endpoints are comparison bounds, not densities asserted
to be achieved by older packings.

The recent 88D parent's full computational proof has now been freshly
[reproduced and bound to our catalogue Gram](research/priority_resolution_20261115/REPLAY88.md).
Research remains capped at **96D**, including parent inputs.
The [remaining-dimensions review](research/covering_and_deep_cosets_20261115/REMAINING_DIMENSIONS.md)
accounts for all 49 dimensions: **40 registered best-known claims** and
**nine dimensions without a new registered record claim**.
In **55D**, the comparison target is the full **116,783.149646…** printed in
[Leterrier (2024), Proposition 3.2.22](https://infoscience.epfl.ch/record/307330/files/EPFL_TH10417.pdf).
Our [best certified construction](research/cross_coset55_20261116/PROOF.md)
has minimum10 and center density **91,333.867143…**, which is **below this target**.
The source's displayed construction has a norm6 period vector, so its printed
claim is a challenge target rather than a proved input. Our former 55D
minimum10/local-optimality claim is withdrawn. **No 55D SOTA claim is made**.
The [new MW56 projection search](research/mw56_projection55_20261116/PROOF.md)
uses that full target. Exact certificates exclude all coset subsets in
24 tested projection families; this is a search exclusion, not a new bound.

The [complete 48D–96D audit](research/range48_96_20261115/TABLE.md)
now covers **every dimension in that interval**, with named comparisons.
It adds [explicit 49D–52D minimum-6 lattices](research/range48_96_20261115/P48_PROOF.md),
with complete compatible-coset enumerations above P48, and a stronger
**81D** torus bound above Shioda's MW80.
These improve the published-input comparisons named in the table above.
The audit also includes
Leterrier's published **54D** guarantee and stronger nearby section
controls. [Proofs, source scope and research decisions](research/range48_96_20261115/).

The [remaining-target continuation](research/cm58_20261116/FRONTIER.md) certifies
local optimality against lattice Gram deformations in 54D, 56D, 80D and 92D.
The proposed 92D deformation is now **rejected by an exact shortest vector**;
it was never included in the density table.

The [Craig/code continuation](research/simplex_normals_20261116/PROOF.md#a-separate-attack-on-the-nine-missing-dimensions)
checks exact ceilings on Chen's published guarantee formula in all nine
missing dimensions. That formula cannot improve eight of them. Its remaining
92D case requires an even binary **[93,43,24] code**, which is not constructed
here. These formula ceilings are not upper bounds on actual lattice densities.
A [subsequent exact audit](research/remaining_gaps_20261116/PROOF.md) excludes
cyclic and shift-by-three code constructions for that target, proves Ne64's
local lattice optimality, and rejects 2,237,409 primitive positive forms in
a specified 64D coefficient box. It supplies **no new packing record**.

The [preceding adaptive-layer continuation](research/adaptive_layers_20261115/)
strengthens **13 dimensions**. The larger gains over our preceding bounds
are **15.47% in 50D, 11.80% in 51D, 11.34% in 52D, 18.53% in 53D,
15.47% in 77D and 9.99% in 85D**. Seven covering refinements add
0.011–0.040% in 57D, 81D–82D and 93D–96D. The main table above continues
to compare with named published/derived sources; these incremental gains
use the pinned preceding repository checkpoint. All new proofs and exact
comparisons are included, with two full P48 replay report sets.

The exact 36D density is `(16777216/3176523) * pi^18/18!`. Center density is **not** the packing fraction: it omits the unit-ball volume. The exact ratio to KP36 is `139968/117649`. [Machine-readable 36D construction](constructions/dim36/result.json) and [comparison sources](constructions/dim36/#sources-and-comparison-scope) distinguish the verified improvement over KP36 from the non-exhaustive literature review. The [59D result](constructions/dim59/result.json) records its exact rational squared center density and the separate, finite comparison scope.

For 59D, `V₅₉ = pi^(59/2)/Gamma(61/2)` converts center density into the occupied
fraction shown here. The classical κ13 normal gives the exact comparison
`delta=2^22/(9*sqrt(3))`, stronger than both the earlier internal reference
`262144` and the historical guarantee reported by Vardy. The published κ13
embedding and its primitivity are checked from the catalogue matrices;
Nebe's published Γ72 structure supplies the remaining existence and minimum
inputs. [What other works report, derivation and verification](constructions/dim59/#literature-derived-comparison).
No stronger published 59D lower bound was located in the dated review.
Our bound also exceeds Vardy's reported Elkies guarantee; the unknown actual
minimum of that older lattice does not prevent comparing the proved bounds.

As a separate internal check, our earlier Gamma72 section has center density
262144, below the new construction by exactly **4.722057%**. This earlier
internal reference is useful for tracking our progress; the table uses the
stronger κ13 comparison. Its
[certificate and provenance](research/novel_routes_20261113/audits/PROPAGATED_BENCHMARKS.md#where-our-reference-numbers-come-from)
remain available alongside the named published sources.

Additional merged constructions include an explicit **65D minimum-12 lattice** and periodic packings in 60–62D and 66–71D. Their named comparisons are documented in [the certified higher-dimensional research](research/novel_routes_20261113/README.md).

The [cross-dimension lower-bound research](research/cross_dimension_lower_20261115/)
adds counting-based **existence proofs in 73–78D and 81–85D** and explicit
**L80 antipode certificates in 76–79D**. Its tables identify each derived
product comparison and distinguish existence bounds from explicit packings.
The current best-known claims and external comparisons are collected in the
[dated literature review](research/covering_and_deep_cosets_20261115/PRIORITY.md).

## Independently verified unrestricted 53D improvement

Density <6.8841195e-8, 7.06624096% below the previous repository certificate and 7.1441135129% below the attributed external numerical benchmark (not a certified prior bound). [Self-contained proof](verifications/general53_pc29_geomfix_20261117/) includes the full independent 132,889-node, 512-bit replay, exact PSD and angular signs, Fourier/tail proof, 629 radial annuli and independent odd-dimensional Gamma/pi normalization. The earlier 53D certificate and all other results are preserved.

## Independently verified unrestricted 54D upper bound

Density <4.8695034e-8, 5.4884765414% below the attributed external numerical spinless-bootstrap estimate (N=484), not a certified prior bound. [Self-contained proof](verifications/general54_pc29_geomfix_20261117/) includes the full independent 95,299-node, 512-bit replay, exact PSD and angular signs, Fourier/tail proof, 905 radial annuli from contact 24 through cutoff 250, and exact rational density normalization. The fast Arb polynomial shift is independently reviewed and checked against exact rational compositions. Existing results and the weekly README layout are preserved.

## Dimension-specific compact kernels in 43–44D

Independent complete audits establish lattice-only bounds <2.5363822051e-6 in 43D and <1.7676809433e-6 in 44D, respectively 0.7255% and 0.4991% below the previous repository bounds. Separately optimized positive-semidefinite harmonic matrices use compact powers 6, 9 and 12. The 43D proof closes with zero annular budget and needs no mod-2 profile count; 44D uses a small certified mod-2 annular residual. Each dimension has its own complete rank-two geometry proof, exact cap checks, Fourier factors, scalar tail and normalization. The 384-bit independent replays cover 212,423 and 175,733 nodes. Arithmetic and Jet interval evaluation are shared with the producer; independence is at the algebra, cap and trace-coverage level, not an independent arithmetic backend. External comparator attributions are retained. [Portable proof package](verifications/lattice/compact_power6_adapted_20261117/README.md).

## Stronger 51D bound from three intermediate-power radial features

Unrestricted density <1.38612634e-7 improves the preceding rank-two repository certificate by 0.85844683% and lies 9.5298250776% below the attributed external numerical benchmark (not a certified prior bound). [Self-contained proof](verifications/general51_rank3mid_20261117/) uses profiles with powers 12,18,24, exact 3-by-3 PSD blocks and a 4-by-4 origin block, six tensor-Bernstein controls, independent 78,425-node interval coverage, 633 radial annuli and exact odd-dimensional normalization. The compact correction has degree 48; the Gaussian/Fourier polynomials retain degree 71. All earlier proofs, including the rank-two 51D certificate and 43/44D lattice release, remain preserved.
