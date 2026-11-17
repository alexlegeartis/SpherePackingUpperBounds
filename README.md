# Sphere-packing density bounds and constructions

Reproducible certificates for unrestricted sphere packings and lattice-only upper bounds, organized by dimension on the model of KissingNumbers. Upper-bound certificates and packing constructions are separated below. Neither establishes optimality or a kissing-number result. Comparisons are source-attributed; literature searches are not exhaustive priority guarantees.

> **17November verified release:** six periodic packing improvements in85D–90D and lattice-only refinements in36D–48D. [Exact before/after values, methods and reproduction](CHANGES_20261117.md).

<!-- lattice-upper-20261116:start -->
## Lattice-only upper bounds

Comparison sources are named per row. Dimensions36–48 compare with the pinned preceding lattice-only repository certificates, which remain unchanged. These are lattice-only bounds, not unrestricted-packing claims or exhaustive priority guarantees. [Scope and reproduction](LATTICE_RESULTS.md).

| n | Lower | Prior / source | New | Gain ≈ | Gap ≈ | Method |
|---:|---:|---:|---:|---:|---:|---|
| 14 | 0.0216240 | 0.0352393<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.0347128** | 1.49% | 3.87% | [Near-frame index exclusion + weighted mod-2 class profiles](verifications/lattice/dim14/README.md) |
| 15 | 0.0168575 | 0.0256119<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.0255112** | 0.39% | 1.15% | [Near-frame index exclusion + weighted mod-2 class profiles](verifications/lattice/dim15/README.md) |
| 16 | 0.0147081 | 0.0187462<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.0186552** | 0.49% | 2.25% | [Near-frame index exclusion + weighted mod-2 class profiles](verifications/lattice/dim16/README.md) |
| 17 | 0.00881131 | 0.0138840<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.0138325** | 0.37% | 1.02% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim17/README.md) |
| 18 | 0.00592836 | 0.0103705<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.0103439** | 0.26% | 0.60% | [Near-frame index exclusion + weighted mod-2 class profiles](verifications/lattice/dim18/README.md) |
| 19 | 0.00412080 | 0.00782448<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.00780844** | 0.21% | 0.43% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim19/README.md) |
| 20 | 0.00322586 | 0.00595809<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.00595464** | 0.06% | 0.13% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim20/README.md) |
| 21 | 0.00246588 | 0.00456046<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.00455796** | 0.05% | 0.12% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim21/README.md) |
| 22 | 0.00212766 | 0.00348597<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.00348294** | 0.09% | 0.22% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim22/README.md) |
| 23 | 0.00190532 | 0.00261354<br>[Signed logarithmic profile + degree95 radial enrichment — Kravatskiy (2026), unrestricted bound](verifications/dim23/refinement_combined20261115/) | **≤0.00260841** | 0.20% | 0.73% | [Conditioned nested Gram-moment class profiles](verifications/lattice/dim23/README.md) |
| 25 | 6.77212e-4 | 0.00132931<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤0.00132524** | 0.31% | 0.63% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim25/README.md) |
| 26 | 2.69220e-4 | 9.15120e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤9.07210e-4** | 0.86% | 1.22% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim26/README.md) |
| 27 | 1.28675e-4 | 6.31343e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤6.22939e-4** | 1.33% | 1.67% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim27/README.md) |
| 28 | 6.97587e-5 | 4.37001e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤4.29640e-4** | 1.68% | 2.00% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim28/README.md) |
| 29 | 2.78789e-5 | 3.03843e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤2.98026e-4** | 1.91% | 2.11% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim29/README.md) |
| 30 | 1.44286e-5 | 2.12403e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤2.07933e-4** | 2.10% | 2.26% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim30/README.md) |
| 31 | 1.18377e-5 | 1.49409e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤1.46088e-4** | 2.22% | 2.41% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim31/README.md) |
| 32 | 1.10407e-5 | 1.05765e-4<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤1.01512e-4** | 4.02% | 4.49% | [Multiwindow cross-class pair-sum coupling + mod-2 profiles](verifications/lattice/dim32_threewindow/README.md) |
| 33 | 4.14068e-6 | 7.55168e-5<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤7.38272e-5** | 2.24% | 2.37% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim33/README.md) |
| 34 | 1.76697e-6 | 5.43082e-5<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤5.31485e-5** | 2.14% | 2.21% | [Cross-class pair-sum angular coupling + mod-2 profile bounds](verifications/lattice/dim34/README.md) |
| 35 | 9.46190e-7 | 3.93267e-5<br>[Lattice class-profile LP — Zubrilina (2026)](https://arxiv.org/abs/2610.35830) | **≤3.76648e-5** | 4.23% | 4.33% | [Multiwindow cross-class pair-sum coupling + mod-2 profiles](verifications/lattice/dim35_threewindow/README.md) |
| 36 | 7.33034e-7 | 2.74287e-5<br>[Multiwindow cross-class pair-sum coupling + mod-2 profiles — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim36_threewindow/README.md) | **≤2.74025e-5** | 0.10% | 0.10% | [Nonnegative quadratic Fourier multiplier + rank-two PSD kernel + reoptimized mod-2 angular-cut profiles](verifications/lattice/rank_two_20261117/README.md) |
| 37 | 3.21356e-7 | 1.85737e-5<br>[Optimized cyclic PSD kernel + forced collinear triples + conditioned mod-2 class bounds — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim37_cyclic_strong/README.md) | **≤1.85604e-5** | 0.07% | 0.07% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 38 | 2.17566e-7 | 1.30587e-5<br>[Cyclic PSD kernel with forced collinear triples + mod-2 class profiles (tight envelopes) — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim38_tight/README.md) | **≤1.27199e-5** | 2.59% | 2.64% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 39 | 1.68028e-7 | 9.06829e-6<br>[Cyclic PSD kernel with forced collinear triples + mod-2 class profiles (tight envelopes) — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim39_tight/README.md) | **≤8.76307e-6** | 3.37% | 3.43% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 40 | 1.47983e-7 | 6.86970e-6<br>[Cyclic PSD kernel with forced collinear triples + mod-2 class profiles (tight envelopes) — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim40_tight/README.md) | **≤6.15990e-6** | 10.33% | 10.56% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 41 | 8.63683e-8 | 4.53126e-6<br>[Reoptimized cyclic PSD kernel + origin-ratio tuning + lattice class profiles — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim41/README.md) | **≤4.49367e-6** | 0.83% | 0.85% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 42 | 5.75168e-8 | 3.53653e-6<br>[Reoptimized cyclic PSD kernel + origin-ratio tuning + lattice class profiles — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim42_own/README.md) | **≤3.52172e-6** | 0.42% | 0.43% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 43 | 4.01571e-8 | 2.56673e-6<br>[Reoptimized cyclic PSD kernel + origin-ratio tuning + lattice class profiles — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim43/README.md) | **≤2.55492e-6** | 0.46% | 0.47% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 44 | 3.20085e-8 | 1.78395e-6<br>[Reoptimized cyclic PSD kernel + origin-ratio tuning + lattice class profiles — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim44/README.md) | **≤1.77966e-6** | 0.24% | 0.25% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 45 | 2.52315e-8 | 1.18594e-6<br>[Optimized cyclic PSD kernel + forced collinear triples + conditioned mod-2 class bounds — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim45_cyclic/README.md) | **≤1.18079e-6** | 0.43% | 0.44% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 46 | 2.27179e-8 | 7.88008e-7<br>[Optimized cyclic PSD kernel + forced collinear triples + conditioned mod-2 class bounds — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim46_cyclic/README.md) | **≤7.85863E-7** | 0.27% | 0.28% | [Nonnegative quadratic Fourier multiplier + scalar rank-two PSD kernel + conditioned mod-2 profiles](verifications/lattice/rank_two_20261117/README.md) |
| 47 | 2.14660e-8 | 5.27910e-7<br>[Optimized cyclic PSD kernel + forced collinear triples + conditioned mod-2 class bounds — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim47_cyclic/README.md) | **≤5.25880e-7** | 0.38% | 0.40% | [Dimension-adapted PSD angular kernels + universal rank-two lattice sums + full collinear-line domination + mod-2 class profiles](verifications/lattice/rank_two_20261117/README.md) |
| 48 | 2.31782e-8 | 3.54688e-7<br>[Optimized cyclic PSD kernel + forced collinear triples + conditioned mod-2 class bounds — Kravatskiy (2026), preceding lattice-only bound](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/verifications/lattice/dim48_cyclic_strong/README.md) | **≤3.52779E-7** | 0.54% | 0.58% | [Additive Fourier sum-of-squares + nonnegative quadratic multiplier + scalar rank-two PSD kernel + conditioned mod-2 profiles](verifications/lattice/rank_two_20261117/README.md) |

Upper bounds rounded **up**, lower bounds **down**, to6 significant digits. [Exact registered values](LATTICE_RESULTS.json). Gain is the reduction from the named upper; Gap is the reduction in its additive gap to the named lattice lower bound.

Lower bounds in38D–42D use the lattice sections of [Dorofeev–Sun–Wang, Theorem1 and Table1](https://arxiv.org/html/2607.20359v4).43D–44D use Conway–Sloane lattice constructions([1982, Corollary8](https://doi.org/10.2307/2007025)); these are not the larger nonlattice entries in the unrestricted table. Source formulas remain in the JSON.

<!-- lattice-upper-20261116:end -->

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

## Packing constructions: lower bounds

**17November update: six independently audited periodic packing improvements in85D–90D.** Gains over the pinned preceding repository constructions are**50.43%,71.27%,47.73%,38.68%,16.86% and4.69%**, respectively. [Exact source reconstruction, proof, source attributions and independent replay](research/packing_release8590_20261117/). Previous constructions remain available. These comparisons do not assert exhaustive worldwide priority.

**Comparison method / source** names the published work behind each recorded
comparison. Values marked **derived** are calculated from published inputs,
rather than quoted as numerical records in those sources. All densities in
this table are occupied-volume fractions.

The[historical literature review dated16November2026](research/covering_and_deep_cosets_20261115/PRIORITY.md) records published guarantees and derived comparisons. The six new85D–90D rows instead compare with their explicitly pinned repository predecessors; that historical review is not a fresh exhaustive priority survey for this release.

| n | Comparison density / method | Our density / method | Increase |
| ---: | --- | --- | ---: |
| 36 | ≈6.1615e-7<br>[KP36 code-based lattice — Kschischang–Pasupathy (1992)](https://doi.org/10.1109/18.119683) · [Catalogue parameters](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KP36.html) | ≈7.3303e-7<br>[Contact-preserving metric deformation; exact Gram matrix and complete minimum certificate](constructions/dim36/) | 18.97% |
| 49 | ≈4.128837170599e-9 (**derived**)<br>P48q × A1 — [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | ≈7.151355755657e-9<br>[Explicit deeper P48q half-coset; full minimum certificate](research/covering_and_deep_cosets_20261115/DEEPER_PROOF.md) | **73.205081%** |
| 50 | ≈8.408155459659e-10 (**derived**)<br>P48q × A2 — [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | ≈2.522446637898e-9<br>[Deep plane + fitted metric; explicit minimum-6 lattice](research/adaptive_layers_20261115/P48_PROOF.md) | **200.000000%** |
| 51 | ≈1.798425511734e-10 (**derived**)<br>P48q × A3 — [Bos–Conway–Sloane (1982), §5](https://neilsloane.com/doc/Me103.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | ≈7.969351595567e-10<br>[Deep plane + fitted metric; explicit minimum-6 lattice](research/adaptive_layers_20261115/P48_PROOF.md) | **343.129368%** |
| 52 | ≈1.094979637686e-10 (**derived**)<br>Two Blichfeldt sections of Leterrier54 — [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0); [derivation](research/range48_96_20261115/verify_range.py) | ≈2.879504124e-10<br>[Deep plane + fitted metric; explicit minimum-6 lattice](research/adaptive_layers_20261115/P48_PROOF.md) | **162.973303%** |
| 53 | ≈1.256125809644e-10 (**derived**)<br>One Blichfeldt section of Leterrier54 — [Leterrier (2022), Corollary 1.4](https://doi.org/10.1007/s40993-022-00321-0); [comparison derivation](research/range48_96_20261115/verify_range.py) | ≈1.496705651057e-10<br>[Primitive dual normal; explicit minimum-10 lattice](research/adaptive_layers_20261115/SECTION53_PROOF.md) | **19.152528%** |
| 57 | ≈8.85119018168e-12 (**derived**)<br>MW56 × A1 — [Schütt–Shioda, Theorem 13.35, p=29](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | ≈1.496124208292e-11<br>[Explicit MW56 deep half-coset; complete minimum-10 certificate](research/creative_sections_20261116/PROOF.md) | **69.030851%** |
| 58 | ≈4.474935526505e-12 (**derived**)<br>Γ72 perpendicular to a rank-14 Leech block — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [exact derived section](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈4.784946339509e-12<br>[Explicit MW56 complex-multiplication plane; full minimum-10 certificate](research/cm58_20261116/) | **6.927716%** |
| 59 | ≈2.585163653787e-12 (**derived**)<br>Γ72 perpendicular to √2κ13 — [Nebe (2012), §4.2](https://arxiv.org/pdf/1008.2862) + [classical κ13 catalogue data](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html); [derivation](constructions/dim59/#literature-derived-comparison) | ≈2.637602590720e-12<br>[Six-translate Ne64 antipode; exact period, shifts and attained minimum](constructions/dim59/) | **2.028457%** |
| 60 | ≈[1.668964977674e-12, 1.680573547119e-12)<br>Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 19.04… | ≈2.395250641682e-12<br>[Explicit Ne64 antipode; exact period and shifts](research/novel_routes_20261113/audits/README.md) | **>42.525785%** |
| 61 | ≈[1.508819950024e-12, 1.519314622713e-12)<br>Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 20.54… | ≈2.022035232382e-12<br>[Explicit Ne64 antipode; exact period and shifts](research/novel_routes_20261113/audits/README.md) | **>33.088644%** |
| 62 | ≈[1.440184355054e-12, 1.450201629426e-12)<br>Elkies section guarantee — [Vardy (1999), §4, p.278](https://doi.org/10.1090/S0002-9947-99-02169-8); published log₂ δ prefix 22.13… | ≈1.628191894214e-12<br>[Explicit Ne64 antipode; exact period and shifts](research/novel_routes_20261113/audits/README.md) | **>12.273483%** |
| 65 | ≈3.201378537714e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 7 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈3.896264824634e-13<br>[Explicit Ne64 extension; complete minimum-12 certificate](research/novel_routes_20261113/higher_parent/) | **21.705846%** |
| 66 | ≈2.272528652057e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 6 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈2.431106576448e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **6.978039%** |
| 67 | ≈1.698310457767e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 5 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈2.215613494306e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **30.459863%** |
| 68 | ≈1.454792126225e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 4 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈2.200710220661e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **51.273174%** |
| 69 | ≈1.23719412351e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 3 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈1.908606667819e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **54.268973%** |
| 70 | ≈1.206263329664e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 2 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈1.73198888183e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **43.582984%** |
| 71 | ≈1.238694909866e-13 (**derived**)<br>Γ72/Leech-block section, normal rank 1 — [Nebe (2012), §4.2](https://arxiv.org/abs/1008.2862); [section derivation](research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py) | ≈1.417772488475e-13<br>[Explicit Γ72 antipode; exact normal, period and shifts](research/novel_routes_20261113/higher_parent/) | **14.456956%** |
| 73 | ≈2.132314974725e-14 (**derived**)<br>Γ72 × A1 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | ≈3.482055772649e-14<br>[Compatible half-cosets; existence proof](research/cross_dimension_lower_20261115/COUNTING_PROOF.md) | **63.299316%** |
| 74 | ≈3.575177737935e-15 (**derived**)<br>Γ72 × A2 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | ≈9.533807301161e-15<br>[Compatible half-cosets; existence proof](research/cross_dimension_lower_20261115/COUNTING_PROOF.md) | **166.666667%** |
| 75 | ≈6.315750546658e-16 (**derived**)<br>Γ72 × A3 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | ≈2.750287321247e-15<br>[Compatible half-cosets; existence proof](research/cross_dimension_lower_20261115/COUNTING_PROOF.md) | **335.464843%** |
| 76 | ≈1.279866139755e-16 (**derived**)<br>Γ72 × D4 — [Nebe (2012)](https://arxiv.org/abs/1008.2862) | ≈9.101270327146e-16<br>[Compatible half-cosets; existence proof](research/cross_dimension_lower_20261115/COUNTING_PROOF.md) | **611.111111%** |
| 77 | ≈9.521458141249e-17 (**derived**)<br>Bare L80 section, normal rank 3 — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf); [normal matrix audit](research/cross_dimension_lower_20261115/l80_results.json) | ≈2.992305608205e-16<br>[Adaptive exceptional coset above D5; lattice existence proof](research/adaptive_layers_20261115/ADAPTIVE_PROOF.md) | **214.269681%** |
| 78 | ≈9.523775961404e-17 (**derived**)<br>L80 perpendicular to 2 A2 — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf); [normal matrix audit](research/cross_dimension_lower_20261115/l80_results.json) | ≈1.366254130898e-16<br>[Explicit MW80 antipode; exact period, shifts and minimum attainment](research/creative_sections_20261116/PROOF.md) | **43.457189%** |
| 79 | ≈9.274742674901e-17 (**derived**)<br>L80 perpendicular to a minimum vector — [Bachoc–Nebe (1998)](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) | ≈1.119208813283e-16<br>[Explicit MW80 antipode; exact period, shifts and minimum attainment](research/creative_sections_20261116/PROOF.md) | **20.672762%** |
| 81 | ≈1.595311593607e-17 (**derived**)<br>MW80 × A1 — [Schütt–Shioda (2017 author draft), Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf); [exact comparison](research/covering_and_deep_cosets_20261115/literature_comparisons.json) | ≈2.791096776318e-17<br>[Cubic covering bound + colored layers; periodic existence proof](research/adaptive_layers_20261115/COVERING_PROOF.md) | **74.956215%** |
| 82 | ≈2.541814578499e-18 (**derived**)<br>MW80 × A2 — [Schütt–Shioda, Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | ≈7.514357550704e-18<br>[Cubic covering bound + colored layers; periodic existence proof](research/adaptive_layers_20261115/COVERING_PROOF.md) | **195.629650%** |
| 83 | ≈4.269748590531e-19 (**derived**)<br>MW80 × A3 — [Schütt–Shioda, Theorem 13.35, p=41](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | ≈2.204491680495e-18<br>[Compatible half-cosets; existence proof](research/cross_dimension_lower_20261115/COUNTING_PROOF.md) | **416.304797%** |
| 84 | ≈2.970203932666e-19 (**derived**)<br>Section of BCS88 perpendicular to √2 D4 — [Bos–Conway–Sloane (1982), Theorems 1–2, Eq. (15)](https://neilsloane.com/doc/Me103.pdf); [derivation and exact matrix audit](research/priority_resolution_20261115/CLASSICAL_SECTIONS.md) | ≈7.603722067624e-19<br>[Compatible half-cosets above L80; existence proof](research/cross_dimension_lower_20261115/COUNTING_PROOF.md#4-the-80d-parent-gives-deeper-nonsingular-half-cosets) | **156%** |
| 85 | ≈2.341415812047e-19<br>[Adaptive exceptional-coset lattice — Kravatskiy (2026)](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/research/adaptive_layers_20261115/ADAPTIVE_PROOF.md) | **≥3.52219521096589e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **50.43%** |
| 86 | ≈1.661031863043e-19<br>[Preceding MW92 simplex antipode — Kravatskiy (2026)](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/research/simplex_normals_20261116/PROOF.md) | **≥2.84477266080525e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **71.27%** |
| 87 | ≈1.746637359507e-19<br>[Preceding MW92 simplex antipode — Kravatskiy (2026)](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/research/simplex_normals_20261116/PROOF.md) | **≥2.58037767092447e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **47.73%** |
| 88 | ≈1.746868055945e-19<br>[Preceding MW92 simplex antipode — Kravatskiy (2026)](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/research/simplex_normals_20261116/PROOF.md) | **≥2.42250297455694e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **38.68%** |
| 89 | ≈1.697792690692e-19<br>[Preceding MW92 simplex antipode — Kravatskiy (2026)](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/research/simplex_normals_20261116/PROOF.md) | **≥1.98411302345717e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **16.86%** |
| 90 | ≈1.543460015256e-19<br>[Preceding MW92 simplex antipode — Kravatskiy (2026)](https://github.com/alexlegeartis/SpherePackings/blob/20b638f5d5b0bf00a2a4ce8d70c4e062517a502e/research/mw92_sections_20261115/PROOF.md) | **≥1.61584267273781e-19**<br>[Correlated Mordell–Weil dual-normal sections with optimized simplex translate clusters](research/packing_release8590_20261117/PROOF.md) | **4.69%** |
| 91 | ≈[9.598621803784e-20, 9.59946031323e-20] (**derived**)<br>MW92 / full Rogers-simplex section — [Schütt–Shioda, Theorem 13.35](<https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf>) + [Rogers (1958), via Fejes Tóth §2.2](https://arxiv.org/abs/2202.11358); [derived interval](research/mw92_sections_20261115/PROOF.md#external-comparisons-including-the-full-rogers-consequence) | ≈1.27810250822e-19<br>[Explicit MW92 antipode; complete period and 2 rational shifts](research/mw92_sections_20261115/PROOF.md) | **>33.143162%** |
| 93 | ≈1.702354782363e-20 (**derived**)<br>MW92 × A1 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | ≈3.045563378489e-20<br>[Cubic covering bound + colored layers; periodic existence proof](research/adaptive_layers_20261115/COVERING_PROOF.md) | **78.902977%** |
| 94 | ≈2.534314091589e-21 (**derived**)<br>MW92 × A2 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | ≈7.710930742674e-21<br>[Cubic covering bound + colored layers; periodic existence proof](research/adaptive_layers_20261115/COVERING_PROOF.md) | **204.261053%** |
| 95 | ≈3.980717377502e-22 (**derived**)<br>MW92 × A3 — [Schütt–Shioda, Theorem 13.35, p=47](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf) | ≈2.043145398811e-21<br>[Cubic covering bound + colored layers; periodic existence proof](research/adaptive_layers_20261115/COVERING_PROOF.md) | **413.260602%** |
| 96 | ≈2.795023294827e-22<br>P48q-based lattice — [Bos–Conway–Sloane (1982), Eq. (15), Table 1, p.177](https://neilsloane.com/doc/Me103.pdf); exact center density `3^48/2^24` | ≈6.111066241809e-22<br>[Cubic covering bound + colored layers; periodic existence proof](research/adaptive_layers_20261115/COVERING_PROOF.md) | **118.640977%** |

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

## Registered results

**Recorded upper U** is the historical comparison used for the percentages in that row. **Comparison method / source** names the method behind U and links to its paper or pinned certificate; **Our method / proof** describes the new bound N. A recorded comparison is not a claim about the latest worldwide record. [Full source titles, versions and comparison scope](COMPARISONS.md#table-comparison-methods-and-sources).

| n / Lower L | Recorded U / comparison method | Our upper N / method | Upper decrease | Gap decrease |
| ---: | --- | --- | --- | ---: |
| 4<br>**L:** [≈0.61685](https://cohn.mit.edu/sphere-packing/) | 0.629521681653<br>[Scalar/kernel + count and far-pair refinements — Felpix-Studios, pinned bb97550](https://github.com/Felpix-Studios/sphere-packing/tree/bb97550ba47ad92bacfee567843888be16895a23/preprints/D4) | <0.629520189662<br>[Degree16 far-pair polynomial + unused angular Gram credits](verifications/dim4/refinement_energy20261115/) · [Sources & formulas](COMPARISONS.md) | 0.000237% | 0.011774% |
| 5<br>**L:** [≈0.46526](https://cohn.mit.edu/sphere-packing/) | ≈0.51226<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.50995<br>[Spherical lens caps + long-edge PSD accounting](verifications/dim5/refinement_spherical_lens/) · [Sources & formulas](COMPARISONS.md) | 0.45% | 4.92% |
| 6<br>**L:** [≈0.37295](https://cohn.mit.edu/sphere-packing/) | ≈0.41029<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.40913<br>[Neighbor-count cuts + Schur-saturated cyclic kernel](verifications/dim6/refinement_patch5/) · [Sources & formulas](COMPARISONS.md) | 0.28% | 3.11% |
| 7<br>**L:** [≈0.29530](https://cohn.mit.edu/sphere-packing/) | ≈0.32111<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.32058<br>[Neighbor-count cuts + Schur-saturated cyclic kernel](verifications/dim7/) · [Sources & formulas](COMPARISONS.md) | 0.17% | 2.05% |
| 9<br>**L:** [≈0.14577](https://cohn.mit.edu/sphere-packing/) | ≈0.19099<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.19093<br>[Signed-hinge cyclic profile + average-neighbor count](verifications/dim9/) · [Sources & formulas](COMPARISONS.md) | 0.03% | 0.14% |
| 10<br>**L:** [≈0.099616](https://cohn.mit.edu/sphere-packing/) | ≈0.14340<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.14333<br>[Log-hinge cyclic profile + root-cut radial repair](verifications/dim10/) · [Sources & formulas](COMPARISONS.md) | 0.05% | 0.17% |
| 11<br>**L:** [≈0.066238](https://cohn.mit.edu/sphere-packing/) | ≈0.10672<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.10664<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim11/) · [Sources & formulas](COMPARISONS.md) | 0.08% | 0.21% |
| 12<br>**L:** [≈0.049454](https://cohn.mit.edu/sphere-packing/) | ≈0.079707<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.079503<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim12/) · [Sources & formulas](COMPARISONS.md) | 0.26% | 0.67% |
| 13<br>**L:** [≈0.032014](https://cohn.mit.edu/sphere-packing/) | ≈0.060121<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.059905<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim13/) · [Sources & formulas](COMPARISONS.md) | 0.36% | 0.77% |
| 14<br>**L:** [≈0.021624](https://cohn.mit.edu/sphere-packing/) | ≈0.045010<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.045003<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim14/) · [Sources & formulas](COMPARISONS.md) | 0.02% | 0.03% |
| 15<br>**L:** [≈0.016858](https://cohn.mit.edu/sphere-packing/) | ≈0.033726<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.033672<br>[Quartic radial damping + all-dilation cyclic reduction](verifications/dim15/) · [Sources & formulas](COMPARISONS.md) | 0.16% | 0.32% |
| 16<br>**L:** [≈0.014708](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/BW16.html) | ≈0.024921<br>[Three-point SDP — Tutunov et al. (2025), Table 1](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.024901<br>[Widened shell count + coupled log-hinge radial kernel](verifications/dim16/widened_reconstruction/) · [Sources & formulas](COMPARISONS.md) | 0.08% | 0.20% |
| 17<br>**L:** [≈0.0088113](https://doi.org/10.4153/CJM-1964-065-1) | ≈0.018464<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.018140<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim17/refinement_larger20261115/) · [Sources & formulas](COMPARISONS.md) | 1.75% | 3.36% |
| 18<br>**L:** [≈0.0061679](https://doi.org/10.1023/A:1008733715001) | ≈0.013485<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.013077<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim18/refinement_larger20261115/) · [Sources & formulas](COMPARISONS.md) | 3.02% | 5.57% |
| 19<br>**L:** [≈0.0041474](https://arxiv.org/abs/2505.02394v1) | ≈0.0098180<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0093654<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim19/refinement_larger20261115/) · [Sources & formulas](COMPARISONS.md) | 4.61% | 7.98% |
| 20<br>**L:** [≈0.0040242](https://arxiv.org/abs/2505.02394v1) | ≈0.0071271<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0067426<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim20/refinement_larger20261115/) · [Sources & formulas](COMPARISONS.md) | 5.39% | 12.39% |
| 21<br>**L:** [≈0.0029298](https://arxiv.org/abs/2505.02394v1) | ≈0.0051597<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0049104<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim21/refinement_larger20261115/) · [Sources & formulas](COMPARISONS.md) | 4.83% | 11.18% |
| 22<br>**L:** [≈0.0024510](https://cohn.mit.edu/sphere-packing/) | ≈0.0037259<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0036111<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim22/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.08% | 9.01% |
| 23<br>**L:** [≈0.0019072](https://cohn.mit.edu/sphere-packing/) | ≈0.0026843<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0026135<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim23/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 2.64% | 9.10% |
| 25<br>**L:** [≈0.00067721](https://cohn.mit.edu/sphere-packing/) | ≈0.0013842<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0013611<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim25/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 1.67% | 3.27% |
| 26<br>**L:** [≈0.00026922](https://cohn.mit.edu/sphere-packing/) | ≈0.00099102<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00096119<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim26/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 3.01% | 4.13% |
| 27<br>**L:** [≈0.00015759](https://cohn.mit.edu/sphere-packing/) | ≈0.00070823<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00067374<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim27/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 4.87% | 6.26% |
| 28<br>**L:** [≈0.00010464](https://cohn.mit.edu/sphere-packing/) | ≈0.00050525<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00047150<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim28/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 6.68% | 8.43% |
| 29<br>**L:** [≈0.000034145](https://cohn.mit.edu/sphere-packing/) | ≈0.00035986<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00033305<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim29/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 7.45% | 8.23% |
| 30<br>**L:** [≈0.000021915](https://cohn.mit.edu/sphere-packing/) | ≈0.00025590<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00023903<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim30/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 6.59% | 7.21% |
| 31<br>**L:** [≈0.000011838](https://cohn.mit.edu/sphere-packing/) | ≈0.00018171<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00017253<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim31/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 5.05% | 5.40% |
| 32<br>**L:** [≈0.000011041](https://cohn.mit.edu/sphere-packing/) | ≈0.00012884<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.00012276<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim32/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 4.72% | 5.16% |
| 33<br>**L:** [≈0.0000041407](https://cohn.mit.edu/sphere-packing/) | ≈0.000091236<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000089087<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim33/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 2.35% | 2.47% |
| 34<br>**L:** [≈0.0000017670](https://cohn.mit.edu/sphere-packing/) | ≈0.000064522<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000063254<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim34/refinement_combined20261115/) · [Sources & formulas](COMPARISONS.md) | 1.97% | 2.02% |
| 35<br>**L:** [≈9.4619e-7](https://cohn.mit.edu/sphere-packing/) | ≈0.000045574<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000043911<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim35/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.65% | 3.73% |
| 36<br>**L:** [≈7.3303e-7 (this work)](constructions/dim36/) | ≈0.000032153<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000030495<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim36/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 5.16% | 5.28% |
| 37<br>**L:** [≈3.2136e-7](https://cohn.mit.edu/sphere-packing/) | ≈0.000022659<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000020943<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim37/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 7.57% | 7.68% |
| 38<br>**L:** [≈2.1757e-7](https://arxiv.org/html/2607.20359v4) | ≈0.000015951<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000014518<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim38/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 8.98% | 9.10% |
| 39<br>**L:** [≈1.6803e-7](https://arxiv.org/html/2607.20359v4) | ≈0.000011217<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.000010220<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim39/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 8.89% | 9.02% |
| 40<br>**L:** [≈1.4798e-7](https://arxiv.org/html/2607.20359v4) | ≈0.0000078801<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000073057<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim40/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 7.29% | 7.43% |
| 41<br>**L:** [≈8.6368e-8](https://arxiv.org/html/2607.20359v4) | ≈0.0000055306<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000051253<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim41/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 7.33% | 7.45% |
| 42<br>**L:** [≈5.7517e-8](https://arxiv.org/html/2607.20359v4) | ≈0.0000038781<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000036196<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim42/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 6.67% | 6.77% |
| 43<br>**L:** [≈4.6423e-8](https://arxiv.org/html/2607.20359v4) | ≈0.0000027169<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000026096<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim43/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.95% | 4.02% |
| 44<br>**L:** [≈3.9244e-8](https://arxiv.org/abs/2505.02394v1) | ≈0.0000019018<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000018326<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim44/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 3.64% | 3.71% |
| 45<br>**L:** [≈3.5582e-8](https://arxiv.org/abs/2505.02394v1) | ≈0.0000013301<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈0.0000012656<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim45/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 4.85% | 4.98% |
| 46<br>**L:** [≈2.8610e-8](https://doi.org/10.1007/BF01232378) | ≈9.2952e-7<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈8.7124e-7<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim46/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 6.27% | 6.47% |
| 47<br>**L:** [≈2.2670e-8](https://arxiv.org/abs/2505.02394v1) | ≈6.4908e-7<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈5.9377e-7<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim47/refinement_profile20261114/) · [Sources & formulas](COMPARISONS.md) | 8.52% | 8.83% |
| 48<br>**L:** [≈2.3178e-8](https://doi.org/10.4153/CJM-1971-081-3) | ≈4.5291e-7<br>[Cohn–Elkies LP / modular bootstrap — Afkhami-Jeddi et al. (2020)](https://arxiv.org/abs/2006.02560v2) | ≈4.0844e-7<br>[Signed logarithmic profile + degree95 radial enrichment](verifications/dim48/refinement_larger20261115/) · [Sources & formulas](COMPARISONS.md) | 9.82% | 10.35% |

The 4D comparison uses the pinned Felpix-Studios certificate and its complete full-mode replay; extra digits show the small strict improvement. The 16D strict bound is 0.024900975, compared with the paper’s displayed 0.02492121. The full widened-shell count certificate is included and independently replayed. No worldwide-priority claim is made. Other densities are shown approximately to five significant digits and percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. The 5D, 6D, 7D and 9–16D baselines are pinned paper-reported decimal values, not independently replayed here.

These certificates use the existing [Cohn–de Laat–Salmon continuum framework](https://arxiv.org/abs/2206.15373). Proof notes state the analytic hypotheses and normalization; exact rational witnesses and independent interval verifiers check the computational conditions over full domains. Historical numerical comparisons in package notes are pinned comparisons, not current best-known claims.

## Reproduce

For the upper-bound suite, Python 3.12 or newer and a C++20 compiler are required by the pinned 4D proof. Install dependencies and run:

    python -m pip install -r requirements.txt
    python run_all.py

Never disable assertions (`python -O`). No optimization solver or network is needed after dependency installation. The 48D compact triangle replay may take several minutes. Logs record the original audits; fresh replay is recommended. No GitHub Actions workflow is installed.

The separate 36D lower-bound replay is `python constructions/dim36/verify_all.py`; it requires SymPy and a C++17 compiler and performs the complete minimum enumeration. It is not silently included in the historical upper-bound suite.

## Layout

- `verifications/dim17/`: solid-harmonic cyclic witness, explicit finite-sum proof and independent exact-domain audit.
- `verifications/dim18/`: exact degree-31 continuum witness, independent audit and two full-domain verification paths.
- `verifications/dim5/`: count-enhanced cyclic witness, proof, independent audit and global replay.

- `verifications/dim4/`: cyclic kernel, published neighbor-count inputs, exact witness and global verification.
- `verifications/dim19/`: exact witness, proof and two verification paths.
- `verifications/dim20/`: exact witness, proof and two verification paths.
- `verifications/dim21/`: exact witness, proof and two verification paths.
- `verifications/dim48/`: current cyclic witness, complete proof dependencies and independent angular, triangle-domain and radial audits.
- `RESULTS.json`: machine-readable registry.
- `run_all.py`: full certificate replay, failing on any unsuccessful subprocess.

The referenced continuum theorem and standard positive-definite Gegenbauer addition theorem are mathematical inputs, not independently formalized theorems in this repository.

The [origin-coordinate refinements](verifications/origin_slack/) tighten the existing 4D, 5D, 19D, 20D and 21D certificates. All earlier witnesses remain available; full source replays and independent exact Schur/Arb audits passed. Dimension 18 is newly registered separately below.

## Research tools

The [research index](research/README.md) separates registered certificates,
historical results and unsuccessful searches. The [16 November branch
reconciliation](research/branch_integration_20261116/README.md) records where
the three research branches were incorporated and which results were superseded.

- [Dimension-generic annular three-point search](research/annular_three_point/): global-positivity SDP discovery, exact export and full-domain checks. The degree-31 18D witness is now registered after independent audit. Degree-23 17D/18D and degree-31 17D certificates are weaker than published bounds; 22D/23D candidates are rejected. All earlier experiments remain preserved.
- [Origin-coordinate slack research](research/origin_slack/): exact Schur-complement postprocessing and independent 5D audit. This preserves the initial research snapshot; the subsequent fully replayed refinements are registered in [verifications/origin_slack](verifications/origin_slack/).

Research dependencies are separate from the verification-only root requirements. See each package for reproducibility commands, tested versions, scope and checksums.

The newly audited [18D certificate](verifications/dim18/) proves density <0.013408931, improving the pinned published upper by 0.566611%. It uses the continuum theorem without external neighbor-count bounds. Comparison freshness and worldwide-priority caveats remain in force.

The [17D solid-harmonic cyclic certificate](verifications/dim17/) proves density <0.018335552. Its proof uses an explicit cyclic generalization and periodic PSD summation, not the individual-sign hypothesis of the printed Theorem 1.1. The [search code and rejected experiments](research/solid_cyclic/) are retained separately.

The [18D cyclic refinement](verifications/dim18/refinement_cyclic/) strengthens the earlier bound to <0.013284728. Both cyclic packages retain explicit summation proofs and independent exact-domain audits; earlier witnesses remain available.

The [19D cyclic refinement](verifications/dim19/refinement_cyclic/) proves density <0.009557491 and the [20D cyclic refinement](verifications/dim20/refinement_cyclic/) proves <0.006977702. Both use the explicitly proved cyclic generalization and independently reconstructed exact-domain checks. All older witnesses and Schur certificates remain replayable.

The current 17D bound <0.018204096 uses a [mixed-Gaussian repaired cyclic witness and separate saturation](verifications/dim17/refinement_repaired/). The [common cyclic origin saturation](verifications/cyclic_origin_slack/) sharpens 18D, 19D and 20D to <0.013224473, <0.009521000 and <0.006936042. All prior witnesses remain preserved and replayed. See [methods and experiments](METHODS_AND_EXPERIMENTS.md) for the systematic distinction between valid proof transformations and unverified search directions.

The [47D restriction and origin-Gram certificate](verifications/dim47/) gives <6.16879891e-7 by transferring the fully replayed 48D witness. The pinned upper decreases 4.96%, and the additive gap to the cited 2025 antipode construction decreases 5.14%. The proof and independent verifier explicitly preserve the cyclic theorem and Fourier restriction hypotheses.

The [6D count-enhanced certificate](verifications/dim6/) gives <0.410189987, using externally published near-neighbor bounds. The [45D](verifications/dim45/) and [46D](verifications/dim46/) exact radial/angular certificates give <1.311423575e-6 and <9.14882847e-7. All new proof chains and preserved historical packages are included in the aggregate replay.

The [portable 6D discovery/export package](research/dim6_discovery/) reproduces the successful bounded warm-start search and exports its frozen candidate for full exact verification.

The [improved 6D certificate](verifications/dim6/refinement_mixture99/) strengthens the earlier result to <0.409135399, closing 3.09% of the pinned additive gap. Its physical and Fourier signs hold without Gaussian repair; the same published neighbor-count dependencies remain. The earlier 6D proof is preserved.

The [6D refinement discovery package](research/dim6_refinement_discovery/) reconstructs the frozen 99/100 source exactly from preserved checkpoints and exports it for full verification.

The [44D Gaussian cyclic certificate](verifications/dim44/) gives <1.842868599e-6, about 3.10% below the pinned upper. The [small 6D follow-on](verifications/dim6/refinement_patch5/) gives <0.409128533, with no Gaussian repair and unchanged external count dependencies.

The [7D count-enhanced certificate](verifications/dim7/) gives <0.320584772, closing 2.05% of the pinned additive gap to E7. Its external near-neighbor count inputs are cited explicitly.

The [guarded 5D refinement](verifications/dim5/refinement_guarded/) gives <0.510129239; its unguarded checkpoint is invalid. The [22D tetrahedral aggregate research](research/tetrahedral_aggregate/) proves stronger auxiliary edge/triangle inequalities but makes no density-record claim.

Portable exact-source reconstruction and export companions are available for [guarded 5D](research/dim5_guarded_discovery/), [7D](research/dim7_discovery/) and the [small 6D follow-on](research/dim6_patch5_discovery/). Numerical checkpoints remain distinct from certified witnesses.

The [5D long-edge lens extension](verifications/dim5/refinement_lens/) gives <0.509996147, closing 4.83% of the pinned additive gap. It uses newly audited long-edge PSD accounting, not the old pointwise far-sign condition. Earlier proof packages remain preserved.

The [stronger 5D spherical-lens follow-on](verifications/dim5/refinement_spherical_lens/) gives <0.509953498, closing 4.92% of the pinned additive gap. Its exact transverse spherical-cap proofs and new tail budgets are replayed; the centroid-only version remains available. Portable search/export companions are included for [centroid](research/dim5_lens_discovery/) and [spherical](research/dim5_spherical_lens_discovery/) lenses.

The [22D](verifications/dim22/), [23D](verifications/dim23/), [26D](verifications/dim26/) and [27D](verifications/dim27/) Gaussian cyclic certificates have fresh primary and independent proof replays. [Release checks](research/gaussian_22_27_release/) distinguish those from reused historical audits. Portable discovery and exact export are included for [22/23](research/gaussian22_23_discovery/), [26](research/dim26_gaussian_discovery/) and [27](research/dim27_gaussian_discovery/).

The remaining [29D](verifications/dim29/), [30D](verifications/dim30/), [31D](verifications/dim31/), [33D](verifications/dim33/), [35D](verifications/dim35/) and [37D](verifications/dim37/) entries complete the Gaussian cyclic certificates in dimensions25–40. Their [release checks](research/gaussian_remaining_release/) distinguish fresh primary/export replays, complete independent audits and reused unchanged historical proofs.

The [21D Gaussian cyclic refinement](verifications/dim21/refinement_gaussian/) proves density<0.004985977000, retaining all earlier21D proofs. Its [release evidence](research/gaussian21_release/) distinguishes fresh proof checks from unchanged historical audits.

The [21D discovery companion](research/dim21_gaussian_discovery/) includes the original seven-solve search history, exact inputs and deterministic export replay; numerical checkpoints are not proof certificates.

The [46D](verifications/dim46/refinement_angular/) and [47D](verifications/dim47/refinement_angular/) angular refinements lower the origin cost with rebuilt target-dimensional angular kernels and unchanged radial witnesses. Earlier proofs are preserved; [release checks](research/angular_border_release/) document fresh new and historical46/47 replays and reuse of other unchanged audits.

The subsequent [46D](verifications/dim46/refinement_joint/) and [47D](verifications/dim47/refinement_joint/) joint refits reoptimize the radial Fourier pair for the rebuilt angular kernel, improving both bounds further. [Release evidence](research/angular_radial_refit_release/) includes fresh primary/independent replays, exact export checks and explicit historical-check scope.

## Citation

If you use these bounds or verification packages, please cite:

Alexey Kravatskiy. (2026). *Sphere-packing density: new upper bounds in dimensions 4 through 48, with verification packages*. https://github.com/alexlegeartis/SpherePackingUpperBounds

[CITATION.cff](CITATION.cff) supplies GitHub's citation metadata; [CITATION.bib](CITATION.bib) provides BibTeX. For a reproducible reference to particular results, also record the commit used.

