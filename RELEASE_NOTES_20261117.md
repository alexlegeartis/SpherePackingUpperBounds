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
