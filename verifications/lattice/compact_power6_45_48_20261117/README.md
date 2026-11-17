# Optimized compact-power harmonic matrices and exact dimensional restrictions

**Verified lattice-only upper bounds.** The producer proof and independent 384-bit geometry replay cover all 183451 nodes and 91726 leaves with no unresolved boxes.

Conservative lattice-only upper bounds:

| Dimension | Conservative upper |
| --- | --- |
|45|1.1728551765e-6|
|46|7.721596228e-7|
|47|5.130033534e-7|
|48|3.4559190288e-7|

Method: three-feature compact harmonic positive-semidefinite matrices with optimized radial powers6,9,12, universal rank-two lattice-plane inequalities, and exact Fourier/mod-2 annular closure. Dimensions45–47 use exact target-dimensional harmonic expansions of the48D kernel and separately reoptimized, independently checked scalar Fourier pairs.

Install `requirements.txt` and run `python run_all.py`. Full replay checks all immutable hashes, exact PSD/origin and harmonic restriction identities, exact Fourier positive factors and Laguerre involution, the entire main geometry trace at384bits,153 exact Bernstein/Sturm cap inequalities, all physical tails, annular budgets and exact volume rounding. The runner works in a temporary copy. No optimizer, external network, or unpublished file is required. Assertions must remain enabled. `--dimension 47` selects a target scalar replay while retaining its shared source geometry. `--quick` checks exact algebra/caps only and is not a full certificate replay.

The48D geometry has compact support19/13 in squared radius; only the three reduced primitive lines of a lattice plane can contribute. The proof generalizes the earlier two-feature kernel to three radial features. These are lattice-only results, not unrestricted sphere-packing upper bounds or optimality claims. The stored frozen input's pending labels record its freeze-time state; the complete replay and final audit reports determine certification status. Exact-only checker messages saying source geometry is still required describe that checker’s scope; the full runner supplies this geometry proof.

Underlying framework: [Zubrilina's lattice class-profile method](https://arxiv.org/abs/2610.35830), and [Cohn–de Laat–Salmon's positive-definite/Fourier framework](https://arxiv.org/abs/2206.15373). See `source48/THEOREM_REVIEW.md` for the specific rank-two and collinear derivation and audit scope. Independent checks share Arb arithmetic with the producer; no independent interval backend is claimed.

Large JSON inputs are losslessly gzip-compressed. `STORED_FILES.json` authenticates stored bytes and `MANIFEST.json` authenticates every decompressed logical file; the runner checks both before execution.
