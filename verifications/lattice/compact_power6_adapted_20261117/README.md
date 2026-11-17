# Dimension-adapted three-feature compact harmonic lattice bounds

Complete producer proofs and independent audits establish the following conservative **lattice-only** sphere-packing density upper bounds:

| Dimension | Upper bound | Improvement over preceding repository bound |
|---:|---:|---:|
|43|<2.5363822051e-6|0.7255%|
|44|<1.7676809433e-6|0.4991%|

The method uses compact radial powers6,9,12 with positive-semidefinite harmonic matrices, optimized separately for each ambient dimension. Universal rank-two lattice-plane inequalities and exact Fourier closure are verified independently for each dimension.43D has zero annular budget: its final proof needs no mod-2 profile count.44D's small positive residual budget uses verified mod-2 annular bounds. These are not unrestricted-packing results or optimality claims.

## Reproduce

Install `requirements.txt`, then run `python run_all.py`. The full runner verifies all hashes, materializes inputs in a temporary directory, and runs the independent exact algebra, cap and384-bit complete geometry checks, followed by fresh scalar annulus and tail checks and exact normalization. `--dimension 43` selects one dimension; `--jobs 1` runs serially instead of the default two workers. `--quick` checks exact algebra and caps only and is **not** a full geometry replay. Assertions must remain enabled. No optimizer, network or unpublished file is required.

43D has212423 geometry nodes/106212 leaves and653 scalar annuli;44D has175733 nodes/87867 leaves and633 annuli. Each source includes153 independently checked Bernstein/Sturm cap inequalities. See the per-dimension `THEOREM_REVIEW.md` and `REPORT.json` for audit scope. Exact-only checker messages saying geometry is still required describe that individual checker's scope; the full runner additionally supplies the complete geometry proof.

Large JSON inputs are losslessly compressed. `STORED_FILES.json` authenticates stored bytes and `MANIFEST.json` authenticates decompressed logical files before execution. Independent checks share Arb arithmetic and reviewed Jet interval evaluation with the producer; no independent arithmetic backend is claimed.

The positive-definite/Fourier framework is related to [Cohn–de Laat–Salmon](https://arxiv.org/abs/2206.15373); the class-profile comparison and annular context are [Zubrilina's lattice bounds](https://arxiv.org/abs/2610.35830). The specific compact-kernel and rank-two arguments are stated in the audit reviews. External benchmark attributions remain separate from the internal improvements above.
