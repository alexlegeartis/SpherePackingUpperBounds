# Sphere-packing density upper bounds

Reproducible certificates for unrestricted sphere packings, organized by dimension on the model of KissingNumbers. These are upper bounds on density, not packing constructions, kissing-number results, or optimality claims. No worldwide-priority claim is made.

## Registered results

| Dimension | Known lower density L | Older upper U | New upper N | Upper decrease | Additive gap decrease | Proof and details |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | [≈0.61685](https://cohn.mit.edu/sphere-packing/) | [≈0.63610](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.63343 | 0.42% | 13.89% | [4D proof](verifications/origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 5 | [≈0.46526](https://cohn.mit.edu/sphere-packing/) | [≈0.51226](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.50995 | 0.45% | 4.92% | [5D proof](verifications/dim5/refinement_spherical_lens/) · [Sources & formulas](COMPARISONS.md) |
| 6 | [≈0.37295](https://cohn.mit.edu/sphere-packing/) | [≈0.41029](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.40913 | 0.28% | 3.11% | [6D proof](verifications/dim6/refinement_patch5/) · [Sources & formulas](COMPARISONS.md) |
| 7 | [≈0.29530](https://cohn.mit.edu/sphere-packing/) | [≈0.32111](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.32058 | 0.17% | 2.05% | [7D proof](verifications/dim7/) · [Sources & formulas](COMPARISONS.md) |
| 17 | [≈0.0088113](https://doi.org/10.4153/CJM-1964-065-1) | [≈0.018464](https://arxiv.org/abs/2006.02560v2) | ≈0.018204 | 1.41% | 2.69% | [17D proof](verifications/dim17/refinement_repaired/) · [Sources & formulas](COMPARISONS.md) |
| 18 | [≈0.0061679](https://doi.org/10.1023/A:1008733715001) | [≈0.013485](https://arxiv.org/abs/2006.02560v2) | ≈0.013224 | 1.93% | 3.57% | [18D proof](verifications/cyclic_origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 19 | [≈0.0041474](https://arxiv.org/abs/2505.02394v1) | [≈0.0098180](https://arxiv.org/abs/2006.02560v2) | ≈0.0095210 | 3.02% | 5.24% | [19D proof](verifications/cyclic_origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 20 | [≈0.0040242](https://arxiv.org/abs/2505.02394v1) | [≈0.0071271](https://arxiv.org/abs/2006.02560v2) | ≈0.0069360 | 2.68% | 6.16% | [20D proof](verifications/cyclic_origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 21 | [≈0.0029298](https://arxiv.org/abs/2505.02394v1) | [≈0.0051597](https://arxiv.org/abs/2006.02560v2) | ≈0.0051495 | 0.20% | 0.45% | [21D proof](verifications/origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 44 | [≈3.9244e-8](https://arxiv.org/abs/2505.02394v1) | [≈0.0000019018](https://arxiv.org/abs/2006.02560v2) | ≈0.0000018429 | 3.10% | 3.16% | [44D proof](verifications/dim44/) · [Sources & formulas](COMPARISONS.md) |
| 45 | [≈3.5582e-8](https://arxiv.org/abs/2505.02394v1) | [≈0.0000013301](https://arxiv.org/abs/2006.02560v2) | ≈0.0000013114 | 1.40% | 1.44% | [45D proof](verifications/dim45/) · [Sources & formulas](COMPARISONS.md) |
| 46 | [≈2.8610e-8](https://doi.org/10.1007/BF01232378) | [≈9.2952e-7](https://arxiv.org/abs/2006.02560v2) | ≈9.1488e-7 | 1.57% | 1.62% | [46D proof](verifications/dim46/) · [Sources & formulas](COMPARISONS.md) |
| 47 | [≈2.2670e-8](https://arxiv.org/abs/2505.02394v1) | [≈6.4908e-7](https://arxiv.org/abs/2006.02560v2) | ≈6.1688e-7 | 4.96% | 5.14% | [47D proof](verifications/dim47/) · [Sources & formulas](COMPARISONS.md) |
| 48 | [≈2.3178e-8](https://doi.org/10.4153/CJM-1971-081-3) | [≈4.5291e-7](https://arxiv.org/abs/2006.02560v2) | ≈4.1339e-7 | 8.73% | 9.20% | [48D proof](verifications/dim48/) · [Sources & formulas](COMPARISONS.md) |

Densities shown approximately to five significant digits; percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. The 4D, 5D, 6D and 7D baselines are pinned paper-reported decimal values, not independently replayed here. 

These certificates use the existing [Cohn–de Laat–Salmon continuum framework](https://arxiv.org/abs/2206.15373). Proof notes state the analytic hypotheses and normalization; exact rational witnesses and independent interval verifiers check the computational conditions over full domains. Historical numerical comparisons in package notes are pinned comparisons, not current best-known claims.

## Reproduce

Python 3.11 or newer is required. Install dependencies and run:

    python -m pip install -r requirements.txt
    python run_all.py

Never disable assertions (`python -O`). No optimization solver or network is needed after dependency installation. The 48D compact triangle replay may take several minutes. Logs record the original audits; fresh replay is recommended. No GitHub Actions workflow is installed.

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
