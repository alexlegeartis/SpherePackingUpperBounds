# Sphere-packing density upper bounds

Reproducible certificates for unrestricted sphere packings, organized by dimension on the model of KissingNumbers. These are upper bounds on density, not packing constructions, kissing-number results, or optimality claims. No worldwide-priority claim is made.

## Registered results

| Dimension | Known lower density L | Older upper U | New upper N | Upper decrease | Additive gap decrease | Proof and details |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | [≈0.61685](https://cohn.mit.edu/sphere-packing/) | [≈0.63610](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.63343 | 0.42% | 13.89% | [4D proof](verifications/origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 5 | [≈0.46526](https://cohn.mit.edu/sphere-packing/) | [≈0.51226](https://arxiv.org/html/2512.04829v2#S3.T1) | ≈0.51050 | 0.34% | 3.75% | [5D proof](verifications/origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 17 | [≈0.0088113](https://doi.org/10.4153/CJM-1964-065-1) | [≈0.018464](https://arxiv.org/abs/2006.02560v2) | ≈0.018204 | 1.41% | 2.69% | [17D proof](verifications/dim17/refinement_repaired/) · [Sources & formulas](COMPARISONS.md) |
| 18 | [≈0.0061679](https://doi.org/10.1023/A:1008733715001) | [≈0.013485](https://arxiv.org/abs/2006.02560v2) | ≈0.013224 | 1.93% | 3.57% | [18D proof](verifications/cyclic_origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 19 | [≈0.0041474](https://arxiv.org/abs/2505.02394v1) | [≈0.0098180](https://arxiv.org/abs/2006.02560v2) | ≈0.0095210 | 3.02% | 5.24% | [19D proof](verifications/cyclic_origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 20 | [≈0.0040242](https://arxiv.org/abs/2505.02394v1) | [≈0.0071271](https://arxiv.org/abs/2006.02560v2) | ≈0.0069360 | 2.68% | 6.16% | [20D proof](verifications/cyclic_origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 21 | [≈0.0029298](https://arxiv.org/abs/2505.02394v1) | [≈0.0051597](https://arxiv.org/abs/2006.02560v2) | ≈0.0051495 | 0.20% | 0.45% | [21D proof](verifications/origin_slack/) · [Sources & formulas](COMPARISONS.md) |
| 48 | [≈2.3178e-8](https://doi.org/10.4153/CJM-1971-081-3) | [≈4.5291e-7](https://arxiv.org/abs/2006.02560v2) | ≈4.1339e-7 | 8.73% | 9.20% | [48D proof](verifications/dim48/) · [Sources & formulas](COMPARISONS.md) |

Densities shown approximately to five significant digits; percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. The 4D and 5D baselines are newer paper-reported values (8 and 7 decimal places), not independently replayed here. 

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
