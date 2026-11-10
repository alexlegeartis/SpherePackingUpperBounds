# Sphere-packing density upper bounds

Reproducible certificates for unrestricted sphere packings, organized by dimension on the model of KissingNumbers. These are upper bounds on density, not packing constructions, kissing-number results, or optimality claims. No worldwide-priority claim is made.

## Registered results

| Dimension | Known lower density L | Older upper U | New upper N | Upper decrease | Additive gap decrease | Proof and details |
| --- | --- | --- | --- | --- | --- | --- |
| 19 | [≈0.0041474](https://arxiv.org/abs/2505.02394v1) | [≈0.0098180](https://arxiv.org/abs/2006.02560v2) | ≈0.0096184 | 2.03% | 3.52% | [19D proof](verifications/dim19/) · [Sources & formulas](COMPARISONS.md) |
| 20 | [≈0.0040242](https://arxiv.org/abs/2505.02394v1) | [≈0.0071271](https://arxiv.org/abs/2006.02560v2) | ≈0.0069795 | 2.07% | 4.76% | [20D proof](verifications/dim20/) · [Sources & formulas](COMPARISONS.md) |
| 48 | [≈2.3178e-8](https://doi.org/10.4153/CJM-1971-081-3) | [≈4.5291e-7](https://arxiv.org/abs/2006.02560v2) | ≈4.1339e-7 | 8.73% | 9.20% | [48D proof](verifications/dim48/) · [Sources & formulas](COMPARISONS.md) |

Densities shown approximately to five significant digits; percentages to two decimals. Rounded display values are not substitute rigorous bounds. Full values and certificate bounds remain in [RESULTS.json](RESULTS.json); percentages use those unshortened values. 

These certificates use the existing [Cohn–de Laat–Salmon continuum framework](https://arxiv.org/abs/2206.15373). Proof notes state the analytic hypotheses and normalization; exact rational witnesses and independent interval verifiers check the computational conditions over full domains. Historical numerical comparisons in package notes are pinned comparisons, not current best-known claims.

## Reproduce

Python 3.11 or newer is required. Install dependencies and run:

    python -m pip install -r requirements.txt
    python run_all.py

Never disable assertions (`python -O`). No optimization solver or network is needed after dependency installation. The 48D compact triangle replay may take several minutes. Logs record the original audits; fresh replay is recommended. No GitHub Actions workflow is installed.

## Layout

- `verifications/dim19/`: exact witness, proof and two verification paths.
- `verifications/dim20/`: exact witness, proof and two verification paths.
- `verifications/dim48/`: current cyclic witness, complete proof dependencies and independent angular, triangle-domain and radial audits.
- `RESULTS.json`: machine-readable registry.
- `run_all.py`: full certificate replay, failing on any unsuccessful subprocess.

The referenced continuum theorem and standard positive-definite Gegenbauer addition theorem are mathematical inputs, not independently formalized theorems in this repository.
