# Sphere-packing density upper bounds

Reproducible certificates for unrestricted sphere packings, organized by dimension on the model of KissingNumbers. These are upper bounds on density, not packing constructions, kissing-number results, or optimality claims. No worldwide-priority claim is made.

## Registered results

| Dimension | Known lower density L | Older upper U | New strict upper N | Upper decrease | Additive gap decrease |
| --- | --- | --- | --- | --- | --- |
| 20 | 0.0040241650630681997943 | 0.0071270536033763 | 0.006979502 | 2.070303% | 4.755298% |
| 48 | 2.3178295310541237588e-8 | 4.529067791E-7 | 4.13385739314e-7 | 8.726087% | 9.196747% |

[References, exact lower-bound formulas and comparison conventions](COMPARISONS.md). [20D proof](verifications/dim20/) · [48D proof](verifications/dim48/).

Both use the existing [Cohn–de Laat–Salmon continuum framework](https://arxiv.org/abs/2206.15373). Proof notes state the analytic hypotheses and normalization; exact rational witnesses and independent interval verifiers check the computational conditions over full domains. Historical numerical comparisons in package notes are pinned comparisons, not current best-known claims.

## Reproduce

Python 3.11 or newer is required. Install dependencies and run:

    python -m pip install -r requirements.txt
    python run_all.py

Never disable assertions (`python -O`). No optimization solver or network is needed after dependency installation. The 48D compact triangle replay may take several minutes. Logs record the original audits; fresh replay is recommended. No GitHub Actions workflow is installed.

## Layout

- `verifications/dim20/`: exact witness, proof and two verification paths.
- `verifications/dim48/`: current cyclic witness, complete proof dependencies and independent angular, triangle-domain and radial audits.
- `RESULTS.json`: machine-readable registry.
- `run_all.py`: full certificate replay, failing on any unsuccessful subprocess.

The referenced continuum theorem and standard positive-definite Gegenbauer addition theorem are mathematical inputs, not independently formalized theorems in this repository.
