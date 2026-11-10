# Annular three-point discovery across dimensions

## Update

The fine degree-31 follow-up now gives a locally verified 18D bound <0.013408931, improving the pinned comparator by about 0.5666%. See [PROOF_DIM18.md](PROOF_DIM18.md). Independent review passed; the certificate is registered at ../../verifications/dim18/. Literature freshness caveats remain. Earlier coarse runs produced noncompetitive bounds. See [GLOBAL_POSITIVITY_FOLLOWUP.md](GLOBAL_POSITIVITY_FOLLOWUP.md). The section below documents the original sampled-only experiment.

## Original result

**No new certified density upper bound was obtained in this run.** This package is reproducible discovery and verification infrastructure, the original experiment did not add an entry to RESULTS.json. It tests genuine three-point constraints, unlike the earlier two-point Gaussian experiments.

The source framework is Cohn–de Laat–Salmon, *Three-point bounds for sphere packing*, Theorem 1.1, https://arxiv.org/abs/2206.15373. The exact verification routines extend the dimension-20 implementation in alexlegeartis/SpherePackingUpperBounds at commit 85e395353a937369e5cff6400f66337967812270. Reference files are included unchanged for regression testing. The kernel is a sum of radial-hat PSD blocks times normalized Gegenbauer kernels, with an isolated origin coordinate. It is a general sphere-packing density framework, not a kissing-number bound.

## Experiments and failures

Targeted dimensions were 17, 18, 22, 23, which had no registered certificate in the repository. All runs used physical/Fourier degree 39, five radial nodes, angular degree 16, exclusion squared n/8, and truncation squared 1.2 times exclusion squared.

The initial sampled SDP gave numerical candidate densities 0.0179942506553, 0.0129719244189, 0.00322108794540, and 0.00211673735276 respectively. **Every candidate failed exact PSD verification. These numbers are not upper bounds.** Clarabel returned optimal_inaccurate and maximum constraint violations between 6.4e-7 and 4.7e-6.

Larger strict feasibility margins produced candidates 0.0205364607919, 0.0159257576320, 0.00920420547679, 0.00963733647824. Their exact matrices were positive definite, but **every candidate failed exact Fourier positivity**. These are also not bounds. The outputs and failure logs are retained to prevent accidental promotion of solver estimates.

The generic exact verifier successfully replays the existing 20D certificate, including the Fourier identity, global Fourier and exterior signs, all 561 angular corner polynomials, and the complete annular inequality. The verified existing density is <0.006979502. This regression is not a new result.

## Files

- `search_initial.py`: initial sampled SDP discovery.
- `search.py`: strict-margin sampled SDP discovery; arguments allow dimensions 3–48.
- `rationalize.py`: exact rational candidate export, expressly unverified.
- `verify_global.py`: exact LDL, Fourier identity, Sturm signs, full angular rectangles, and interval objective.
- `verify_annular.py`: full annular inequality using Arb Taylor bounds.
- `reference/`: original 20D witness and verification sources, with proof attribution.
- `candidates*.json`, `witnesses*/`, `audit*.log`: preserved experiments. Only the audited degree-31 18D witness is registered; successful17D and coarse18D witnesses are noncompetitive; all other claimed failures remain rejected.

Dimensions 1 and 2 already have known sharp density bounds. No strict improvement is possible there; similarly dimensions 3, 8, and 24 are solved. Search input is restricted to dimensions 3–48 because the displayed Gegenbauer recurrence assumes n>2. This package does not invent records in solved dimensions.

## Reproduce

Install `requirements.txt` in a virtual environment. From this directory:

    python search_initial.py --dimensions 17 18 22 23 --output candidates.json
    python search.py --dimensions 17 18 22 23 --output candidates_margin.json
    python rationalize.py candidates_margin.json witnesses_margin
    python verify_global.py witnesses_margin/dim17.json

The last command is expected to fail at Fourier positivity. To replay the valid reference:

    python verify_global.py reference/certificate.json
    python verify_annular.py reference/certificate.json

Both verifiers must succeed on the same exact witness before any certificate is accepted. Never run with Python optimization (-O). The verifier currently requires strict positivity of PSD pivots and sign polynomials, so failure does not establish mathematical infeasibility. Numerical solver outcomes may differ across platforms; retained exact witnesses make the reported rejection replayable.

## Research next step

Replace sampled Fourier/exterior constraints with globally positive SOS or Bernstein certificates, or add exact root-isolation cutting planes until the resulting rational polynomials pass. Strengthen angular coverage using full fine-cell corner intervals rather than just sampled admissible angles. The present exact rejection mechanism is intentionally fail-closed; a dense sample is never proof of a global sign. The broader 1–48 benchmark table and newer 2025/2026 bounds should be managed separately, and only a fully certified improvement over the current comparator should enter the result registry.
