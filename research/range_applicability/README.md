# Range applicability: verified, noncompetitive diagnostics

These examples test the cyclic method outside dimensions 17–20. They are not new records and do not change the repository results registry.

- **16D:** strict upper bound <0.025150156, worse than the published 0.02492121 in [arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2).
- **22D:** strict upper bound <0.003908297, worse than the pinned published 0.0037259419689206 from [Cohn's table](https://cohn.mit.edu/sphere-packing/). This pinned table is not an exhaustive current literature claim.

Both complete producer proofs and independent exact/interval audits passed. The original positive-definite witnesses are preserved. Exact origin Schur saturation is derived separately and never substituted into a strict-positive-definite source verifier.

## Reproduction

From this directory, install `python -m pip install -r requirements.txt` and run `python verify_all.py`. Do not run Python with `-O`. Verification requires NumPy and python-flint, not an optimization solver. It replays both source continuum proofs and independent proofs, checks frozen certificate hashes, and validates the strict saturated bounds. Full rational subdivision may take several minutes.

See [RANGE_APPLICABILITY.md](RANGE_APPLICABILITY.md) for scope and exact mathematical inputs, [independent16/REVIEW.md](independent16/REVIEW.md) and [independent22/REVIEW.md](independent22/REVIEW.md) for independent proof details. The periodic cyclic argument and mixed-Gaussian lemma are documented in [the registered 17D source proof](../../verifications/dim17/PROOF.md) and [the repaired proof](../../verifications/dim17/refinement_repaired/PROOF.md). This is a cyclic generalization, not a claim that the kernel satisfies the literal pointwise-sign premise of Theorem 1.1.

16D uses 77 rational domains with 5 exact impossible-triangle exclusions; 22D uses 30 domains with 50 such exclusions. The exclusions prove a strict triangle-inequality contradiction by exact rational arithmetic, not sampled geometry. The independent saturation checks include exact rational objective rounding.

## Search code

`search.py`, `search_conditioned.py`, and `search_repair_ready.py` plus `export.py` and `export_guarded.py` preserve the numerical search and rational export machinery. Install `requirements-search.txt` only for searches; invoke each script with `--help` where supported and inspect its command-line parser. Solver output alone is never a certificate. These finite-basis experiments do not establish improvement in 16D, 22D, or every dimension in the implementation range.
