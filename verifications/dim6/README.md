# Six-dimensional unrestricted sphere-packing density

This frozen computer-assisted certificate proves density strictly below **0.410189987** for arbitrary sphere packings in R^6, conditional on the standard Fourier/PSD-kernel facts and the four published neighbor-count bounds stated in the proof.

The full producer replay and independent source-reviewed audit pass. This is not an optimality theorem, a claim of external peer review, a formal proof-assistant verification, or an exhaustive worldwide-priority claim.

The comparator is the paper-reported **0.4102905** in [arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1). Its original witnesses have not been independently reproduced here. The new strict upper decimal is smaller by 0.000100513, approximately 0.02450% of that upper bound.

## Replay

Python 3.11 or newer is recommended. Install the verification dependencies and run:

    python -m pip install -r requirements.txt
    python verify.py

No optimization solver or network is needed after installing dependencies. Never disable assertions; every executable verifier explicitly refuses `python -O`. The tested versions are NumPy 2.3.5 and python-flint 0.9.0.

The replay works from another working directory when given an absolute path to `verify.py`. It reconstructs all rational matrices and Fourier polynomials, validates the Gaussian repair, proves the exact Schur factorization, checks complete pair/far/cyclic domains, and verifies outward decimal rounding. All serialized success flags are recomputed.

For the separately implemented polynomial/Fourier/origin and antipodal arithmetic check:

    python audit_independent.py .

That secondary command alone is not the complete density verifier. See [PROOF.md](PROOF.md) and [INDEPENDENT_AUDIT.md](INDEPENDENT_AUDIT.md).

## Frozen result

The complete outward density enclosure is centered at

    0.4101899868311988677407604362994285017515171298713618040883764151908393153314

with radius below 4.5e-77. The strict published decimal 0.410189987 is rounded upward.

Checks cover 197 pair boxes, 3,097 far boxes and 153 compact cyclic boxes. The original matrices are strictly positive definite; the final saturated origin block is positive semidefinite by an exact rational factorization. Published average-neighbor inputs are mathematical dependencies, not independently recomputed spherical SDP certificates.

`discovery_snapshot.json` preserves the numerical search checkpoint for provenance only. The proof is defined by the exact finite-decimal coefficients, rational matrices and Gaussian parameters in `witness.json`, `witness_psd_exact.json` and `prerequisites.json`, checked by `verify.py`.

The repository `verify.py` wrapper checks the registry entry and frozen mathematical input hashes, then runs the primary complete proof and the separate independent audit. The original primary entry point is preserved as `verify_primary.py`. See `independent/` for the separate audit and the root comparison document for exact lower-density formulas and percentage conventions.
