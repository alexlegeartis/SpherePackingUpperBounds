# Dimension 45 sphere-packing density upper bound

Frozen exact witness: **Delta_45 < 0.000001311423575**.

This is approximately 1.40% below the pinned published comparator 0.0000013300905665. It is not an optimality or worldwide-priority claim. Read [PROOF45.md](PROOF45.md) for the complete new angular/cyclic/radial construction, normalization and source caveats.

## Reproduce

Python 3.11+ and python-flint 0.9.0:

    python -m pip install -r requirements.txt
    python verify.py

Do not use Python -O. No optimizer or network is needed. The entrypoint verifies the new 45D angular polynomial, the full cyclic triangle domain and tails, all 112 new Fourier coefficients, both infinite-domain signs, the entire coupled interval, PSD and the origin Gram relation, and exact rational final rounding.

The new angular polynomial is in `continuum_three_point/cyclic_angular40_exact.json`; the new radial pair and M are in `certificate.json`. Historical 48D proof notes in the subdirectory supply only general analytic lemmas and must not be read as substitutes for the new witness verification.

The complete independent audit is included under `independent/`. It reconstructs the new angular proof, checks radial signs with exact Sturm sequences, and uses a distinct second-derivative coupling enclosure (12 leaves, depth 8), explicit harmonic formulas and rational final rounding. Run `python independent/verify_all.py` and read `independent/REVIEW.md`. Both independent audit components passed.

The repository `verify.py` wrapper checks the registry entry and frozen mathematical input hashes, then runs the primary complete proof and the separate independent audit. The original primary entry point is preserved as `verify_primary.py`. See `independent/` for the separate audit and the root comparison document for exact lower-density formulas and percentage conventions.
