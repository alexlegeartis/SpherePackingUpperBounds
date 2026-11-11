# Dimension 46 density upper bound

Frozen exact witness: **Delta_46 < 0.000000914882847**.

Read [PROOF46.md](PROOF46.md) for the full construction, global hypotheses, normalization and comparison caveats. The bound is about 1.57% below the pinned published comparator 0.0000009295151556. It is not an optimality or worldwide-priority claim.

## Reproduce

Python 3.11+ and python-flint 0.9.0:

    python -m pip install -r requirements.txt
    python verify.py

Do not use Python -O. No optimizer or network is needed. This checks the frozen source angular/cyclic proof, the new full dimension-46 Fourier identity, both infinite-domain signs, the complete coupled interval, dimension-specific angular PSD, exact origin Gram factorization and rigorous final decimal rounding.

`certificate.json` defines the exact witness. `radial_perturb46_degree15.json` and `audit_perturb.py` preserve its rational export provenance; they are not needed for the standalone frozen verifier. Full inherited source proofs and dependencies are retained under `continuum_three_point/`.

An independently implemented complete audit is included under `independent/`. It uses exact rational Sturm sequences for both global signs, a distinct full-function second-derivative coupling enclosure (10 leaves, depth 8), explicit-formula angular reconstruction, and exact rational final rounding. Its full source replays also passed. Run `python independent/verify_all.py` and read `independent/REVIEW.md`.

The repository `verify.py` wrapper checks the registry entry and frozen mathematical input hashes, then runs the primary complete proof and the separate independent audit. The original primary entry point is preserved as `verify_primary.py`. See `independent/` for the separate audit and the root comparison document for exact lower-density formulas and percentage conventions.
