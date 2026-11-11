# 47D sphere-packing density upper bound

The exact dimension-transfer and angular-origin saturation construction proves

**Delta_47 < 0.000000616879891.**

It improves the pinned published comparator 0.0000006490757338 by about 4.96%. It is not an optimality or worldwide-priority claim. Read [PROOF47.md](PROOF47.md) for the complete transfer proof, assumptions, normalization, provenance and comparison caveat.

## Reproduce

Requires Python 3.11+ and python-flint 0.9.0. From this directory:

    python -m pip install -r requirements.txt
    python verify.py

Do not use Python -O. No optimizer or network is needed for verification. The complete replay checks frozen source blobs, all source48 sign/coupling/cyclic domains, exact dimension47 Fourier identities and harmonic coefficients, the exact rank-one origin Gram factorization, and rigorous decimal rounding. The inherited source48 scripts write their audit outputs into their local copied directory.

`verify_transfer.py` alone checks only the new exact algebra; it does not replace the full inherited proof replay. `certificate47_transfer.json` and `transfer_exact_audit.json` store derived rational data and intervals. Source48 data and proof dependencies are preserved under `continuum_three_point/`.

The `*_screen.py` scripts are numerical research tools, not proof gates. Their parameter scans found no competitive direct transfer for dimensions 25–46. They use NumPy, SciPy and mpmath. No numerical optimization is part of the registered 47D proof.

The registry wrapper also runs the independent full inherited-source and transfer audit in `independent/`, validates frozen certificate hashes and checks the registered threshold. See [independent review](independent/REVIEW.md). `verify_all.py` is the producer-only full chain.
