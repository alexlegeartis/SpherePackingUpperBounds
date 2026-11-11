# Dimension 44 Gaussian cyclic density certificate

The full primary verifier proves Delta_44 < 0.000001842868599, versus the pinned Cohn-table comparator 0.0000019017703144 (about 3.097% lower). Current worldwide priority is not exhaustively established.

Status: primary and independent exact/interval audits PASS.

Run with Python 3.11+ and python-flint 0.9.0:

    python -m pip install -r requirements.txt
    python verify.py

No optimizer, network access, or unpublished input is required. The replay checks fresh 44D angular PSD and signs, the complete cyclic domain, all 112 Fourier coefficients, global Fourier/physical signs, radial coupling, origin Gram normalization, and exact rational upward rounding. Assertions must remain enabled.

Read PROOF44.md for the theorem and normalization. certificate.json is the exact radial witness. continuum_three_point/cyclic_angular40_exact.json is the exact new 44D angular witness. The copied historical 48D notes are background only; their old numerical parameters and bounds are not substituted into this proof.

This is a computer-assisted upper bound, not a packing construction, optimality proof, peer review, or proof-assistant formalization.

The registry wrapper runs `verify_primary.py` and the separate independent audit, with frozen hashes and the registered threshold checked. See [the independent audit](independent/AUDIT.md) for the distinction between full source replay and separately implemented algebra/pair checks.
