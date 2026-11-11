# Dimension 20: cyclic refinement

Strict unrestricted density upper bound: **<0.006977702**. Earlier proof and exact witness remain in the parent directory; the earlier origin-slack refinement remains in ../../origin_slack/.

This uses a cyclic generalization, not the literal pointwise-sign assumption of Theorem 1.1. PROOF.md and REVIEW.md supply the finite PSD triangle-summation proof, periodic Poisson normalization and unrestricted periodic approximation. No neighbor-count input is used.

After installing root requirements.txt:

    python verifications/dim20/refinement_cyclic/verify.py

Both producer and independent paths must pass. Independent domain checks use exact rational Bernstein subdivision for all 20 compact and 10 far cells, without geometric exclusions. The exact radii and metadata reconciliation are documented in REVIEW.md. Do not use Python -O.
