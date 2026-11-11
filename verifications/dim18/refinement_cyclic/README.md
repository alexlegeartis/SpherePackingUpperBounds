# Dimension 18: cyclic refinement

Strict unrestricted density upper bound: **<0.013284728**. Earlier 18D proof and exact witness remain in the parent directory.

This uses a cyclic generalization, not the literal pointwise-sign assumption of Theorem 1.1. The finite PSD triangle-summation proof, periodic Poisson normalization and unrestricted periodic approximation are supplied in PROOF.md and REVIEW.md. No neighbor-count input is used.

After installing root requirements.txt:

    python verifications/dim18/refinement_cyclic/verify.py

Both producer and independent paths must pass. Independent domain checks use exact rational Bernstein subdivision for all 20 compact and 10 far cells, with no geometric exclusions. The exact exclusion square is 9/4; truncation square is 2.6999999999999997, not 27/10. Do not use Python -O.
