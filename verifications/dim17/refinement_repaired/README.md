# Dimension 17: mixed-Gaussian repair and origin saturation

Strict unrestricted density upper bound: **<0.018204096**.

The preserved positive-definite source certificate proves <0.018250915 using two Gaussian scales. Its raw polynomial-Gaussian tail alone is invalid: verify_repaired.py and slow_gaussian.py are required. The independent repair_checks.py then proves exact origin Schur saturation, yielding the registered stronger bound. Do not replace the source origin entry with a singular saturated matrix and run the strictly-positive-pivot source verifier.

The source certification-status label records its pre-audit history. Its mathematical inputs and bytes are frozen; completed independent status is documented in REVIEW.md.

    python verifications/dim17/refinement_repaired/verify.py

Run from repository root after installing requirements.txt, without Python -O. Both producer and independent domain paths must pass. The independent proof covers all 56 compact and 21 far cells; exactly three impossible-triangle boxes are excluded using strict rational inequalities, never floating-point geometry. Exact truncation square is 29749999999999996/10^16, not 119/40.

The cyclic generalization and periodic summation remain essential; this is not the pointwise-sign hypothesis of the literal Theorem 1.1. See PROOF.md and REVIEW.md for corrected Fourier normalization, full mixed-tail repair and exact PSD saturation.
