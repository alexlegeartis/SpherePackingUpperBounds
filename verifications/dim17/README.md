# Dimension 17: solid-harmonic cyclic bound

Strict unrestricted sphere-packing density upper bound: **<0.018335552**.

This witness uses a cyclic generalization, not the literal individual-kernel-sign assumption of Theorem 1.1. PROOF.md and REVIEW.md supply the finite PSD triangle-summation argument, periodic Poisson normalization and periodic approximation. No external neighbor-count input or lattice-only restriction is needed.

From the repository root, install requirements.txt and run:

    python verifications/dim17/verify.py

The wrapper checks the frozen witness and registry, then runs the producer and independent full-domain verifiers. Do not use Python -O. Exact exclusion squared is 17/8 and truncation squared 51/20.

The independent computation reconstructs the Fourier transform and solid Gegenbauer polynomials separately and uses exact rational Bernstein subdivision for all 20 compact and 10 far cells. It proves the pair inequality with 384-bit Arb and the outward decimal with a rational Machin bound. Source hashes, metadata reconciliation and audit limitations are documented in REVIEW.md.
