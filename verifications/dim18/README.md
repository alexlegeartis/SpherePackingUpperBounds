# Dimension 18

Strict unrestricted sphere-packing density upper bound: **<0.013408931**.

The exact rational certificate passes both producer full-domain verification and an independent verifier using a different Fourier reconstruction, explicit Gegenbauer formula, integer-square-root angular caps, Machin bounds and 384-bit Arb annular checks. See PROOF.md and REVIEW.md.

From the repository root after installing requirements.txt:

    python verifications/dim18/verify.py

Do not use Python -O. The exact truncation square is 2.6999999999999997, not 27/10. This improves the pinned published comparator 0.0134853404450862, without claiming worldwide priority. No external neighbor-count bounds are needed.
