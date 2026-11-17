# Certified joint 51D power-six upper bound

Unrestricted sphere-packing density < **1.37012082e-7**.

Run `python run_all.py` after installing requirements. The replay checks exact block PSD and origin Schur complement, Fourier positivity and infinite tail, all support-crossing angular controls, all 633 annuli and their spherical-code caps, and the complete 3D cyclic partition with 53,227 leaves. Typical replay takes several minutes.

The verifier uses the standard Gegenbauer addition theorem, spherical-code LP bound, Poisson summation, and periodic approximation as mathematical inputs, not formally verified theorems. See THEOREM_REVIEW.md. Independent geometry uses nine entries and 512-bit Arb arithmetic, differing from the producer's symmetry-grouped six-entry contraction. This package is the exact rounded joint candidate at public research commit e23d47fa756568c9e301313ddf0481181bc0d6e4, not the earlier rejected cap-guarded candidate.

No claim of global priority is made here. It improves the previously published project certificate 1.37781396e-7 by about 0.5584%.
