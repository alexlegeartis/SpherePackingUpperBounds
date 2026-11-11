# Angular discovery provenance for dimension 47

The frozen source_candidate.json is the fifth angular LP checkpoint after exact angular-critical-point separation. It is numerical discovery data, not a proof.

Run angular_refit.py --dimension 47 --degree 19 --alpha 40 --rounds 5 to reproduce the search design. This needs NumPy, SciPy and python-flint. The included dimension-48 reference determines column conditioning only. A new optimizer run need not reproduce identical floating-point coefficients.

export_angular.py rationalizes the saved checkpoint with denominator bound 10^12, subtracts exactly 1/1000 from its constant harmonic coefficient, then renormalizes that coefficient to 1. It writes the exact angular witness. This export has been replayed; no new optimizer replay is implied.

The parent verifier uses only the frozen exact witness. This folder records the angular stage. The subsequent radial reoptimization and its exact export are documented in radial/. The full package uses those new radial coefficients, not the earlier unchanged-radial substitution.
