# Dimension 34 Gaussian cyclic certificate

Strict packing-density upper bound: **0.00006394047941**. Full primary and separate independent audits passed. This improves the pinned comparator 6.45221967438e-05; no exhaustive worldwide-priority assertion is made.

Within the repository run `python verify.py` to check registration and both complete chains. A standalone copy runs `python verify_primary.py` and `python independent/verify_all.py`. Python 3.11+ and python-flint 0.9.0 are required; do not use optimization flags. No optimizer, network or external neighbor-count input is needed.

See PROOF.md for the summed cyclic Gaussian mechanism, COMPARATORS.json for the pinned upper and independent/comparison.json for the newer lower formula where applicable. The theorem uses cyclic triangle sums, not the literal individual-sign condition of Theorem 1.1. [Discovery and deterministic export](../../research/dim34_gaussian_discovery/) are packaged separately.
