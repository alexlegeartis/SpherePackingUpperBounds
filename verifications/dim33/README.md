# Dimension 33: Gaussian cyclic density certificate

Strict unrestricted sphere-packing density upper bound **<0.00009009682341**.

The central idea is a dimension-specific angular PSD kernel with Gaussian damping, paired with a freshly optimized radial Fourier pair. A rank-one origin Gram reduces the objective; grouping triangle contributions into cyclic sums relaxes individual sign constraints. See [PROOF.md](PROOF.md) and [independent/AUDIT.md](independent/AUDIT.md).

Run `python verify.py` within the repository for registration, hash binding and both full primary/independent chains. Standalone copies run `python verify_primary.py` and `python independent/verify_all.py`. Python3.11+ and python-flint0.9.0 are required, with assertions enabled. No optimizer, network or external near-neighbor count bound is needed.

The independent audit reconstructs the Fourier identity and harmonic PSD, checks all cyclic domains and analytic tails, proves global radial signs by exact rational Sturm and verifies the full coupled interval and rational outward rounding. No unverified numerical samples stand in for proof domains.

[Discovery and exact export](../../research/dim33_gaussian_discovery/) reconstruct the frozen witnesses. No fresh optimizer rerun is claimed. Comparators are pinned published inputs, not exhaustive worldwide-priority assertions; lower constructions are not independently rebuilt.
