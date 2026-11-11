# 21D Gaussian cyclic refinement

Strict unrestricted density upper bound **<0.004985977000**. This replaces the registered origin-slack value0.005149549312; all earlier21D proofs remain preserved.

The central idea is a freshly optimized dimension-specific angular PSD kernel with Gaussian radial damping, coupled to a Fourier-positive radial pair. A rank-one origin Gram keeps its objective small; triangle contributions are grouped into cyclic sums instead of requiring each term to be nonpositive. See [PROOF.md](PROOF.md) and [independent/AUDIT.md](independent/AUDIT.md).

Within this repository run `python verify.py` for registration, hash gates and complete primary/independent chains. A standalone copy runs `python verify_primary.py` and `python independent/verify_all.py`. Python3.11+ and python-flint0.9.0 are required; never disable assertions. No optimization solver, network or external near-neighbor bound is required by proof replay.

The new bound reduces the pinned published upper by3.366179% and its additive gap by7.789082%. Relative to the prior repository bound, the reductions are3.176439% and7.369051%. Published lower/upper comparisons are cited inputs, not independently rebuilt packings or exhaustive worldwide-priority assertions.
