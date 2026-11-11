# Dimension 26: Gaussian cyclic certificate

Strict unrestricted density upper bound **<0.0009746292767**. The central idea is a dimension-specific positive-semidefinite angular kernel, Gaussian radial damping and a coupled Fourier-positive radial function. Triangle terms need only have a nonpositive cyclic sum; an optimized rank-one origin Gram keeps the objective small. See [PROOF.md](PROOF.md).

Run `python verify.py` inside this repository for registration, input hashes and both complete proof chains. For a standalone primary proof run `python verify_primary.py`; the independent audit is in independent/verify_all.py. Assertions must remain enabled. Install requirements.txt; no optimizer or network is needed for proof replay.

The cited upper is a pinned comparison, not an exhaustive worldwide-priority assertion. No external near-neighbor count bound is needed. Discovery code is separate from proof inputs. Earlier research and certificates remain preserved.
