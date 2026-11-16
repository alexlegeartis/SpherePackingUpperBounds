# 23D lattice-only upper bound

Certified density **≤0.002608401581**. Method: **Conditioned nested Gram-moment class profiles**. See [the proof](METHOD.md).

Run `python run_all.py` with Python 3.10+ and SymPy; optionally install public `python-flint` for faster global signs. Assertions must remain enabled. The runner pins the certificate SHA256, dimension and bound, and repeats all 480 radial envelopes, 385 profile duals, conditioned and unconditional rank inequalities, global signs and odd-dimensional normalization. No prior density bound or arithmetic quotient theorem is assumed.

From a repository clone, use `python verifications/lattice/replay.py --dimension 23`; it checks and decompresses the stored files in a temporary directory. JSON files are compressed deterministically; MANIFEST.json records original and stored hashes.

Comparison only: Cohn's lattice table gives prior upper 0.0026135387306483, inherited from [Kravatskiy's unrestricted signed logarithmic profile + degree95 radial enrichment](../../dim23/refinement_combined20261115/). That package states the upward-rounded strict bound 0.0026135387307. The new improvement is lattice-only; the unrestricted result is unchanged. Known lattice lower comparison: 0.0019053281934260 (truncated Cohn display), not the larger nonlattice lower bound. Source: https://cohn.mit.edu/sphere-packing/ , fetched 2026-11-16. No exhaustive priority claim.
