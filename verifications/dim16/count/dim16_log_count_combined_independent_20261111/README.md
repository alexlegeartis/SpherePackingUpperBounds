# Reproduce the complete independent 16D audit

Run `python run_all.py` from any working directory. The script rejects optimized Python and uses a temporary copy for all external count replays, leaving their frozen inputs untouched.

The following sibling directories must be preserved:

- dim16_log_hinge_geometry_independent_20261111
- dim16_log_count_radial_independent_20261111
- density_kissing16_reweight_frozen_20261111
- density_kissing16_independent_pair_20261111
- density_kissing16_independent_cyclic_20261111

The three count directories are the same frozen count components supplied with the primary combined certificate. Mathematical hashes are pinned in run_all.py. Its temporary count replay independently reconstructs factors and polynomials before checking every exact covering event; it does not merely trust a cached success JSON.

A successful run rewrites COMBINED_RESULT.json and logs. AUDIT.md explains the theorem linkage, the scope of each independent calculation, the comparison, and the limitations. Ten negative controls are included. The full aggregate run takes approximately one minute in the tested single-thread environment.

SHA256SUMS and manifest.json cover immutable code, inputs and documentation. Timing-dependent generated logs and results are deliberately excluded, so replay does not invalidate source provenance.

The verified strict density bound is 0.024920988. This is a modest improvement over the published displayed 0.02492121, not an optimality claim.
