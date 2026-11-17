# Verified 82D mixed-depth binary-coset gluing

Ordinary density ≥8.56691236634380e-18; minimum 8; determinant 55/16; exact center-density square 2^86/55. This improves the pinned preceding repository result by 14.0072495691%. The derived MW80×A2 example is an attributable external comparator, not asserted strongest; improvement over that example is 237.0392332631%.

Run `python run_all.py` for all exact input/parent/local-coverage/automorphism/Gram/density checks and archived-search integrity. `python run_all.py --full` additionally compiles and replays every primary and independent exclusion for 101,108 and 203, with inherited 81 controls. Requirements: Python 3, numpy, python-flint and g++; assertions enabled. Replays run on temporary copies.

Primary historical PENDING_AUDIT labels are retained intact. Complete independent audits now PASS, including reverse-order searches of 108 (5,099,300,880 nodes / 656 leaves) and 203 (2,489,483,280 nodes / 0 leaves). The primary 108 traversal has 100 leaves because its variable order differs; equal node counts do not imply equal leaf counts.

See l80_code_decoder_20261117/CASE82_PROOF.md and l80_case108203_independent_20261117/AUDIT.md. Parent minimum 8 is inherited from Bachoc-Nebe Theorem 4.7 with its constructive hypotheses verified. No claimed catalogue identification or optional local-extremality assertion is needed.
