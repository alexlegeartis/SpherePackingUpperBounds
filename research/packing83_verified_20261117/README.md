# Verified 83D mixed-depth binary-coset gluing

Ordinary density ≥ 2.56962386473426e-18; minimum 8; determinant 23/4; exact center-density square 2^85/23. Improvement over the pinned preceding repository result is 16.5631010300%. The current derived external comparator is the BCS88 section perpendicular to √2 D5, giving improvement **>371.814259%**. [Exact audit and machine-readable comparison](../bcs83_comparator_20261117/). Archived audit files retain their original MW80 × A3 comparison (501.820882483%) for reproducibility; it is no longer the leading table comparison. Neither comparator is asserted strongest worldwide.

Run `python run_all.py` for exact input, parent, local coverage, automorphism, parity geometry, Gram and density checks plus archived search integrity. `python run_all.py --full` additionally compiles and replays all primary and independent exclusions for cases 101, 108, 203, 301 and 302, including inherited controls. Requirements: Python 3, numpy, python-flint and g++; assertions enabled. Replays use temporary copies.

Historical primary PENDING_AUDIT labels remain intact. All independent audits now PASS. The new independent traversals used 7,416,859 nodes for 301 and 150,439,952 nodes for 302, zero leaves, all 16 glues. The primary 302 traversal used 149,419,024 nodes; variable order differs.

See CASE83_PROOF.md in the primary decoder folder and both new independent audit directories. Parent minimum 8 inherits Bachoc–Nebe Theorem 4.7 with constructive hypotheses checked. No catalogue identification is claimed.
