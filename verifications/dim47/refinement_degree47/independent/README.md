# Independent 47D degree-47 radial-enrichment audit

Run `python verify_all.py` with Python 3.11+ and python-flint 0.9.0. Do not use
Python -O. No optimizer, NumPy, SciPy, network, or files outside this package
are needed. Read `AUDIT.md` for the full upper-bound proof review and the
saturation-floor comparison. `SATURATION_EXTENSION.md` proves the scalar
free-amplitude extension used in that comparison.

All source inputs and embedded lower proofs are pinned by `INPUT_SHA256.json`.
The aggregate result is `audit_result.json`; all earlier packages remain
unchanged.
