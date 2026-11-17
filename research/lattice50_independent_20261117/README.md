# Reviewed50D lattice-only certificate

Strict density bound: **1.75102693e-7**. See AUDIT.md for theorem, exact scope and source limitations. No SOTA or external percentage comparison is claimed.

Run `python -B verify.py` in this directory with python-flint installed. It checks the finite632 annuli, global spatial/Fourier identity and exact normalization. The previously published scalar rank-two geometric theorem is an inherited dependency, pinned to SpherePackings commit af863e2014561e060af7ee2b0eecffc468f4eca3, path verifications/lattice/rank_two_20261117/. The copied angular source hash is verified. Its full historical geometry replay is not repeated by this wrapper.

The producer snapshot includes earlier exploratory reports; only final annular certificate, exact repaired polynomial, final Schur guard and the checks selected by run_all.py contribute to this audit. The original annular result is separately preserved; replay writes an independent1024-bit result. See exact_result.json and tail_transform.log for completed checks. Shared FLINT/Arb is explicitly not an independent arithmetic backend.

## Research-only publication

This package is intentionally outside the registered results and comparison tables. The verified lattice-only calculation is published for reproducibility, without an established strongest external50D comparator, percentage gain or global record claim. No result for unrestricted packings is asserted. INPUTS_MANIFEST.json binds immutable proof inputs for repeatable checks; MANIFEST.json preserves the original audit snapshot, including its historical output reports.
