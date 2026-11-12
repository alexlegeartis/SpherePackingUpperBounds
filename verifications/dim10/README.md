# Reconstructed 10D upper-bound certificate

Strict upper density: **0.143330282**. This is a freshly reconstructed
certificate, stronger than the previously reported but unavailable
0.143376013 package. It is not byte-identical artifact recovery.

Run `python verify_all.py` after installing `python-flint==0.9.0`.
Python 3.11+ is required; do not use `-O` or `-OO`. No optimizer or network is
needed for proof replay. See PROOF.md for the cyclic PSD finite-sum argument,
normalization and exact proof boundaries. The initial primary whole-domain
angular proof used 15,220 accepted rectangles and 341 excluded rectangles;
the exact radial proof used 31 intervals. The aggregate command also runs independent closed-form Fourier, exact Sturm-ray,
first-derivative pair, analytic-gradient geometry and full binary-cover audits, plus
five mathematical mutation tests and five optimized-Python rejection controls.

The pinned external comparison is the printed 0.14340271 in
https://arxiv.org/html/2512.04829v2#S3.T1. This is a 0.05051% decrease against
that displayed decimal, not a worldwide-priority claim. Exact coefficients
and full-precision results are retained in the JSON files.

The [discovery and deterministic export package](../../research/log_hinge10_reconstruction_20261112/)
regenerates every mathematical coefficient from frozen optimizer outputs;
the export was independently compared against this proof package.
