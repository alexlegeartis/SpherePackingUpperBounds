# Seven-dimensional unrestricted sphere-packing density

This frozen, independently audited computer-assisted certificate proves
density strictly below **0.320584772** for arbitrary sphere packings in R^7.
It lowers the paper-reported **0.32111475** in
[arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1)
by approximately 0.165043% and closes approximately 2.052835% of its additive
gap to E7 density pi^3/105. No exhaustive priority or optimality claim is made.

    python -m pip install -r requirements.txt
    python verify.py
    python audit_independent.py .

No solver or network is needed after dependencies are installed. Never use
Python optimization flags. The secondary audit requires the package directory
argument. Programs can be launched from another working directory.

Full producer replay passes exact PSD, both global polynomial signs, Schur
saturation, 212 pair boxes, 6,487 far boxes and 81 compact cyclic boxes.
The independent audit also passes a distinct pair checker with all five
threshold endpoint allowances. Gaussian repair weights are both zero.

The outward density enclosure is centered at

    0.3205847715830073938695391744713529466970446527646514319059435920441989479863

with radius below 6.1e-77. The strict decimal is rounded upward. The volume
factor is pi^3/840 and the Laguerre parameter is 5/2. See PROOF.md and
INDEPENDENT_AUDIT.md for scope and dependencies.

The five published spherical SDP neighbor-count inputs are external
mathematical dependencies. This is not a proof-assistant formalization or a
clean-room implementation of every continuum checker.

The repository wrapper checks frozen source hashes and the exact registered threshold, then executes the standalone primary proof and separate independent full audit. The original primary entrypoint is preserved as `verify_primary.py`.
