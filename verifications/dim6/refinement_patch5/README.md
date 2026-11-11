# Direct six-dimensional follow-on certificate

This separately frozen computer-assisted certificate proves density strictly
below **0.409128533** for arbitrary sphere packings in R^6. It follows the
audited 99/100 mixture certificate <0.409135399 and is kept separate so it does
not silently change that integration target.

Against paper-reported 0.4102905 in
[arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), it lowers
the upper by approximately 0.283206% and closes approximately 3.111610% of
the additive gap to E6 density pi^3/(48 sqrt(3)). No optimality or exhaustive
worldwide-priority claim is made.

    python -m pip install -r requirements.txt
    python verify.py
    python audit_independent.py .

No optimization solver or network is needed after installing dependencies.
Do not use Python optimization flags. The secondary audit requires this
package's directory argument, even when launched from another working directory.

The complete producer replay and independent full audit pass. Both Gaussian
weights are zero; global signs are certified by exact Sturm arithmetic.
The source is a directly prepared checkpoint with exact origin saturation,
not a mixture. All off-origin diagonal repair shifts are zero.

The producer's density enclosure is centered at

    0.4091285329176186823329127475902823329974333183903733542684807075086836562255

with radius below 5.5e-77. The strict decimal is rounded upward. Standard
Fourier/PSD facts and the published spherical SDP neighbor-count inputs in
PROOF.md remain dependencies. This is not a proof-assistant formalization or
a clean-room reimplementation of every continuum checker.

The repository wrapper checks exact source hashes and the registered threshold, enforces zero Gaussian guards and runs both the standalone primary proof and separate independent audit. The original producer entrypoint is `verify_primary.py`.
