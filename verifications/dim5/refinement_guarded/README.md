# Improved five-dimensional unrestricted density

This frozen, independently audited computer-assisted certificate proves
density strictly below **0.510129239** for arbitrary sphere packings in R^5.
It improves the prior repository result 0.510501094478. Against the pinned
paper-reported **0.5122645** in
[arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), it
lowers the upper by approximately 0.416828% and closes approximately
4.542443% of the additive gap to D5 density pi^2/(15 sqrt(2)).
No exhaustive worldwide-priority or optimality claim is made.

    python -m pip install -r requirements.txt
    python verify.py
    python audit_independent.py .

No solver or network is needed after installing dependencies. Never use
Python optimization flags. The secondary audit requires this directory as
its argument. The full producer and independent audit both pass.

A tiny exact PSD kernel correction repairs the unguarded checkpoint's far
sign lobe without changing the pair inequality. The final source is checked
directly on every complete domain. Gaussian repair weights are both zero.
The unguarded checkpoint must not be described as a bound.

The density enclosure is centered at

    0.5101292386566609905645356651355475179736476519821396988966199837323083554957

with radius below 8.7e-77. The strict decimal is rounded upward. The volume
factor is pi^2/60 and the Laguerre parameter is 3/2.

Standard Fourier/PSD facts and the cited published spherical SDP count inputs
remain dependencies. This is not formal proof-assistant verification or a
clean-room reimplementation of every continuum checker. See PROOF.md and
INDEPENDENT_AUDIT.md.

The repository wrapper checks frozen source hashes and the exact registered threshold, then executes the standalone primary proof and separate independent full audit. The original primary entrypoint is preserved as `verify_primary.py`.
