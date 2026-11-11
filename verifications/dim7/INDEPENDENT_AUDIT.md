# Independent seven-dimensional audit

PASS: density strictly below **0.320584772**. Frozen source SHA-256:
1cf76290693c7461ad4da8e5389db283ea5bce02c152ea7ab0115bbf88948bf1.

Separate exact arithmetic reconstructed half-integer Laguerre coefficients and
the full Fourier identity, all 66 positive rational PSD pivots, both global
Sturm signs, the Schur factorization, pi^3/840 normalization and strict outward
rounding. No Gaussian repair or off-origin diagonal shift is used.

A distinct pair checker passed 2,727 leaves and explicitly checked the smaller
allowance at all five half-open count thresholds. Fresh full continuum replay
passed 6,487 far boxes and 81 cyclic boxes. Table 2's dimension-seven counts
136 through 140 and epsilon thresholds were independently checked against
arXiv:1206.2608; each witness radius is reduced by 1e-6 from the publication.
The comparator 0.32111475 was checked against arXiv:2512.04829v2 Table 1.

Standard analytic facts and the original spherical SDP certificates remain
external dependencies. This is not formal proof-assistant verification or a
clean-room reimplementation of every continuum algorithm.
