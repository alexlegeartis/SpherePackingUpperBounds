# 14D curved-profile cyclic certificate

**Complete primary replay PASS: density expression < 0.045003383.**
The expression is 0.045003382347615314738618663339172944… . Independent
review is still required before treating this as an independently audited
record. The supplied newer comparator is 0.04501031 from Table 1 of
https://arxiv.org/html/2512.04829v2.

Run `python verify.py` with Python 3.11+ and python-flint 0.9.0, without
Python -O or -OO. The verifier needs no optimizer, network or external source
folder. Input hashes are checked before all mathematical gates. Both proof
entrypoints reject disabled assertions; four negative tests are recorded.

The profile is h(u)=exp[-14(u/s-1)-8(u/s-1)^2], with squared separation s=7/4.
The angular polynomial was reoptimized and reconstructed in dimension 14;
no 15D Gegenbauer coefficients were reused as a certificate. The kernel uses
nonnegative scalar angular harmonics and an exact rank-one origin extension.
Its actual angular polynomial degree is 12, within the degree-19 search span.

Exact rational cap/sign-transition and full cyclic geometry pass: 5,537
accepted rectangles, 371 sorting exclusions, maximum axis depth 15, analytic
unbounded tails and the all-dilation monotone-ratio lemma. The radial function
has degree 63. Two exact Fourier reconstructions agree; complete root isolation
proves the Fourier ray and physical negativity from X=12. The full pair
interval [pi*7/4,12] passes with 65 leaves, depth 13. The remaining ray follows
from physical negativity and 0<h<=1.

The exact positive Fourier normalization is retained, rather than rounded to
one. The radial multiplier guard is 1000001/1000000. The final density formula
is ((pi*s/4)^7/7!)*(P(0)+MB)/H(0); a rational Machin upper bound for pi proves
the strict target without relying on floating rounding. See THEOREM.md.

The angular guard is 1/10000, followed by renormalization to constant harmonic
coefficient one. `source_candidate.json` is the frozen numerical angular lead;
`angular_exact.json` is deterministically reconstructed and checked, without
rewriting it. `certificate.json` is the final exact radial data.

The discovery directory preserves all three two-round profile screens,
all rejected root-guided checkpoints, the narrow pair counterexamples, and
the final coefficients and cuts. A fixed radial cutoff X=12 eliminated the
need to rely on tiny far-pair cancellations; it was verified afresh in 14D.
The last pair lobe near X=9.474376 was repaired by adding its numerical maximum
and nearby points; only the subsequent complete interval proof establishes
feasibility. Nine radial repair solves were used after the initial cut-12
candidate. No rejected snapshot is a density bound.

The final certificate and mathematical inputs are frozen for independent
review. `INPUT_SHA256.json` pins them; `SHA256SUMS` covers this package.
`full_replay_result.json` and `full_replay.log` contain the complete primary
result. Earlier 15D certificates and 14D Chebyshev experiments were unchanged.
