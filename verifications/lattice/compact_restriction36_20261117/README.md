# 36D lattice upper bound by dimension restriction and free-origin optimization

**Lattice-only density <2.727705e-5.** Full independent audit passed, including all 216,069 source-geometry nodes, 153 exact cap tests, 653 scalar annuli, exact Fourier factors and the entire infinite tail.

The method restricts the certified 37D three-feature compact harmonic kernel to 36D, optimizes its origin coupling subject to the exact Schur PSD condition, and rebuilds the scalar Fourier witness and mod-2 annular closure. This improves the preceding repository certificate 2.7402464246e-5 by 0.457675%; the cited external Zubrilina Table 1 bound 2.86680609e-5 by 4.852128%. No unrestricted-packing or exhaustive worldwide-priority claim is made.

Install `requirements.txt`, then run `python run_all.py`. The package is self-contained; no optimizer or network is required. The runner authenticates and decompresses inputs into a temporary directory, replays the complete source geometry and all scalar checks, and rejects unsuccessful subprocesses. `--quick` checks only exact algebra and cap inequalities; it is not a full proof replay. Do not run Python with `-O`.

[Theorem and audit scope](THEOREM_REVIEW.md) explains dimension transfer, the free-origin argument, lattice-plane coverage, annular counts and analytic dependencies. Source geometry is the published 37D package, copied here with its frozen SHA256. Interval replay uses Arb, not a second independent interval backend. Large inputs are losslessly compressed; manifests authenticate stored and decompressed bytes.
