# Compact harmonic-matrix lattice bounds in46D and47D

Certified lattice-only upper bounds: **46D<7.790880151e-7**, **47D<5.186745468e-7**. These do not bound unrestricted packings.

Method: exact target-dimensional PSD harmonic restriction of the compact two-feature48D kernel, with reoptimized origin coupling, positive-factor Fourier witnesses and mod-2 annular profiles. The changed constant harmonic matrix is used in the origin Schur complement. [Independent audit](AUDIT.md).

Run `python run_all.py` with `python-flint` and `sympy` installed, assertions enabled. The runner checks hashes, invokes the sibling `compact_matrix48_20261117` complete geometry proof, and replays both target scalar certificates in temporary copies. `--quick` runs exact algebra using archived interval receipts; it does not freshly certify interval signs. The entire repository (or both sibling packages) is required. Portable output-directory handling is the only verifier change; mathematical inputs are unmodified.

External comparisons are Zubrilina(2026), [Table1](https://arxiv.org/html/2610.35830v1#S4.T1):46D8.80676408e-7 and47D6.07126117e-7. Prior repository bounds are internal history, not the principal comparison. No global-priority or optimality assertion is made.
