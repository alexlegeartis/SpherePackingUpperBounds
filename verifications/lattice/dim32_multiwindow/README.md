# Stronger 32D lattice bound: two source windows

**Density ≤0.000102834508.** This supersedes the earlier 0.000104134823 bound in [the preserved original32D package](../dim32/README.md).

Run `python run_all.py` with Python3.10+ and SymPy. Public `python-flint` is optional acceleration; the same exact checks work without it. Do not use `python -O`. The default complete independent verifier checks all five angular cuts, each cut's declared source-prefix union and endpoint, 240 radial intervals, all193 class branches, the Fourier global sign by exact root count, the tail by exact rational Bernstein subdivision followed by positive shifted coefficients, and the final exact rational density. No sign condition is skipped.

There are two nested source windows, not a disjoint partition. Each global population inequality holds separately, so multiplying by nonnegative weights and adding remains valid despite overlap. The verifier explicitly uses the entire prefix of bins indicated by each cut, rather than only the first bin. No prior-density or arithmetic trace-gap assumption is used. [METHOD.md](METHOD.md) and [ANGULAR_PROOF.md](ANGULAR_PROOF.md) explain the continuous geometric proof; the window-specific indicator is in METHOD.md.

SHA256: `7e79ec20cbcdda94e569193379b5197f230a24da0fff7e34c6b748ff3b1b7c87`.

Optional independent producer replay: `python verify_sympy.py certificate.json` (exact tail root counts can be slower). Audit reports are provenance, not executable-proof substitutes. The compressed replay materializes and hash-checks a temporary directory. The original32D certificate remains unmodified.

The original comparison certificate n032.json is preserved from Nina Zubrilina's https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb , accompanying *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . It is a benchmark, not a proof assumption. The Cohn-snapshot comparison is documented in the top-level lattice table; no optimality or exhaustive-priority claim is made.
