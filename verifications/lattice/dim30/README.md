# 30D lattice-only density bound

**Density ≤ 0.000207932946.** This is an unconditional lattice bound, not an unrestricted sphere-packing result.

Run `python run_all.py` with Python 3.10+ standard library. Do not use `python -O`. The complete independent verifier reconstructs the Fourier polynomial, all120 radial intervals, continuous angular inequalities, projective caps through sextic degree, every class budget, and the exact density normalization. Global Fourier and tail signs are proved by exact compactified-half-line Bernstein subdivision. Odd-dimensional normalization is an exact squared-density comparison using a rational lower bound for pi. No floating tolerance or prior density assumption enters the proof.

Optional second implementation: `python verify_sympy.py certificate.json` requires SymPy; optional public `python-flint` speeds it up, although exact tail root counts may remain slow. It is not required for the complete default proof. The verifier writes a local JSON audit report; the compressed release runner uses a temporary directory and checks stored and decompressed hashes.

Certificate SHA256: `dd5bb853132587d1a7169ab7bc3afbfdd077e4966336d1073e31d3de81b129db`.

## Proof and attribution

See [METHOD.md](METHOD.md), [ANGULAR_PROOF.md](ANGULAR_PROOF.md), and [PROJECTIVE_CAP_PROOF.md](PROJECTIVE_CAP_PROOF.md). These combine Zubrilina's modulo-two class framework with global short-class angular coupling. The full standard-library proof reconstructs each cap's Gegenbauer expansion. `gegenbauer.py` preserves the independent auditor's exact recurrence. Established Delsarte/Gegenbauer and projective polynomial techniques are not claimed novel.

`comparison_certificate.json` is unchanged upstream `certificates/n030.json`, from https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb , attributed to Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . It is a named historical benchmark, not a proof dependency. Optional source replay: `python ../prior/taylor.py comparison_certificate.json`.

`independent_audit.json` preserves the pre-packaging full audit for the immutable candidate. Audit reports are provenance, not substitutes for executable proof. No optimality or exhaustive worldwide-priority claim is made.
