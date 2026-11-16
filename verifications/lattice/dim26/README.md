# 26D lattice-only density bound

**Density ≤ 0.000907209997.** This is an unconditional lattice bound, not an unrestricted sphere-packing bound.

Run `python run_all.py` with Python 3.10+ using only the standard library. Do not use `python -O`. This complete verifier reconstructs the Fourier polynomial, all 120 radial intervals, continuous angular cuts, exact projective code caps through sextic degree, every class budget, and the final rational density. Both global half-line signs are proved by exact Bernstein subdivision after rational compactification. No sign test is skipped and no floating tolerance is used. The separate optional replay `python verify_sympy.py certificate.json` uses Python+SymPy and exact root counts; optional `python-flint` accelerates it, but tail root counts may still take much longer. It is not necessary for the complete standard-library proof.

The verifier writes its JSON result locally. The compressed release runner materializes a temporary copy and checks stored/decompressed hashes. Certificate SHA256: `f83d66a3746811c8887a0095a3df112ff0b322be24981634c02e493ae022778a`.

## Proof and attribution

[METHOD.md](METHOD.md), [ANGULAR_PROOF.md](ANGULAR_PROOF.md), and [PROJECTIVE_CAP_PROOF.md](PROJECTIVE_CAP_PROOF.md) provide the full radial/global angular argument and the continuous multiplicity bounds. There is no arithmetic trace gap, determinant lemma, or prior-density assumption. The exact Gegenbauer recurrence in `gegenbauer.py` is extracted unchanged from the independent auditor's module; `verify_complete.py` reconstructs and checks each relevant cap polynomial rather than assuming the formula. Delsarte positivity and the polynomial-bound family are established tools, not claimed novel.

The named previous upper bound is 0.000915119176, from Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . The unchanged `comparison_certificate.json` is `certificates/n026.json` from https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb . It is a benchmark only, not an input assumption. Optional source replay uses `python ../prior/taylor.py comparison_certificate.json`.

`independent_audit.json` preserves the separately produced full audit report of the exact certificate SHA; it is supporting provenance, not a substitute for the executable proof. No worldwide-priority or optimality claim is made.
