# 47D lattice-only upper bound

**Density ≤0.000000591010.** Immutable certificate SHA256: `d1ea47251eb7d77ab1a7d3b119be1e8ae6515a0cd25db8c511e41b1495031ef8`.

This result uses global pair-sum angular coupling plus modulo-two class budgets. It applies only to lattice packings. The comparison upper **0.0000005937695380** is the earlier unrestricted-packing result in this repository, which also applies to lattices; it is not attributed to Zubrilina. See `COMPARISON.json` for both Cohn rows and the historical proof path. The comparison is not an input to this proof.

## Full exact replay

```sh
python run_all.py
```

Python3.10+ and the public SymPy package are required. The public `python-flint` package is optional acceleration. When it is available, the runner selects the exact Fourier Bernstein path; otherwise exact SymPy root counting is used and can be much slower. Tail Bernstein arithmetic has a standard-library Fraction fallback. The runner pins the certificate SHA256, dimension and reported ceiling, then executes the complete verifier with an explicit zero-exit check. Do not use `python -O`.

The replay checks all5 continuous angular cuts, the generalized monomial projective caps (r=0…8), 120 radial intervals, all empty and conditioned class budgets, both global half-line signs, and odd-dimensional density normalization using a rational lower bound for pi. No numerical sign tolerance or sampled-domain claim substitutes for these checks. Optional independent producer implementation: `python verify_sympy.py certificate.json`.

`METHOD.md` and `ANGULAR_PROOF.md` assign classes to **half-open** bins; closed intervals are used only for conservative continuous envelopes. `GENERAL_PROJECTIVE_MONOMIAL_CAP.md` proves the complete cap family. This extends Zubrilina's mod-2 class framework ([paper](https://arxiv.org/abs/2610.35830)) by global pair-sum angular constraints. Standard Delsarte identities are not claimed as new. No prior-density, near-frame arithmetic, optimality or exhaustive worldwide-priority assertion is needed.

The saved earlier audit report is provenance, not evidence accepted in place of rerunning the verifier. All required code and mathematical proof notes are local to this directory; no private filesystem imports or remote data retrieval are used.
