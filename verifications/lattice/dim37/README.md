# 37D lattice-only upper bound

**Density ≤0.000020580470.** SHA256: `ee54744cac8fdb49c1764ebf2bde57e0b01a0cb967992282ab64779d0cc1ee1a`.

This result uses cross-class pair-sum angular coupling plus mod-2 profile bounds. It is lattice-only, not an improvement to the unrestricted-packing bound. The predecessor in Cohn's table is Alex Kravatskiy's unrestricted upper bound 0.0000209426592470, which also applies to lattices; it is not attributed to Zubrilina. See COMPARISON.json.

Run `python run_all.py` with Python 3.10+ and SymPy; public `python-flint` is optional acceleration. The complete independent verifier reconstructs all continuous angular bounds and projective caps, Fourier sign by exact root count, radial tail by exact Bernstein subdivision and a positive coefficient shift, all 120 radial intervals and 97 class budgets, and odd-dimensional normalization with a rational lower bound for pi. Do not use `python -O`. Optional producer replay: `python verify_sympy.py certificate.json`.

The proof builds on Nina Zubrilina's mod-2 class framework (*Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830), using additional global pair-sum angular inequalities. It has no prior-density or near-frame arithmetic assumption. METHOD.md and ANGULAR_PROOF.md explicitly assign classes to half-open bins; continuous envelope tests on closed intervals are conservative. Three angular inequalities here use one source window. Projective cap identities are standard Delsarte polynomial bounds, not claimed as new. No exhaustive worldwide-priority or optimality claim is made.
