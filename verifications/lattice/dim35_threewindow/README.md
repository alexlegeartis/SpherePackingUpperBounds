# 35D refinement: 3 source windows

Certified target density ≤.000037664745. Candidate SHA256: `e6ba7b62cc3a42d5e4e28282b0e69d068536e34768a37b9d9f43c03e2ac367e8`.

This stronger certificate supersedes the previous public bound. Earlier packages remain available unchanged. Full portable independent replay passed.


Run `python run_all.py`. Python 3.10+ and SymPy are needed for the generic independent verifier; `python-flint` is optional acceleration. The 36D default complete verifier uses only the standard library and exact compactified half-line Bernstein sign proofs. Do not run with `-O`.

The generic default checks every angular source prefix, global Fourier sign by exact root count, tail sign by exact rational Bernstein subdivision and positive shift, all radial intervals and class budgets, and exact density. Nested source windows need not be disjoint: valid inequalities combine with nonnegative weights. No prior-density or arithmetic near-frame assumption is used. Optional producer replay: `python verify_sympy.py certificate.json`.

The comparison certificate comes from Nina Zubrilina's `AstroNinja389/lattice-packing-bound`, pinned commit `340a0e93bb95942c2b1b80092be43149504970eb`, accompanying *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830. It is a comparison only, not a proof assumption. The method and angular/projective proof files retain attribution. No optimality or exhaustive priority claim is made.
