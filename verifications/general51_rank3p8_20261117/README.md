# Stronger unrestricted sphere-packing upper bound in dimension 51

**Every sphere packing in R^51 has density less than 1.37781396 × 10^-7.**

This three-feature compact cyclic kernel uses radial powers 8,12,16, thirty exact 3-by-3 PSD harmonic blocks, and a 4-by-4 origin block. Six tensor-Bernstein controls establish support-crossing signs. The compact correction has degree 32; the separate Gaussian and Fourier polynomials have degree 71. No lattice or antipodal assumption is used.

The complete independent 512-bit replay covers 91,139 nodes and 45,570 accepted leaves, summing all nine matrix entries rather than the producer's symmetry-grouped expressions. Exact PSD/origin conditions, all angular controls and edges, Fourier positivity, 633 radial annuli and infinite tail, code-cap budgets and independently reconstructed odd-dimensional normalization pass. Frozen inputs remained byte-identical throughout the audit.

## Reproduce

Install `python-flint` and `sympy`; keep assertions enabled.

```sh
python run_all.py
python run_all.py --full
```

The runner safely extracts the archive after verifying its digest, validates every stored file and runs the checks. Full mode additionally regenerates all radial annuli and independently replays the entire three-dimensional partition. Quick mode checks exact bindings/PSD/profile identities, six angular controls, all edges, Fourier/tail conditions, cap budgets, normalization and polynomial-shift validation.

[proof.tar.xz](proof.tar.xz) contains the proof source, exact data, partition trace, seven audit logs and receipts, `COUNTING_THEOREM.md`, `THEOREM_SCOPE.md` and `AUDIT.md`. Standard kernel/Poisson theorems are mathematical inputs, not formally verified. Shared scalar algorithms were reviewed and rerun, not wholly independently reimplemented. Historical producer status notes are explicitly named PRODUCER_* and do not override the completed independent audit.

## Comparison

The decrease from the preceding repository certificate <1.38612634e-7 is **0.59968415%**; that proof remains preserved. Relative to the [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini numerical spinless-bootstrap estimate](https://doi.org/10.1007/JHEP12%282020%29066) ≈1.5321362440043e-7 (N=457), the decrease is **10.0723603797%**. [Original dataset](https://hdl.handle.net/1721.1/125646). This external comparison is numerical, not a certified prior upper bound or an exhaustive worldwide-priority guarantee.
