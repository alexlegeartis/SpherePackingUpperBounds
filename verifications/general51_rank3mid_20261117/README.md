# Stronger unrestricted sphere-packing upper bound in dimension 51

**Every sphere packing in R^51 has density less than 1.38612634 × 10^-7.**

The degree-29 angular kernel now has three radial profiles with powers 12,18,24, 30 exact 3-by-3 PSD harmonic blocks and a 4-by-4 origin block. Six tensor-Bernstein controls prove support-crossing angular signs. The compact scalar correction has degree 48, distinct from the degree-71 Gaussian/Fourier polynomials. No lattice or antipodal premise is used.

The complete compact triangle domain was independently replayed at 512 bits: 78,425 nodes and 39,213 accepted leaves. The independent replay sums all nine matrix entries, rather than the producer's six grouped terms. Exact PSD/principal minors, origin Schur complement, all six angular controls, twelve edges, 633 radial annuli, infinite tail, code-cap budgets and odd-dimensional Gamma/pi normalization pass. Frozen proof inputs are byte-identical to the independent audit.

## Reproduce

Install `python-flint` and `sympy`; assertions must remain enabled.

```sh
python run_all.py
python run_all.py --full
```

The runner validates and safely extracts the self-contained archive, checks every file digest, and runs the verifier. Quick mode checks exact bindings, feature/control identities, PSD, angular signs, edges, Fourier/tail conditions, code caps, normalization and the polynomial-composition validation. Full mode also regenerates all 633 radial annuli and independently replays the entire partition.

[proof.tar.xz](proof.tar.xz) includes the complete proof source, original audit logs/receipts, `COUNTING_THEOREM.md`, `THEOREM_SCOPE.md` and `AUDIT.md`. `STORED_FILES.json` lists digests. Standard positive-definite kernel and Poisson theorems are mathematical inputs, not formally verified here. Shared scalar algorithms were reviewed and rerun, not wholly independently reimplemented.

## Comparison

This improves the previous repository certificate <1.39812853e-7 by **0.85844683%**; that rank-two proof remains preserved. Against the [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini numerical spinless-bootstrap estimate](https://doi.org/10.1007/JHEP12%282020%29066) ≈1.5321362440043e-7 (N=457), the decrease is **9.5298250776%**. [Original dataset](https://hdl.handle.net/1721.1/125646). The external comparison is numerical, not a certified prior upper bound, and does not establish exhaustive worldwide priority.
