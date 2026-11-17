# Unrestricted sphere-packing upper bound in dimension 51

**Every sphere packing in R^51 has density less than 1.39812853 × 10^-7.**

This degree-29 compact cyclic harmonic-matrix certificate uses two radial features, exact rational PSD blocks and a repaired degree-71 Fourier polynomial. It assumes neither lattice structure nor antipodality nor mod-2 restrictions.

The complete compact triangle domain was independently replayed at 512-bit precision: 98,507 nodes and 49,254 accepted leaves. All angular signs, the infinite scalar tail, 633 contiguous radial annuli, exact code-cap budgets and safe density rounding were verified. Odd-dimensional normalization was independently checked using exact squared-density inequalities with a Machin-formula enclosure of pi, alongside high-precision Arb Gamma evaluation.

## Reproduce

Install `python-flint` and `sympy`; keep assertions enabled.

```sh
python run_all.py
python run_all.py --full
```

The runner checks the archive digest, safely extracts to a temporary directory, validates every unpacked file, and executes the verifier. The first command checks PSD/bindings, angular corners, all twelve edge intervals, Fourier/tail conditions, exact cap budgets and normalization. Full mode additionally regenerates all 633 annuli and independently replays the complete three-dimensional partition; it can take tens of minutes on one CPU.

[proof.tar.xz](proof.tar.xz) contains the full proof source, exact inputs, trace, seven independent audit logs, AUDIT.md and completion receipts. `STORED_FILES.json` lists their digests. All frozen audited inputs are byte-identical to the independent audit. The scalar arithmetic algorithms were reviewed and rerun, not wholly independently reimplemented. This is the validated geomfix kernel, not an earlier rejected numerical candidate or an unproved rank-three experiment.

## Comparison and scope

The previous repository certificate was <1.52983244e-7, so this is an 8.60904152% reduction. That earlier proof remains preserved. The [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini (2020) numerical spinless-bootstrap benchmark](https://doi.org/10.1007/JHEP12%282020%29066), truncation N=457, is approximately 1.5321362440043e-7; this result lies 8.7464619761% below that numerical value. [Original dataset](https://hdl.handle.net/1721.1/125646). The external value is not a certified sphere-packing upper bound; this comparison does not certify it or establish exhaustive worldwide priority.
