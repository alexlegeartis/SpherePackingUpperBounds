# Unrestricted sphere-packing upper bound in dimension 53

**Every sphere packing in R^53 has density less than 6.8841195 × 10^-8.**

This degree-29 compact cyclic harmonic-matrix certificate uses two radial features, exact rational PSD blocks and a repaired degree-71 Fourier polynomial. It assumes neither lattice structure nor antipodality nor mod-2 restrictions.

The complete compact triangle domain was independently replayed at 512-bit precision: 132,889 nodes and 66,445 accepted leaves. All angular signs, the infinite scalar tail, 629 contiguous radial annuli, exact code-cap budgets and safe density rounding were verified. Odd-dimensional normalization was independently checked using exact squared-density inequalities with a Machin-formula enclosure of pi, alongside high-precision Arb Gamma evaluation.

## Reproduce

Install `python-flint` and `sympy`; keep assertions enabled.

```sh
python run_all.py
python run_all.py --full
```

The runner checks the archive digest, safely extracts to a temporary directory, validates every unpacked file, and executes the verifier. The first command checks PSD/bindings, angular corners, all twelve edge intervals, Fourier/tail conditions, exact cap budgets and normalization. Full mode additionally regenerates all 629 annuli and independently replays the complete three-dimensional partition; it can take tens of minutes on one CPU.

[proof.tar.xz](proof.tar.xz) contains the full proof source, exact inputs, trace, seven independent audit logs, AUDIT.md and completion receipts. `STORED_FILES.json` lists their digests. All frozen audited inputs are byte-identical to the independent audit. The scalar arithmetic algorithms were reviewed and rerun, not wholly independently reimplemented. This is the validated geomfix kernel, not an earlier rejected numerical candidate or an unproved rank-three experiment.

## Comparison and scope

The previous repository certificate was <7.40755520e-8, so this is a 7.06624096% reduction. That earlier proof remains preserved. The [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini (2020) numerical spinless-bootstrap benchmark](https://doi.org/10.1007/JHEP12%282020%29066), truncation N=469, is approximately 7.4137674631496e-8; this result lies 7.1441135129% below that numerical value. [Original dataset](https://hdl.handle.net/1721.1/125646). The external value is not a certified sphere-packing upper bound; this comparison does not certify it or establish exhaustive worldwide priority.
