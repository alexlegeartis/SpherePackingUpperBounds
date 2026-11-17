# Unrestricted sphere-packing upper bound in dimension 52

**Every sphere packing in R^52 has density less than 9.7440293 × 10^-8.**

This degree-29 compact cyclic harmonic-matrix certificate uses two radial features, exact rational PSD blocks, and a repaired degree-71 Fourier polynomial. It assumes neither lattice structure nor antipodality nor mod-2 restrictions.

The complete compact triangle domain was independently replayed at 512-bit precision: 79,185 nodes and 39,593 accepted leaves. Full angular signs, the infinite scalar tail, 629 contiguous radial annuli, exact code-cap budgets and safe density rounding were verified. See AUDIT.md inside [the proof archive](proof.tar.xz) for the theorem, coverage, boundary handling and precise independence scope. This is the fresh geomfix kernel; earlier rejected experimental kernels are not included.

## Reproduce

Install `python-flint` and `sympy`; keep Python assertions enabled.

```sh
python run_all.py
python run_all.py --full
```

The first command checks package integrity, exact PSD and bindings, angular corner signs, all twelve edge intervals, Fourier/tail conditions, code caps and exact density rounding. The full command additionally regenerates all radial annuli and independently replays every accepted three-dimensional box. It can take tens of minutes on one CPU.

All audit source files, arithmetic inputs, partition trace, seven original audit logs and completion hashes are included in `proof.tar.xz`. `run_all.py` checks its digest, safely extracts it to a temporary directory, verifies every unpacked file, and runs the checks there. `STORED_FILES.json` lists the unpacked content digests. The scalar arithmetic algorithms were reviewed and rerun, not wholly independently reimplemented; the separate full-domain replayer and normalization checker provide additional independent checks.

## Comparison and scope

The previous repository certificate was <1.06516626e-7, so this is a 8.52104816% reduction. That earlier proof remains preserved. The [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini (2020) numerical spinless-bootstrap benchmark](https://doi.org/10.1007/JHEP12%282020%29066), truncation N=466, is approximately 1.0661233909836e-7; this result lies 8.6031749945% below that numerical value. [Original dataset](https://hdl.handle.net/1721.1/125646). The external value is not a certified sphere-packing bound, and this comparison does not certify it or establish exhaustive worldwide priority.
