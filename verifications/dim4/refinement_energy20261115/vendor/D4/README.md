# A global sphere-packing upper bound in dimension four

The result is **upper density < 0.629521681653 for every packing of congruent balls in R4**. It is an upper bound, not a proof of D4 optimality.

## Contents

- [paper.pdf](paper.pdf) and [paper.tex](paper.tex): the preprint and editable source.
- [data/](data/README.md): exact scalar/kernel coefficients, 124 census floors, and both kernel-box certificates.
- [verify/](verify/README.md): all verification code and inherited certificate dependencies, with two additional kernel cross-checks.
- [FORMALIZATION_STATUS.md](FORMALIZATION_STATUS.md): the status of the unreleased Lean work.
- [REVIEW.md](REVIEW.md): the release replay and trust review for the available proof.
- [PROOF_DEPENDENCIES.md](PROOF_DEPENDENCIES.md): the mathematical dependency chain and the precise boundary between machine checks, written arguments, and cited theorems.

No executable certificate input comes from outside this folder. Historical absolute filenames inside some data files are provenance strings. Standard software dependencies are installed separately.

## Reproduce

Run from this folder with Python 3.12 or later, the pinned packages, and a C++20 compiler:

```sh
python3 -m pip install -r verify/requirements.txt
python3 verify/verify_all.py --out /tmp/d4-full-check --mode full --workers 2
```

Use a new output directory. Do not use `python -O` or `PYTHONOPTIMIZE`: inherited checkers use assertions. `--mode quick` omits the 112 triangle subdivisions and is **not a fresh verification of the global theorem**.

The command above is the complete computer-assisted proof replay. It reads only
`data/` and `verify/` in this folder, plus the listed Python, FLINT, Arb and C++
dependencies. The separate Lean formalization is unreleased; see
[FORMALIZATION_STATUS.md](FORMALIZATION_STATUS.md).

To rebuild the preprint:

```sh
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
```
