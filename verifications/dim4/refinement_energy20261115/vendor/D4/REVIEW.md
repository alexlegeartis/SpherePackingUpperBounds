# D4 publication verification review

## Release boundary

The published result is the non-Lean computer-assisted proof described in
`paper.tex`, using the exact inputs in `data/` and the five replay stages in
`verify/`. The Lean formalization is unreleased and is not needed for this
claim; see [`FORMALIZATION_STATUS.md`](FORMALIZATION_STATUS.md).

## Reproduction commands

Run from this D4 folder with Python 3.12 or later:

```sh
python3 -m pip install -r verify/requirements.txt
python3 verify/verify_all.py --out /tmp/d4-full --mode full --workers 2
python3 -m unittest discover -s verify -p 'test_*.py' -v
python3 verify/kernel_independent/negative_controls.py --out /tmp/d4-independent-controls
python3 verify/stage5_kernelbox/negative_controls.py --work /tmp/d4-box-controls --workers 2
```

All output directories must be new. Quick mode omits the 112 cyclic
subdivisions and is not a fresh global-certificate replay. Do not use `-O` or
`PYTHONOPTIMIZE`; inherited checkers use assertions.

## Replay scope

The full driver checks every exported manifest, runs the five stages in order,
passes fresh outputs downstream, reruns the independent kernel reconstruction
and exact kernel link, and recomputes the final ledger in exact arithmetic.
Stage 4 also replays `data/projected_eight.json` with an exact rational
Legendre/Bernstein checker. It derives the degree bound from that checker and
fails closed if the certificate is missing, altered, incomplete, or does not
give an objective below nine.

The recorded full replay established the following decisive coverage:

| Stage | Result |
|---|---|
| V7 signs | all 112 cyclic branches |
| Origin credit | 165,503 count vectors |
| Census | 124 group floors |
| Graph slack | annular count 151, projected degree 8, ordered capacity 1208 |
| Far-pair boxes | 286 and 4,663 leaves |
| Independent kernel | 112 branches, 7 pair pieces, exact coefficient comparisons |

The final exact objective remains strictly below `629521681653/10^12`. The
stage-4 improvement remains `3/8000000`; integrating the projected certificate
changes the dependency boundary, not the numerical result.

## Trust and attribution

The exact checks use Python integer/rational arithmetic, FLINT exact
polynomials, Arb outward-rounded balls, and a C++20 subdivision engine built
from source. The written local-energy reductions, spherical positive
definiteness, interval lower-envelope soundness, and density-transfer argument
remain explicit proof inputs. The classical nine-point theorem is retained as
historical attribution and an alternative written route; the available replay
uses the supplied projected-code certificate instead.

Independent assembly and symbolic-link checks compare implementations exactly
but do not constitute formal semantics verification. Hash manifests identify
the exported files only.
