# D4 certificate data

Inputs for the bound Δ₄ < 0.629521681653 in [`../paper.pdf`](../paper.pdf). Numbers are exact rationals stored as strings `"p/q"`.

| File | What it is | Paper |
|---|---|---|
| `v7/scalar.json` | The scalar function: coefficients of `F` (degree 49) and `G`, with `f(x) = exp(-2|x|²) F(2|x|²)` | Section 3 |
| `v7/kernel.json` | The harmonic Gram kernel: radial breaks, degree-six Bernstein pieces, and the 19 blocks as integer rays with a common denominator | Section 3 |
| `v7/origin.json` | Origin coupling, pointed caps, tents and allocations, conditional rays and count polynomials, spherical energies and the base objective | Section 3 |
| `v7/graph_energies.json` | Weights of the graph energies | Table 2 |
| `census_floors.json` | Exact raw objective, the 124 group floors after the finite-count refinement (`profile`), and the refined bound 0.629526945573 | Sections 4 and 5 |
| `kernel-box/kbox_j1.json`, `kbox_j5.json` | The kernel-box subdivision trees proving `-K ≥ λ q` on `B₁` (λ = 3/10⁶, 286 leaves) and `B₅` (λ = 1/10⁶, 4663 leaves) | Section 6 |
| `kernel-box/ledger.json` | The final ledger: raised floor, binding groups and the bound | Section 6 |

The `status` field in `v7/origin.json` is a historical label from the producer. Some JSON files also contain old absolute paths as provenance text. No program opens them.

## How these files are checked

[`../verify`](../verify) holds every program that checks these files, with the further certificates and code they need. One command reruns all of it:

    python3 verify/verify_all.py --out /tmp/d4-check --mode full

The stage 1 checker recomputes the signs and the raw objective from `v7/`. Stages 2 and 3 recompute the count census and reproduce `census_floors.json` exactly. Stage 5 re-evaluates every leaf of the two `kernel-box` trees and reproduces `ledger.json` exactly. See [`../verify/README.md`](../verify/README.md).

## Projected-code certificate

`projected_eight.json` supplies the exact nonnegative Legendre coefficients
used for the projected code bound in stage 4. The replay program
`verify/stage4_graph_slack/src/verify_projected_eight.py` reconstructs the
polynomial, checks its sign on the complete interval
`[-1, 223/735]` by rational Bernstein subdivision, and verifies the objective
is strictly below nine. Its numerical search provenance is retained in the
JSON as history only; it is not used as proof input.
