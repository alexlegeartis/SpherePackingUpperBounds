# Stage 4: broad graph slack

This stage proves that the broad antipodal graph of the base certificate has at most 1208 ordered edges, not 1358. That lowers the density bound by exactly 3/8000000.

## Statement

Take any packing of unit-separated centres in R^4 and any centre x. Its vertices are the neighbours y with |y - x| in [33/28, 7/4]. Two vertices are joined when the cosine of their directions from x lies in [-1, -7/8]. Both intervals are closed. Let e be the number of ordered edges.

- PROVED: e <= 1208. There are at most 151 vertices (exact polynomial certificate), and each vertex has degree at most 8. The projected S^2 code bound is replayed from `data/projected_eight.json` by an exact rational Legendre/Bernstein checker.
- Consequence: the graph energy delta (1358 - e) with delta = 1/50000000 is at least 150 delta = 3/1000000 at every anchor. The census floors use only that it is at least 0. So every floor rises by 3/1000000, and the density bound drops by 3/1000000 / 8 = 3/8000000.

The full proof is in [proofs/graph_slack.md](proofs/graph_slack.md).

## Projected-code certificate

The checker reconstructs the supplied nonnegative Legendre combination and
proves its sign on the entire closed interval `[-1, 223/735]` by exact
rational Bernstein subdivision. Its exact objective is
`89804130521/9999900000 < 9`, so the addition-formula argument gives
`A(3,223/735) <= 8`. The classical nine-point result remains a cited
alternative and historical attribution; it is not an input to this replay.

The exact geometric gap used by the projection is still
`1/3 - 223/735 = 22/735 > 0`; it is checked independently by `gap_check.py`.

## Run

Use CPython 3.12 or later (the release-wide requirement; stage 4 itself uses
only the standard library):

    python3 run.py --out OUT

Optional: `--census DIR` and `--stage2 DIR` take fresh outputs of stages 3 and 2. Each folder must hold a `RESULT.json` with status PASS for that stage, or this stage fails. The census file is `census_floors.json` or `FINITE_COUNT_ACCEPTANCE.json` (same schema); the stage 2 file is `ACCEPTANCE.json`. The defaults are `../../data/census_floors.json` and `../stage2_origin/reference/ACCEPTANCE.json`.

Exit 0 and `OUT/RESULT.json` with `"status": "PASS"` mean every check passed. Any failed assertion, nonzero child exit or timeout gives `"status": "FAIL"` and exit 1.

Expected decisive values in `RESULT.json`:

| Key | Value |
|---|---|
| `annular_count` | 151 |
| `exact_count_upper` | 151196972764393/999999000000 |
| `bernstein_nodes`, `bernstein_leaves`, `bernstein_depth` | 93, 47, 12 |
| `projected_eight.degree` | 27 |
| `projected_eight.bernstein_nodes`, `projected_eight.bernstein_leaves`, `projected_eight.bernstein_depth` | 23, 12, 7 |
| `projected_eight.objective` | 89804130521/9999900000 (< 9) |
| `projected_eight.count_bound` | 8 |
| `corner_cosines` | 697/1089, 41/49, 41/49 |
| `projected_cosine_max` | 223/735 |
| `gap_to_one_third` | 22/735 |
| `ordered_edge_capacity` | 1208 |
| `graph_origin_charge` | 4769/25000000 |
| `slack_in_energy` | 3/1000000 |
| `density_saving` | 3/8000000 |

Runtime: 0.2 s. Peak RSS: 26 MB. Measured on an Apple Silicon Mac with Python 3.14.6.

## What it checks

1. `src/verify_annular_count_151.py` with `src/annular_count_151_certificate.json`. It rebuilds the degree-20 polynomial p from its harmonic coefficients, proves p <= 0 on [-1, 41/49] by exact Bernstein subdivision, and computes p(1)/a_0 < 152. Recomputed.
2. `src/verify_projected_eight.py`: the exact Legendre coefficients, full-interval Bernstein sign, endpoint, and objective-to-count conversion. Recomputed.
3. `src/gap_check.py`: Q(7/8, 7/8) = 223/735, 1/3 - 223/735 = 22/735, 7/8 > 41/49, 151 * 8 = 1208, (1358 - 1208)/(8 * 50000000) = 3/8000000. Recomputed.
4. `run.py` itself: the corner cosines of the annulus, the projected cosine bound, the exterior distance 16335/3136 > 49/16, and the place of the graph charge in U_raw. It checks U_raw = objective(`data/v7/origin.json`) + 4769/25000000/8 against `census_floors.json` and the stage 2 `ACCEPTANCE.json`, and it checks `exact_objective` = U_raw - new_beta/8 - 3/8000000 in the census file. Recomputed.

All arithmetic is exact (Python `fractions`). No floating point enters a decision.

## Trust base

- Python 3 standard library: `fractions.Fraction` and integer arithmetic.
- The addition formula for normalized Legendre kernels, as stated in the written proof.
- Stages 2 and 3 for the census floors and for the fact that they use only delta (1358 - e) >= 0.

## Files

- `run.py`: stage driver.
- `src/verify_annular_count_151.py`, `src/verify_projected_eight.py`, `src/verify_exterior_graphs.py` (exact Bernstein helpers), `src/annular_count_151_certificate.json`, `src/gap_check.py`: exact certificate checkers.
- `proofs/graph_slack.md`: the written proof.
- `CHANGES.md`, `MANIFEST.json`, `requirements.txt`.
