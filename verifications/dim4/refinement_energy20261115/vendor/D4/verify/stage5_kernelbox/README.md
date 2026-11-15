# Stage 5: far-pair kernel box and the final ledger

This stage proves the far-pair increment H_R >= lambda_j max(Q_N, 0) and combines it with the census floors into the final bound: every sphere packing in R^4 has density at most U_new < 0.629521681653.

## Statements

Notation: a = 133/250, q(t) = t^2 (t^2 - a^2), K is the raw V7 pair kernel of the base certificate (the polarization of `src/raw_v7.py`), and R_1 = 501473/500000, R_5 = 102597/100000 are kernel breaks. B_j is the set of pairs (y, z) with |y|, |z| in [1, R_j] and cos(y, z) = t in [-1, -a].

1. CERTIFIED (interval arithmetic). -K(y, z) >= (3/10^6) q(t) on B_1 and -K(y, z) >= (1/10^6) q(t) on B_5. Proof: `proofs/kernel_box_method.md`. Data: `../../data/kernel-box/kbox_j1.json` (94 roots, 286 leaves) and `kbox_j5.json` (2350 roots, 4663 leaves).
2. PROVED (written). For every anchor whose N near neighbours lie in [1, R_j], H_R >= lambda_j max(Q_N, 0) with Q_N = N (13561 N - 179244)/250000, lambda_1 = 3/10^6 and lambda_j = 1/10^6 for j = 2..5. Proof: `proofs/far_pair_lemma.md`.
3. PROVED (written). With K_R = K + H_R, a floor a on every anchored energy E_R gives density <= U_raw - a/8. Proof: `proofs/conditional_theorem.md`.
4. COMPUTED EXACT. a_new = min over the 124 census groups g = (N, j) of floor_g + lambda_j max(Q_N, 0). The minimum is attained only at (24, 1), a_new = a_plus + (3/10^6) Q_24 = 0.0038258783318289675... Then

       U_new = U_raw - a_new/8 - 3/8000000 = 0.629521681652091506... < 629521681653/10^12.

## Run

    pip install -r requirements.txt
    python3 run.py --out OUT [--workers 4]

Optional: `--census DIR` and `--stage2 DIR` take fresh outputs of stages 3 and 2. Each folder must hold a `RESULT.json` with status PASS for that stage, or this stage fails. The census file is `census_floors.json` or `FINITE_COUNT_ACCEPTANCE.json` (same schema); the stage 2 file is `ACCEPTANCE.json`. The defaults are `../../data/census_floors.json` and `../stage2_origin/reference/ACCEPTANCE.json`. `--lam1` and `--lam5` exist only for negative controls.

Exit 0 and `OUT/RESULT.json` with `"status": "PASS"` mean every step passed. A failed assertion, nonzero child exit, missing output or timeout (1 hour per child) gives `"status": "FAIL"` and exit 1.

Expected decisive values in `RESULT.json`:

| Key | Value |
|---|---|
| `kernel_box.j1_roots`, `j1_leaves_replayed` | 94, 286 |
| `kernel_box.j5_roots`, `j5_leaves_replayed` | 2350, 4663 |
| `kernel_box.j1_min_leaf_lower`, `j5_min_leaf_lower` | 2.2e-8, 1.6e-9 |
| `kernel_box.kformula_vs_raw_v7_points` | 72 |
| `kernel_box.q_gegenbauer` | 1/16, 14593/125000, 13561/250000 |
| `raw_v7_spot_checks` | additivity true, 88 band points, 0 failures |
| `binding_group`, `Q_binding` | [24, 1], 43866/3125 |
| `next_group` | [26, 5] |
| `graph_saving` | 3/8000000 |
| `U_new_decimal` | 0.629521681652092 |
| `ceiling` | 629521681653/1000000000000 |
| `ledger_equals_committed` | true |

`a_plus`, `a_new`, `U_raw`, `U_new` and `next_raised_floor` are exact fractions with several hundred digits.

Runtime: 285 s to 332 s with 4 workers. Peak summed RSS of all processes: 0.18 to 0.29 GB. Measured on an Apple M2 Pro with Python 3.14.6 and python-flint 0.9.0. The tree replay is about 150 s, the kformula check 50 s and the raw_v7 spot checks 90 s.

Negative controls (about 10 minutes):

    python3 negative_controls.py --work SCRATCH [--workers 2]

It copies the D4 folder twice into SCRATCH and breaks one thing in each copy: (1) it replays `kbox_j1.json` at lambda_1 = 4/10^6, which is false at the dip; (2) it merges two sibling leaves of the dip root in `kbox_j5.json`. Both must end with status FAIL and a "negative leaf" error.

## What it checks

`run.py` runs five scripts as child processes and then checks the ledger itself.

1. `src/check_kernelbox.py`, five gates.
   - Gate 1: the Gegenbauer expansion of q with nonnegative coefficients, s_j = 1 - 1/(2 R_j^2) <= a for j = 1..5, and 2 + 2a > 49/16. Exact.
   - Gate 2: `src/validate_kformula.py`. The explicit kernel formula `src/kformula.py` must meet the independent raw_v7 enclosure at 72 rational points and q0 at 9 radii.
   - Gate 3: kcert point boxes reproduce kformula at 6 points.
   - Gate 4: replay of both trees. The roots must tile B_j exactly, every split must be a midpoint bisection, and every leaf lower bound of -K - lambda q is recomputed by `src/kcert.py` and must be >= 0. No stored leaf value is used.
   - Gate 5: `src/ledger.py` with the certified lambdas.
2. `src/hrlemma_check.py`: the exact facts of the far-pair lemma (Gegenbauer identity, caps, far-pair distance, Q_N table).
3. `src/coord_check_kernelbox.py`: raw_v7 only, without kformula or kcert. Additivity of the anchored energy on four random rational stars, and -K - (3/10^6) q >= 0 at 88 rational points of B_1.
4. `src/ledger.py --out`: the recomputed ledger must equal `../../data/kernel-box/ledger.json` key by key.
5. `src/exact_values.py`: a_new, U_new and the next floor (26, 5) as exact fractions. It also checks that U_raw in the census file equals U_raw in the stage 2 file.
6. `run.py` recomputes a_new from `census_floors.json` with its own code. It checks the 124 distinct groups, U_raw against stage 2, graph_saving = 3/8000000, U_new < 629521681653/10^12, and the 10^-12 grid ceiling.

Recomputed: every leaf bound, the kformula and raw_v7 comparisons, the far-pair facts, a_new, U_new and the ceiling. Replayed rather than searched: the split trees come from `src/kbox_run.py` and are checked, not regenerated. To regenerate them (about 40 s and 290 s with 4 workers):

    python3 src/kbox_run.py --j 1 --lam 3/1000000 --workers 4 --out TREE1.json
    python3 src/kbox_run.py --j 5 --lam 1/1000000 --workers 4 --out TREE5.json

## Trust base

- python-flint 0.9.0 (FLINT and Arb ball arithmetic), at 256 bits. All box bounds are Arb enclosures. A comparison that Arb cannot decide counts as failure, except where the code takes the safe branch documented in `src/kcert.py`.
- Python `fractions.Fraction` for every rational step (raw_v7, the ledger, the identities).
- The written proofs in `proofs/`.

## Inputs and hypotheses from other stages

- `../../data/v7/*.json` (stage 1): the base certificate. Its signs, including the cyclic sign Gamma <= 0, are verified in stage 1, not here.
- `../stage2_origin/reference/ACCEPTANCE.json` (stage 2): `exact_source_objective` = U_raw.
- `../../data/census_floors.json` (stages 2 and 3): the 124 group floors `new_floor`, the group definition, `new_beta` and `graph_saving`.
- Stage 4: the graph saving 3/8000000.

## Not included

- The run receipts and logs of the original tree search. The trees are replayed instead.
- Independent third-party review. The kernel box had author checks only before this package.

## Files

- `run.py`, `negative_controls.py`: drivers.
- `src/`: the checker code. `paths.py` holds all input and output paths. `CHANGES.md` lists every edit to the original code.
- `proofs/`: `kernel_box_method.md`, `far_pair_lemma.md`, `conditional_theorem.md`.
- `MANIFEST.json`, `requirements.txt`.
