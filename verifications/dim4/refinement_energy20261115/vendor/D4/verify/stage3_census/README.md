# Stage 3: finite-count refinement of the census floors

This stage reruns "Finite counts in the spherical energies" in Section 5 of
[`../../paper.pdf`](../../paper.pdf). It raises the 124 group floors of
stage 2 and produces [`../../data/census_floors.json`](../../data/census_floors.json).

## What it proves

Notation is that of [stage 2](../stage2_origin/README.md): an admissible
neighbourhood `X`, its energy `E(X)`, its counts `n_1..n_5`, `N = n_5`, and its
group `(N, j)` with `j` the first index such that `n_j = N`.

Stage 2 bounds the spherical part of `E(X)` with three repaired spherical
certificates (index 0: `N = 24`, radii at most `R1`; index 1: every
`N <= 26`; index 2: `N = 25`, radii at most `R4`). Each certificate proves
`P <= eta` for its pair polynomial and `Gamma <= 3 epsilon` for its cyclic
polynomial. Summing the unrepaired anchored certificates over a spherical code
of `N` points gives

```
H_g <= N B_raw + eta N(N-1) + epsilon N(N-1)(N-2).
```

The repaired bound used in stage 2 charges `(M-1)(epsilon H + eta)` per point,
with `M = N` for the fixed-count certificates and `M = 26` for the bounded
one. The stage establishes, with exact arithmetic:

1. For each of the 124 groups `(N, j)` and each certificate that applies to
   it, the spherical surplus can rise to
   `N Bbar - W (N B_raw + eta N(N-1) + epsilon N(N-1)(N-2))`. Within a group
   this is constant, so the new floor is the stage 2 floor minus the old
   surplus plus the largest available surplus (old or new). Hence
   `E(X) >= new_floor(N, j)` for every admissible `X` in group `(N, j)`.
2. The 124 new floors equal `profile[].new_floor` in
   `data/census_floors.json`. Their minimum is
   `a_+ = 0.003783766971828967494...`, in group `(24, 1)` at counts
   `(24, 24, 24, 24, 24)`. Exactly 7 groups have floor below
   `b = 8 (U_raw - 629/1000)`.
3. With the broad-graph slack `3/8000000` (the lemma in Section 5 that lowers the capacity 1358 to 1208,
   checked in stage 4), Theorem 2.1 gives density at most
   `U_raw - 3/8000000 - a_+/8 = 0.6295269455720915060...`, below
   `629526945573/10^12`.

The derivation of the inequality for `H_g` is in
[`finite_count_gain.md`](finite_count_gain.md) and Section 5 of the paper. This
stage checks its inputs and recomputes the arithmetic.

## Command

From this folder, with Python 3.12 or later and `numpy` (see
`requirements.txt`):

```
python3 run.py --out /path/to/new/dir
python3 run.py --out /path/to/new/dir --stage2 /path/to/stage2/out
```

`--out` must not exist. Without `--stage2` the stage reads
`../stage2_origin/reference/`, the committed stage 2 output. With `--stage2`
it reads a fresh stage 2 output directory and requires its `RESULT.json` to
say `PASS`. Do not run Python with `-O`.

## Output

Exit code 0 if and only if `DIR/RESULT.json` has `"status": "PASS"`. The run
fails unless the fresh `DIR/FINITE_COUNT_ACCEPTANCE.json` equals
`data/census_floors.json` in every field except `source` (the input folder).
In particular the stage compares all 124 group floors, `new_beta`,
`exact_raw_objective`, `exact_objective` and `strict_upper_bound` as exact
rationals. Expected values:

| Field | Value |
|---|---|
| `new_beta` | `a_+`, about 0.003783766971828967494 |
| `old_beta` | stage 2 `beta`, about 0.0036901023146053 |
| `exact_raw_objective` | `U_raw`, about 0.6300002914435702 |
| `exact_objective` | about 0.6295269455720915060 |
| `strict_upper_bound` | `629526945573/1000000000000` |
| `group_floors` | 124 |

`DIR/finite_count_gain.json` is the output of the separate arithmetic script
`multiplicity/finite_count_gain.py`, rerun on the same input. The checker
requires both scripts to find the same `a_+`.

Measured on an Apple M2 Pro (macOS) with Python 3.14.6: under 1 second, peak
RSS 49 MB. A full chain is stage 1, then stage 2 with `--stage1`, then this
stage with `--stage2`.

## Recomputed or replayed

| Claim | How |
|---|---|
| Integer Gram rays and the fixed-count restriction of the three certificates | Recomputed (`fixed_cardinality_spherical_algebra.gram_blocks`) |
| `B_raw` of each certificate | Recomputed from the certificate and compared with the stage 2 angular report |
| Same spherical polynomial `gbar`, weight `W` and `Bbar` as `origin.json` | Recomputed |
| Repair accounting `B_repaired = B_raw + (M-1)(epsilon H + eta)` | Recomputed |
| Pair and triangle sign proofs of the three certificates | From stage 2 `ANGULAR_0..2.json`. Recomputed when `--stage2` points to a fresh run; read from `reference/` otherwise |
| The 124 old group floors and old surpluses | From stage 2 `LOCAL_ENUMERATION.json`, same rule |
| The 124 new floors, `a_+`, the objective and outward rounding | Recomputed, twice: by the checker and by `finite_count_gain.py` |
| Broad-graph slack `3/8000000` | Constant here; proved in stage 4 |

The output fields `inherited_sign_checks_rerun` and `new_count_census` are
kept as the original checker writes them (`false`), so the output stays equal
to `data/census_floors.json`. They describe this checker alone. When run with
`--stage2` on a fresh stage 2 output, the signs and the census were recomputed
upstream; `RESULT.json` records this as `"stage2_fresh": true`.

## Trust base

- The written derivation in Section 5 of the paper and Theorem 2.1.
- Stage 2 (and through it stage 1).
- Stage 4 for the slack `3/8000000`.
- Python `fractions.Fraction`, integer arithmetic and NumPy integer arrays.

## Files

- `run.py`: the driver.
- `multiplicity/finite_count_gain.py`: the separate arithmetic script,
  unchanged.
- `finite_count_gain.md`: the written derivation of the finite-count
  inequality.
- `original/verify_finite_count.py`: the private checker, for `diff`.
- The checker also reads, from the stage 2 folder: `source/target625/origin_completion/FINAL_ORIGIN_REGISTRY.json`,
  `source/target625/geometry/fixed_cardinality_spherical_algebra.py` (run),
  three certifier modules (identity only) and the three certificate folders;
  and `../../data/v7/origin.json`. Each is pinned by the digests in the stage 2
  `ACCEPTANCE.json`.
- `EXPORT_MANIFEST.json`: path, size and SHA-256 of every file in this folder.
- `CHANGES.md`, `requirements.txt`.
