# Stage 2: origin credit and the count census

This stage reruns Section 4 of [`../../paper.pdf`](../../paper.pdf). It
recomputes the energy floor from neighbour counts and the 124 group minima that
stage 3 refines.

## What it proves

Notation follows the paper. `K` is the V7 kernel in
[`../../data/v7/`](../../data/v7/), supported on radii in `[1, 7/4]`. A
neighbourhood `X` is admissible if it is a finite subset of `R^4 \ {0}` and
`X ∪ {0}` is unit separated. Its energy `E(X)` is the anchored energy of
Section 2. The radii are

```
R1 = 501473/500000, R2 = 1005643/1000000, R3 = 50379/50000,
R4 = 50767/50000,   R5 = 102597/100000.
```

Let `n_i(X)` be the number of points of `X` with `|x| <= R_i`. Put `N = n_5`,
and let `j` be the first index in `1..5` with `n_j = N`. The pair `(N, j)` is
the group of `X`. There are 124 nonempty groups.

The stage establishes, with exact arithmetic:

1. Count domain. For every admissible `X`, the vector `n(X)` lies in the set
   `C` of integer vectors with `0 <= n_1 <= ... <= n_5`,
   `n <= (24, 24, 24, 25, 26)`, `(23 - n_1)(n_5 - n_1) >= 0` and
   `(23 - n_2)(n_4 - n_2) >= 0`. The caps come from three spherical code
   certificates for `S^3` (at most 24, 25 and 26 points with cosine at most
   203/400, 103/200 and 21/40). The two implications come from two affine
   latitude proofs. `C` has 165503 vectors. Realizability is not assumed.
2. Group floors. For every admissible `X` in group `(N, j)`,
   `E(X) >= floor(N, j)`. Here `floor(N, j)` is the minimum over the vectors
   of `C` in that group of `max(Q(n), L_1(n), ..., L_12(n))`, plus the three
   integer count energies, plus the spherical surplus that applies to the
   group. `Q` is the exact projection bound (4.1) of the paper. `L_1..L_12`
   are twelve affine floors from the ray decomposition of the degree-zero Gram
   block. The surplus is the largest of three local spherical proofs: `N = 24`
   with all radii at most `R1`, `N = 25` with all radii at most `R4`, and
   every `N <= 26` at radius `R5`. All other kernel components are
   nonnegative and are dropped.
3. Minimum. The smallest of the 124 floors is
   `beta = 0.0036901023146053...`, in group `(24, 1)` at counts
   `(24, 24, 24, 24, 24)`. It replaces the old origin credit
   `0.0000066644185634...`.
4. Bound. With `U_raw = 0.6300002914435702...` from stage 1, Theorem 2.1 of
   the paper gives density at most `U_raw - beta/8 = 0.6295390286542445...`,
   which is below `125907805731/200000000000 = 0.629539028655`.

The exact floors of all 124 groups, with a minimizing count vector for each,
are `first_containing_radius_profile` in `DIR/LOCAL_ENUMERATION.json`.

## Command

From this folder, with Python 3.12 or later and `numpy` (see
`requirements.txt`):

```
python3 -m pip install -r requirements.txt
python3 run.py --out /path/to/new/dir
python3 run.py --out /path/to/new/dir --stage1 /path/to/stage1/out
```

`--out` must not exist. `--stage1` takes the output directory of
[`../stage1_v7/run.py`](../stage1_v7/) in full mode. Then the stage checks
that stage 1 reported `PASS` with all 112 cyclic branches rerun, that its
exact raw objective equals the one recomputed here, that its record holds the
same four V7 data files, and it takes the old credited objective from the
stage 1 composer. Without `--stage1`, it uses the committed stage 1 values in
`reference/INHERITED_FULL_V7_ACCEPTANCE.json`, and `RESULT.json` says
`"stage1_checked": false`. Do not run Python with `-O`: the checkers use
`assert`.

## Output

Exit code 0 if and only if `DIR/RESULT.json` has `"status": "PASS"`.
Expected values (all exact, as `p/q` strings in `RESULT.json`):

| Field | Value |
|---|---|
| `exact_source_objective` | `U_raw`, about 0.6300002914435702 |
| `new_origin_credit` | `beta`, about 0.0036901023146053 |
| `objective` | about 0.6295390286542445 |
| `strict_upper_bound` | `125907805731/200000000000` |
| `count_vectors_checked` | 165503 |
| `group_minima` | 124 |

`DIR` also holds the outputs of the original programs: `ACCEPTANCE.json`,
`PROJECTION.json`, `COUNT_GEOMETRY.json`, `counts/`, `RADIAL_PLANES.json`,
`ANGULAR_0..2.json`, `LOCAL_ENUMERATION.json`, `STAGE1_INHERITED.json`, and
`tree/`, the working copy of the source that was run. The run fails unless
every one of these reports equals its copy in `reference/`, except for timing
fields. Stage 3 reads `DIR` through its `--stage2` option.

Measured on an Apple M2 Pro (macOS) with Python 3.14.6 and NumPy 2.5.2:
100 to 145 seconds, peak RSS 210 MB, one core.

## Recomputed or replayed

| Claim | How |
|---|---|
| Projection bound `Q`: 245 range, 49 origin and 25 inverse equations, exact LDL pivots | Recomputed from `kernel.json` and `origin.json`. The projection certificate supplies the solved matrices; every equation is checked exactly. |
| 165503-vector enumeration of the projection plus count energies (27 profiles) | Recomputed |
| Spherical count caps 24, 25, 26 | Recomputed: the three certificates in `source/research/.../spherical_sdp/` are rechecked by exact Bernstein subdivision |
| Affine latitude implications | Recomputed from the radii and cutoffs, then compared field by field with the records in `origin.json` |
| Nested and conditional count polynomials, nonnegative on `C` | Recomputed, including the 165503-vector check |
| Twelve affine floors `L_i` | Recomputed; each must equal its stored record |
| Three angular surplus proofs (pair and triangle subdivisions) | Recomputed; each must equal its stored `VERIFIED.json` |
| Group floors, the 124 group minima and `beta` | Recomputed by enumerating all 165503 vectors; `beta` must equal `consolidation/FINAL_RESULT.json` |
| Exact objective and outward rounding | Recomputed |
| V7 signs and `U_raw` | Stage 1. Checked from its output with `--stage1`, else taken from `reference/` |

Stored results are used only as expected values that fresh results must
equal. No optimizer, SDP solver or native binary runs.

## Trust base

- The written local-energy theorem (Theorem 2.1 of the paper) and the
  composition in Section 4.
- Stage 1: the V7 sign conditions and `U_raw`.
- Python `fractions.Fraction` and integer arithmetic; NumPy object and int64
  arrays in the Bernstein subdivisions. The angular engine bounds every
  coefficient below 2^59 before casting, so midpoint sums stay inside int64.
  A box is excluded only when an upper bound is strictly negative; an
  unresolved box fails the run.
- The original checker code, unchanged. [`CHANGES.md`](CHANGES.md) lists the
  driver changes, which are path plumbing and the stage 1 hand-off.

## Files

- `run.py`: the driver.
- `PROOF.md`: the written argument that the driver checks.
- `source/`: 32 files of the private source tree, byte for byte, at their
  original relative paths. `source/MANIFEST.json` is the producer's manifest
  of the full tree; the driver checks every file it binds against it. The four
  V7 data files are not copied here: `run.py` places `../../data/v7/*.json`
  into the working copy, and the manifest pins their bytes.
- `reference/`: the accepted output of the private run of 2026-10-20.
- `EXPORT_MANIFEST.json`: path, size and SHA-256 of every file in this folder.
- `CHANGES.md`, `requirements.txt`.

## Not included

- Programs that searched for the certificates (SDP and LP solvers). They are
  not needed: every certificate is checked exactly here.
- Three older angular proofs. The bounded proof here dominates them, and the
  credit check shows that the smaller set suffices.
- `affine_latitude_certificates.json`, cited in `origin.json`. The implications
  are recomputed instead.
- The spherical certificates `spherical_b35ca313fbedf710.json` and
  `spherical_abf37dac4759864f.json`, cited in `origin.json`. They prove the
  base spherical energies nonnegative, which is part of stage 1. This stage
  proves the surpluses with its own three angular certificates.
