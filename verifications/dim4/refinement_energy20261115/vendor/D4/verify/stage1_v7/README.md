# Stage 1: the V7 base certificate

This folder reruns the base certificate of Section 3 of [`../../paper.pdf`](../../paper.pdf).
It recomputes every sign condition from the exact data in [`../../data/v7/`](../../data/v7/)
and outputs the raw objective

```
U_raw = (F(0) + A)/8 = 0.630000291443570...
```

as an exact fraction. Stage 2 and later stages start from `U_raw`. It equals
`exact_raw_objective` in [`../../data/census_floors.json`](../../data/census_floors.json),
and `run.py` checks that equality.

## Command

From this folder, with CPython 3.12 or later, the packages in
`requirements.txt` and a C++20 compiler:

```
python3 -m pip install -r requirements.txt
python3 run.py --out /path/to/new/dir --workers 2          # full: 30 to 50 min
python3 run.py --out /path/to/new/dir --quick              # quick: under 2 min
```

`--out` must not exist. `--workers` is 1 to 4 (default 2); it only sets how
many cyclic branches run at once. `--cxx` picks the compiler (default
`clang++`). Do not run Python with `-O`: the checkers use `assert`.

The exit code is 0 only if every step passed. The output directory holds:

| Path | Content |
|---|---|
| `RESULT.json` | `status` (`PASS` or `FAIL`), `mode`, `branches_rerun` (112 or 0), `raw_objective` as `p/q`, step exit codes and times |
| `record/` | working copy of the exported record, with all checker outputs written by the original programs |
| `record/current/candidates/round9_six_graphs_degree18/VERIFIED.json` | full mode: the coordinator's complete report (all pairs and all 112 branches) |
| `record/replay/v7_full_origin_composition.json` | full mode: the composer's report |
| `checks/` | outputs of the independent checkers |
| `replay.log` | output of every step |

Expected full-mode result (2026-11-13 run, see below):

```
"status": "PASS", "mode": "full", "branches_rerun": 112,
"raw_objective": "369462646652...p/q, 783 digits over 783 digits",
"raw_objective_decimal": "0.630000291443570",
"credited_objective_decimal": "0.629999458391250"
```

## What stage 1 proves

Notation follows the paper. Minimum distance is 1, `f(x) = exp(-2|x|^2) F(2|x|^2)`,
`fhat(xi) = (pi^2/4) exp(-s) G(s)` with `s = pi^2 |xi|^2 / 2`, and `K` is the
V7 kernel, supported on radii in `[1, 7/4]`. Breaks of the seven radial pieces:
`1, 501473/500000, 1005643/1000000, 50379/50000, 50767/50000, 102597/100000, 3/2, 7/4`.
All arithmetic is exact (Python integers, `fractions.Fraction`, FLINT `fmpz_mat`,
or integer intervals with outward rounding).

1. **Scalar signs.** `F` and `G` have degree 49 with rational coefficients.
   The identity `G_j = sum_{k>=j} F_k (-1)^j (k!/j!) C(k+1, k-j)` holds exactly,
   and `G(0) = 1`. Sturm sequences prove `G(s) > 0` for all `s >= 0` and
   `F(t) < 0` for all `t >= 49/8`. So `fhat >= 0`, `fhat(0) = pi^2/4`, and
   `f(x) < 0` for `|x| >= 7/4`.
2. **Gram positivity.** The 19 angular blocks (degrees 0 to 18) are sums of
   `w v v^T / D` with integer `w > 0`, integer vectors `v` and one common
   denominator `D > 0`. There are 628 rays. So every block is positive
   semidefinite by construction.
3. **Origin identity.** With `M` the degree-zero block without its origin
   entry, `M u = v` and `u.v = kappa >= 0` hold exactly (49 equations). The
   physical pair row `pair_origin_v` equals `v` plus the exact linear origin
   shifts of the 98 caps, 4 tents, 2 allocations, 2 signed spherical rays and
   2 conditional radial rays.
4. **Local energies.** Each kernel component has nonnegative anchored energy
   on every finite unit-separated neighbourhood, so condition (1) of
   Theorem 2.1 holds with `beta = 0`. The checks behind this are:
   - the spherical count certificates `s05075_taper12_exact.json` (at most 24
     neighbours with radius at most `50379/50000`), `s0515_count25_exact.json`
     (at most 25 up to `50767/50000`) and `s0525_count26_exact.json` (at most
     26 up to `102597/100000`), each verified by `certify.py`;
   - the nested count quadratic and the two conditional count polynomials,
     checked on every admissible count vector;
   - the affine latitude implications `(23-n_1)(n_5-n_1) >= 0` and
     `(23-n_2)(n_4-n_2) >= 0`, recomputed by `affine_latitude_counts.check`;
   - the two signed spherical energies, replayed from
     `spherical_energy_sources/spherical_b35ca313fbedf710.json` and
     `spherical_abf37dac4759864f.json`. Each carries the uniform repair of
     Section 5 of the paper: its scalar constant is lowered by
     `epsilon (sum q / q_0 - 1) + eta` with `eta = 1/1000000` and
     `epsilon = 1/250000` and `1/500000` respectively;
   - the graph capacities 862, 1102 (weight 0), 1358, 120, 604 and 1880, and
     the moderate capacity 1880 again by an independent Sturm checker.
5. **Pair inequality.** For every `x` with `|x| >= 1`,
   `q(x) = f(x) + 2K(x,0) + K(x,x) <= 0`.
   For `|x| >= 7/4` this is item 1. On each of the 7 closed pieces of
   `[1, 7/4]` an exact Bernstein subdivision proves the pair polynomial
   nonpositive. The polynomial includes the pointed self terms, the nested
   and conditional count pair terms and the signed spherical diagonal
   (`pair_convention` in `origin.json`). The full pair count gauge
   `w = -q/(1-psi)` is attached only after these fresh proofs, so it is
   nonnegative.
6. **Triangle inequality.** For all `x, y` with `|x|, |y|, |x-y| >= 1`,
   `Gamma(x,y) = K(x,y) + K(-x, y-x) + K(-y, x-y) <= 0`.
   The domain splits into 84 triples of radial regions and 28 exterior
   cases where one side exceeds `7/4`. On each case the checker builds the
   cyclic polynomial `P` in the side lengths and the cosine (multiplied by
   the positive factor `(2 r_1 r_2 r_3)^18`, or `(2 r_1 r_2)^18` on exterior
   cases), shifted by the exact constant
   `eta_b = -(nested count term + conditional count polynomial terms + conditional ray terms)`
   for that case, with the signed spherical terms inserted into `P`. It then
   proves, on every box of an exact Bernstein subdivision of the closed case,
   that every upper Bernstein coefficient of `P - eta_b` is at most the
   positive factor times a lower bound for the sum of the nonnegative
   corrections (caps, tents, allocations, the pair count gauge and the graph
   decrements). Caps with `weight * max(w) < 10^-11` are replaced by their
   fixed upper bound on the whole case. A box passes only if this holds for
   every coefficient, directly or by the correlated-factor test in
   `factor_correlated.cpp`; a box with a positive corner value fails the run.
   Six exact polynomial identities per case (672 in all) check the
   construction of `P`.
7. **Objective.** `U_raw = (F(0) + kappa + c_0 + conditional constants)/8 + 4769/200000000`,
   where `4769/200000000` is the total graph charge. The full composer then
   also proves the radial origin credit `g` and gives
   `D_final = U_raw - g/8 = 0.62999945839124970...`
   `< 78749932299/125000000000`. Later stages replace this credit, so stage 1
   reports `U_raw` as `raw_objective` and `D_final` as `credited_objective`.

Items 1 to 6 are the hypotheses of Theorem 2.1 of the paper with `beta = 0`,
and give `Delta_4 <= U_raw`. The theorem itself (Poisson summation and
periodization) is not checked here.

## Components

[`COMPONENTS.md`](COMPONENTS.md) lists every V7 local-energy component, the
lemma that makes its energy nonnegative, where that lemma is proved
(`proofs/` or a cited theorem), the program that checks it, and whether a run
recomputes the check. It also lists what rests only on written proofs.

## Read set and what is in this folder

`record/` holds exactly the files that the full run opens for reading:
96 files. Two of them (`record/geometry/...`) are read only by
`checks/replay_graphs.py`. They are byte copies of the private source tree at the
same relative paths. I found them by tracing a full run with a
`sys.addaudithook` hook that logged every file open, directory listing and
subprocess launch, including in the forked workers and the compiler calls.
Files that the run writes (outputs, the fresh spherical replays, the compiled
libraries) are excluded. The four V7 data files are read from `../../data/v7/`.
`MANIFEST.json` lists path, size and SHA-256 of every exported file;
`run.py` refuses to run if any file differs or an unlisted file is present, except for non-executable Finder `.DS_Store` metadata.

Most of the exported Python files are hash-bound older coordinators and
helpers: the V7 coordinator pins its whole ancestry by SHA-256 and checks it
at start and end. They are exported because the binding reads them.

Not exported, because the run never reads them: the rest of the 2,957-file
research archive, the historical binaries, the historical execution records
of the 111 V6 branches, and `affine_latitude_certificates.json`. `origin.json`
names that last file only in a `proof_file` provenance string; the affine
latitude proofs are stored inline in `origin.json` and recomputed.

## What is recomputed and what is replayed

| Item | Full mode | Quick mode |
|---|---|---|
| Scalar Fourier identity and Sturm signs | recomputed (coordinator, composer, two independent checkers) | recomputed (composer, quick parent checks, two independent checkers) |
| 628 Gram rays, origin identity | recomputed | recomputed |
| Spherical count certificates | recomputed (`certify.py`) | recomputed (quick parent checks) |
| Signed spherical energies, affine latitude, count polynomials | recomputed | recomputed |
| Graph capacities | recomputed (coordinator, `checks/replay_graphs.py`, independent moderate checker) | recomputed (same) |
| Broad19 rays inside the Gram kernel | recomputed (`checks/check_broad19_rays.py`) | recomputed |
| 7 pair pieces | recomputed | recomputed |
| 112 cyclic branches | recomputed | **not run** |
| Origin credit and `U_raw` | recomputed (composer) | recomputed (composer `--prepare`) |

Nothing is taken from a stored success flag. `ORIGIN_CREDIT_VERIFIED.json`
is an input of the composer: it holds the credit vectors `a` and `b`. The
composer rebuilds `M` from the Gram rays, checks `M b = a`, `M u = v`,
`g = 2 u.a - a.b > 0` and the shifted square exactly, compares every stored
number with its own result, and proves `a.B(r) >= 0` on each radial piece
from nonnegative Bernstein coefficients, an exact Sturm count, or a
factorization it multiplies out and compares coefficient by coefficient.

## Runs on 2026-11-13

All runs used a fresh copy of the whole `D4` folder in a temporary directory,
on macOS arm64 (10 cores, 16 GB, shared with other jobs; load average 16 to 38).
An audit hook traced every file open: no run read anything outside the copy
except the Python installation and `/System/Library/CoreServices/SystemVersion.plist`.

| Run | Exit | Wall time | Peak RSS | Result |
|---|---|---|---|---|
| full, `--workers 4` | 0 | 1665 s | 1.98 GB summed over the process group | `PASS`, 112 branches, 391,420 nodes, 195,766 leaves, 50,934 correlated-factor leaves, 672 identities, 7 pair pieces |
| quick, `--workers 1` | 0 | 82 s | 0.20 GB | `PASS`, 0 branches |

The full run above predates `checks/check_broad19_rays.py` and
`checks/replay_graphs.py`. Both run after the coordinator, identically in
both modes, and passed in the later quick runs.

The full run spends almost all its time in the coordinator (1660 s); the
exterior case (5,5) alone takes about 13 minutes. Peak memory grows with
`--workers`: about 0.2 GB for the parent plus up to about 0.45 GB per worker.
Use `--workers 2` on a machine with less than 4 GB free.

Both runs give the same `U_raw`, equal to `exact_raw_objective` in
`census_floors.json`, and the full run's credited objective equals the
accepted 2026-10-19 value exactly. The branch, node and leaf counts also
match the 2026-10-19 run.

Negative controls, each on its own fresh copy, in quick mode:

| Change | Result |
|---|---|
| `pair_origin_v[0]` in `data/v7/origin.json`: numerator +1, manifest unchanged | `FAIL`, exit 1, at the manifest check |
| same change, manifest hash updated | `FAIL`, exit 1: the coordinator rejects the physical pair row identity |
| weight of ray 0 of angular block 1 in `data/v7/kernel.json` multiplied by 100, manifest hash updated | `FAIL`, exit 1: pair piece 1 `[501473/500000, 1005643/1000000]` gives `VIOLATION` at `r = 501473/500000` |
| `claimed_even_ordered_edge_bound` 1880 changed to 1878 in `record/current/six_exterior_graph_route/moderate_exterior_graph_certificate.json`, manifest hash updated | `FAIL`, exit 1: the coordinator's source hash pin rejects it. Run directly on the changed file, `verify_moderate_exterior_graph.py`, `checks/verify_moderate_independent.py` ("Unexpected capacity") and `checks/replay_graphs.py` each exit 1 |

After `check_broad19_rays.py` and `replay_graphs.py` were added, a quick run
on a fresh copy passed again (exit 0, 70 s, 0.23 GB). Its traced read set is
exactly the 96 files in `record/`.

## Trust base

- CPython 3.14.6, `numpy` 2.5.2 (object and int64 arrays; no decision
  uses floating point), `python-flint` 0.9.0 with FLINT 3.6.0 for exact
  integer and rational arithmetic.
- The C++ helpers `fast_split.cpp` (integer midpoint subdivision) and
  `factor_correlated.cpp` (correlated-factor acceptance), compiled fresh by
  each run with `-std=c++20 -O3 -march=native`. The tested compiler was
  Apple clang 21.0.0 (clang-2100.3.34.2). No prebuilt binary is used.
  Both helpers hash their source and the built library before and after
  the run.
- The reviewed Python sources in `record/` and `checks/`, and the hardware.

This is not a formal proof. The per-leaf subdivision trees are not stored;
rerunning is the way to check them.

## Provenance

Private source: `work/d4/six-graph-v7-bound/source/` (driver `replay.py`,
`record/`). The accepted 2026-10-19 replay used the same coordinator and
composer and gave the same objective. `CHANGES.md` lists every file that
differs from its original.
