# Verification of the D4 bound

This folder holds every program needed to check that every sphere packing in ℝ⁴ has density below 0.629521681653. No mathematical certificate input or project source is read from outside the D4 folder; standard software dependencies are installed separately. The programs read the certificate data in [`../data`](../data) and the code and certificates stored here.

## Run everything

From the D4 folder:

    pip install -r verify/requirements.txt
    python3 verify/verify_all.py --out /tmp/d4-check --mode full

The driver first verifies every file listed by the stage manifests, then runs the five certificate stages in order and passes each stage's fresh output to the next. It also runs the independent kernel reconstruction and the exact stage-1/stage-5 formula link. It then recomputes the final ledger from those outputs in exact rational arithmetic. Finally it checks

    U_raw - a_new/8 - 3/8000000 < 629521681653/10^12.

In full mode, exit 0 and the last line `PASS: every sphere packing in R^4 has density < 0.629521681653 (full mode)` mean every check passed. `--out` must be a new folder. `SUMMARY.json` there holds the exact values, and `logs/` holds one log per stage.

| Mode | What it skips | Time on an Apple M2 Pro | Peak memory |
|---|---|---|---|
| `--mode full` | nothing | base five stages about 31 minutes; identity checks add about 10 minutes | about 1 GB |
| `--mode quick` | the 112 triangle subdivisions of stage 1; stage 2 then takes the committed stage 1 values | base stages about 9 minutes, plus identity checks | about 0.3 GB |

`--workers 4` speeds up stages 1 and 5 and needs about 2 GB. The integrated independent comparison uses one worker. The Lean formalization is unreleased and is not a prerequisite for this replay; see [`../FORMALIZATION_STATUS.md`](../FORMALIZATION_STATUS.md).

**Quick mode does not freshly verify the global theorem.** Its final message and `SUMMARY.json` explicitly mark the omitted sign replay. Do not use `-O` or `PYTHONOPTIMIZE`; inherited checks use assertions.

You need Python 3.12 or later, numpy and python-flint (versions in `requirements.txt`), and a C++20 compiler (`clang++` or `g++`) for stage 1. The stage READMEs give the versions used for the recorded runs. NumPy 2.5.2 is pinned, so older Python interpreters are unsupported.

## The stages

Each stage folder has `run.py --out DIR`, writes `DIR/RESULT.json` and exits 0 only on `PASS`. Each README states the exact theorem the stage certifies, what it recomputes, and its trust base. `CHANGES.md` lists every edit made to the original programs; the inherited exports preserve their mathematical checks. The new kernel comparisons and top-level acceptance checks are described separately below. `MANIFEST.json` or `EXPORT_MANIFEST.json` lists the files with their SHA-256.

| Stage | Folder | What it certifies | Paper |
|---|---|---|---|
| 1 | [`stage1_v7`](stage1_v7) | Signs of the scalar function, the 628 Gram rays (positive semidefinite), the pair inequality on 7 radial pieces, the triangle inequality on all 112 cyclic branches (interval subdivision, fresh C++ build), the local-energy components, and the raw objective U_raw | Section 3 |
| 2 | [`stage2_origin`](stage2_origin) | The radial and angular origin credit: all 165503 count vectors, the spherical count certificates, the 12 radial planes, 3 angular proofs and the 124 group minima | Section 4 |
| 3 | [`stage3_census`](stage3_census) | The finite-count refinement of the 124 group floors; output equals `data/census_floors.json` | Section 5 |
| 4 | [`stage4_graph_slack`](stage4_graph_slack) | The broad-graph capacity 1208 and the density saving 3/8000000 | Section 5 |
| 5 | [`stage5_kernelbox`](stage5_kernelbox) | The far-pair lemma, the kernel-box bounds (every leaf of both split trees re-evaluated with Arb), and the final ledger | Section 6 |

The passage from these finite statements to the density of an arbitrary packing is Theorem 2.1 of the paper. Its periodic-patch proof is written in `paper.tex` and `stage5_kernelbox/proofs/conditional_theorem.md`; it is part of the non-Lean proof and does not depend on Lean. Stage 1 independently checks the scalar signs and exact certificate data used by that theorem.

## Kernel identities

The driver additionally requires:

- [`kernel_independent`](kernel_independent/README.md): a separate FLINT implementation rebuilds the Gram/count/spherical polynomial parts, all seven pair pieces, and the origin/objective. Exact coefficient comparisons cover all 84 full and 28 exterior branches. Nonpolynomial component identities are written out in [COMPONENT_IDENTITIES.md](kernel_independent/COMPONENT_IDENTITIES.md), with machine-checked parameter and boundary contracts. Samples are additional diagnostics, not identity proofs.
- [`kernel_link`](kernel_link/README.md): exact symbolic execution and canonical function atoms link the actual stage-5 formula to stage 1 on all 25 ordered radial cells of the far-pair domain, including all 25 singular angular endpoint faces. The certified lower envelope at the endpoint is explicitly distinguished from the exact kernel formula. All deliberately corrupted identity controls must fail.

The independent code is separate arithmetic/assembly code, but its author inspected stage 1 to resolve conventions; it is not a clean-room implementation. Neither identity check substitutes for stage 1's sign subdivisions or stage 5's interval replay.

## What is cited, proved in writing, or trusted

[PROOF_DEPENDENCIES.md](../PROOF_DEPENDENCIES.md) records the dependency map for the available proof. The finite certificate obligations are discharged by the five stage drivers, the independent kernel reconstruction, and the exact kernel link; the arbitrary-packing transfer is the written theorem cited above.

The graph proof uses the exact certificate in `../data/projected_eight.json`.
Stage 4 reconstructs its Legendre polynomial and proves the sign on the entire
interval `[-1,223/735]` with rational Bernstein subdivision, then derives
`A(3,223/735) <= 8` from the objective. The classical nine-point references
remain attribution and an alternative written route. The component identities
and allocation arguments used by the proof are documented in the stage proof
files and checked by the stage drivers.

The computational trust base includes Python integer/rational arithmetic, FLINT exact polynomials, Arb outward-rounded balls, and the C++20 subdivision engine compiled from source during each run. Symbolic adapters translate the documented nonpolynomial functions; they are not a formal semantics proof of the Python code. The stage-5 interval lower-envelope argument remains in [kernel_box_method.md](stage5_kernelbox/proofs/kernel_box_method.md).

## Failure controls

```sh
python3 -m unittest discover -s verify -p test_driver.py -v
python3 verify/kernel_independent/negative_controls.py --out /tmp/d4-independent-controls
python3 verify/stage5_kernelbox/negative_controls.py --work /tmp/d4-box-controls
```

The kernel-link controls run inside its normal command. The independent controls reject incomplete coverage, perturbed coefficients and malformed exact inputs. The driver tests reject fractional/float counts, duplicate or missing groups, stale manifests, inexact floors, and optimized Python. Failure controls test rejection behavior; they do not establish a density bound.

[The release review](../REVIEW.md) records the available replay, remaining
trust assumptions, and the boundary between exact checks and written proofs.
Hash manifests identify exported artifacts; hashes alone are not mathematical
verification.
