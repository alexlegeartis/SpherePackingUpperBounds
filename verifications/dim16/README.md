# Dimension 16: log-hinge cyclic kernel plus a widened spherical count

The current registered refinement is **0.024900975**, with the complete portable proof in [widened_reconstruction/](widened_reconstruction/). Run `python widened_reconstruction/verify_all.py`.

The original certificate preserved below has strict upper target **0.024920988**.

The exact density expression begins 0.0249209879146585408. It is a small
improvement over the pinned published 0.02492121: about 0.000890807% of
that upper bound, or 0.00217369% of its gap to the achieved Barnes-Wall
density pi^8/645120. No claim of a large gap closure is made.

## Full portable replay

Run `python verify_all.py` with NumPy and python-flint installed. The
script explicitly rejects Python -O/-OO, pins every mathematical input,
checks that the radial and geometry components use identical angular and
profile data, replays all three proof components, verifies the exact
cap-to-shell conversion, and recomputes the final density with Arb.
There is no optimizer, network, external repository or undistributed
mathematical input in this replay. `RESULT.json` records its outcome.

The components are:

- `count/`: three preserved sibling directories containing the widened
  spherical-code certificate, exact pair/Sturm audit, exact cyclic
  polynomial reconstruction and two independent integer-only Bernstein
  traversal replays. They certify cardinality at most7802 for inner
  products at most5003/10000.
- `geometry/`: exact angular PSD/cap/sign-transition checks, the complete
  log-hinge cyclic rectangle covering, the all-scale argument and the
  infinite tail. The count and radial function are not assumed here.
- `radial/`: exact Fourier reconstruction, both global sign rays,
  complete piecewise pair coverage, one deterministic repair, exact
  normalization and rounding, and negative controls. Its preserved
  component metadata calls the count an external assumption because
  that component alone does not replay it; this top-level entrypoint
  explicitly supplies and verifies that assumption.

Run `python verify_independent.py` for the complete separate independent
audit. Its three directories are preserved alongside the count directories
under `count/`, sharing the same immutable count inputs. It uses an
independent geometry derivative calculation, Fourier implementation,
rational Sturm rays, pair covering and rounding, plus ten negative controls.
No external directory is needed. Its outcome is `INDEPENDENT_RESULT.json`.

See `PROOF.md` for the theorem linking the components and `BENCHMARKS.md`
for the comparison sources. The independent audit is not replaced by the
primary replay.

## Discovery is not certification

The initial radial polynomial had Fourier/physical sign defects. After
exact-root cuts it still had a tiny positive pair value, which is preserved
as a negative control. Only the final explicitly repaired rational function
passes all radial gates. The contact spherical witness was not simply
declared valid at a wider angle: its blocks were reweighted and both new
global inequalities verified directly. The original enormous contact SOS
proof is not a hidden dependency.

The separate richer factor-row numerical experiment is not part of this
bound. All three of its promising-looking candidates were rejected at
exact feasible cyclic counterexamples. They must not replace these frozen
inputs or be entered as records.

## Discovery provenance

The reproducible research companions, relative to the repository root,
are listed in `DISCOVERY.md`. Replaying proof inputs requires none of the
discovery solvers. Optimizer timings and floating solutions are historical
observations, not promises of bitwise rediscovery on a new solver version.
