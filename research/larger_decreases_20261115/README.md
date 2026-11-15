# Larger decreases, 15 November 2026

Six new certified upper bounds improve the preceding main, commit
`6e7990802b8e045a78e3f54c99669d5fe9b5c437`:

| Dimension | Previous strict upper | New strict upper | Further decrease |
|---|---:|---:|---:|
| 17 | 0.018204096 | **0.018140145562** | **0.35%** |
| 18 | 0.013224473 | **0.013077495192** | **1.11%** |
| 19 | 0.009521000 | **0.0093653907401** | **1.63%** |
| 20 | 0.006936042 | **0.0067425778508** | **2.79%** |
| 21 | 0.004985977000 | **0.0049104292815** | **1.52%** |
| 48 | 4.13385739314e-7 | **4.0843570355e-7** | **1.20%** |

These are packing fractions and unrestricted-packing upper bounds.
Percentages compare with the preceding repository bounds, rather than the
older published comparators in the root table. No exhaustive worldwide
priority assertion is made. `REGISTERED_BOUNDS.json` retains exact inputs.

The successful change applies signed logarithmic radial profiles to
dimensions omitted from the earlier profile campaign: a=-96,b=120 in
17–21D and a=-192,b=216 in 48D. Angular polynomials are reoptimized in the
target dimension. Degree-95 Laguerre corrections enrich the degree-111
radial seed, with root- and pair-guided continuation until the complete
global signs pass. The analytic theorem is unchanged; see the
[shared proof](../../verifications/profile_enrichment_20261114/PROOF.md).

## Verification and exports

Each `verifications/dimN/refinement_larger20261115/` contains frozen
rational inputs and a verifier that executes all four primary and
independent angular/radial proof programs. Both implementations share
FLINT/Arb, with different reconstruction and enclosure formulas. Their
independence is not a separate arithmetic backend or external review.

```sh
python verifications/dim20/refinement_larger20261115/verify.py
python research/larger_decreases_20261115/audit_exports.py 17 18 19 20 21 48
python research/larger_decreases_20261115/test_registration_gates.py 17 18 19 20 21 48
python comparison.py
```

`REGISTERED_REPLAY.json` records the complete registered wrapper runs.
`EXACT_EXPORT_AUDIT.json` reconstructs the exact rational polynomials,
guarded origin costs and objectives from saved binary64 coordinates.
`REGISTRATION_CONTROLS.json` records 42 controls on the real registry/hash
prefix, including altered bounds, witness tampering and historical use.
Source mathematical inputs remain unchanged. Changed historical wrappers
were checked through their actual registry/hash prefixes, and both shared
origin-algebra programs were rerun. The complete historical aggregate,
including the unchanged long 4D base replay, was not rerun for this release.

## 4D: no stronger registered result

The 4D bound remains **0.629520189662**. The following searches did not
produce a substantial certified improvement:

- Reoptimizing degree-47 and degree-95 scalar Fourier functions against
  the frozen kernel gave only microscopic numerical changes.
- Fixed-cardinality three-point energy SDPs were tested at radial/angular
  degrees (4,4), (6,8) and (8,8), including a joint spherical/reallocation
  objective. They did not produce an admissible stronger certificate.
- Radial moment-cone QPs, finite-count Gram projections and conditional
  origin projections did not expose a percent-level reserve in the tested
  relaxations. These were numerical diagnostics, not impossibility proofs.

The saved joint-energy proposal violates the cyclic sign at the valid
spherical triple **(-5/16,1/2,1/2)**. Treating its binary64 coefficients
as exact rationals gives a strictly positive value about **1.08058e-6**.
`dim4_research/EXACT_REJECTION.json` stores the exact fraction and Gram
determinant. Reproduce that rejection with:

```sh
python research/larger_decreases_20261115/audit_4d_rejection.py
```

Solver labels such as `optimal` do not establish feasibility on a continuous
domain. None of these 4D objectives enters the result registry. The current
degree-zero origin is `origin.json`'s separately saturated vector `u`;
the historical first coordinate of the kernel ray is not that origin.

## Numerical discovery

Install `requirements-discovery.txt` in a separate environment. Exact tested
versions are in `RUNTIME.json`. Example:

```sh
python research/larger_decreases_20261115/search_angular.py 20 --a=-96 --b=120 --rounds 8 --tag fresh_n20
python research/larger_decreases_20261115/search_radial.py 20 --s 5/2 --a=-96 --b=120 --degree 95 --cut 60 --angular research/larger_decreases_20261115/fresh_n20/angular.json --tag fresh_n20_radial
python research/larger_decreases_20261115/continue_radial.py research/larger_decreases_20261115/fresh_n20_radial/candidate.json --rounds 28
```

Use a new tag; discovery commands modify their output directories. Solver
versions need not reproduce identical coordinates. Frozen exports and
verification use exact retained data. A raw `certified:false` file is not
a density theorem.

The campaign also tested scale changes, degree-31 angular spaces,
degree-127/159 radial spaces and positive mixtures of profile-kernel
candidates. Those tested controls did not give larger usable gains.
Latest successful inputs remain readable. Earlier rounds and unselected
experiments are retained in `discovery_checkpoints.tar.xz`, with a complete
SHA-256 inventory in `ARCHIVE_MANIFEST.json`. Extract from the repository
root; member paths are relative to that root. Every member was read back
and hashed before its loose copy was removed.
