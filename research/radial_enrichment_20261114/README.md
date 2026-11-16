# Radial enrichment: 14 November 2026

This is the historical 14 November release. Its three certified improvements
are now superseded by the profile certificates selected in the
[current registry](../../RESULTS.json) and [leaderboard](../../RESULTS.md).
The complete primary and independent proof chains remain available for each
frozen witness; discovery failures stay in this research package.

| Dimension | Previous repository strict upper | New strict upper | Decrease |
| --- | --- | --- | --- |
| 45 | 0.000001311423575 | 0.000001287489137 | 1.825073% |
| 46 | 0.0000008777229084 | 0.0000008756454737 | 0.236685% |
| 47 | 0.0000005998306274 | 0.0000005972403006 | 0.431843% |

The old angular kernels remain exactly unchanged. The successful changes are
degree-47 corrections with extended pair coupling in 45D and 47D, and a
degree-63 correction with the same extended coupling in 46D. Minimum squared
separation is 6, kappa=17/8, coupling covers X in [6*pi,40], and physical
negativity is imposed only from 40 onward. Every frozen M receives the exact
guard 1000001/1000000. Exact positive points beyond 7*pi show why the previous
physical-tail sign template does not cover these witnesses.

## Reproduce proofs

Install the root verification requirements and run each new package:

```sh
python verifications/dim45/refinement_radial20261114/verify.py
python verifications/dim46/refinement_radial20261114/verify.py
python verifications/dim47/refinement_radial20261114/verify.py
python research/radial_enrichment_20261114/verify_exports.py
python research/radial_enrichment_20261114/test_registration.py
```

`python run_all.py --jobs 4` includes every new proof and all preserved upper
certificates. Assertions must remain enabled. The new proof packages depend
only on Python and python-flint 0.9.0.

## Bounded discovery and failures

`search.py` is a finite sampled LP, not a verifier. Its span is a nonnegative
multiple of the original normalized degree-111 48D source plus normalized
target-dimensional Laguerre Fourier eigenfunctions. `model.py` constructs
the exact source and basis. Degrees 47, 63 and 47 were selected for 45D, 46D
and 47D respectively. No wider dimension sweep or hidden optimizer run is
claimed. `continue_search.py` adds exact-root-guided points and numerical pair
separators, retaining checkpoints and rejection logs under each `history/`.

Initial sampled witnesses failed on narrow Fourier/physical sign islands and
pair peaks between grid nodes. Root cuts repaired those failures. One 45D
HiGHS attempt returned an unknown/infeasible status; retry with the alternate
LP method recovered the feasible continuation. The first rational exports
used a 100ppm M guard; the final frozen exports use 1ppm and pass complete
interval verification. Earlier checkpoints remain research, even when their
numerical objectives are lower than the certified final objective.

For optional fresh discovery install NumPy, SciPy and mpmath alongside
python-flint; tested versions are recorded in `environment.json`. Example:

```sh
python research/radial_enrichment_20261114/search.py 46 --degree 63
```

`export_check.py` performs rationalization and radial proof gates. Its output
alone is not a density certificate: complete angular/cyclic replay is also
required. `verify_exports.py` deterministically reconstructs the frozen final
certificates and proves the new radial span degrees exactly. The proof
packages read frozen rational inputs rather than trusting discovery floats.

## Comparisons and replay scope

The source commit is 72028bb0c5eaaaeede50cd418ca09339ee404285. The first
(unrestricted-packing) table at https://cohn.mit.edu/sphere-packing/ was
successfully fetched on 14 November 2026. Its 45D, 46D and 47D upper entries
match the preceding repository results to the displayed precision. Raw HTML
and row extracts are preserved under `sources/`. The separate lattice-only
table is not a general-packing comparator. External comparison witnesses are
not independently replayed, and this is not an exhaustive worldwide priority
search.

Historical proof inputs remain unchanged. Registration wrappers accept the
newer package while continuing to verify their own frozen bounds. The full
replay also exposed a pre-existing 16D historical registration assertion that
still required the superseded package as current; its compatibility check was
updated without changing its witness, hashes, mathematical gates or bound.
The historical aggregate exercised all 73 scheduled proof invocations and
reported the old 16D registry assertion as its failure. That corrected chain
passed a separate complete replay, as did the final shared origin and 48D
rounding checks. A clean aggregate rerun from scratch is not claimed.
`REPLAY_STATUS.json` and the retained logs distinguish this scope from the
successful complete new-certificate replays. Fresh generated historical
reports are archived with an inventory, while their previously committed
versions remain preserved in their original locations.
