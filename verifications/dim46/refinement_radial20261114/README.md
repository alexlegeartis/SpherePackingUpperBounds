# Dimension 46: radial continuation of 14 November 2026

The exact witness proves **density <0.0000008756454737**, improving the previous
repository bound 0.0000008777229084 by **0.236685%**. This is an unrestricted
packing-density upper bound. The first table of [Cohn's live sphere-packing
page](https://cohn.mit.edu/sphere-packing/), retrieved on 14 November 2026,
still displays the preceding upper value. The preserved source snapshot and
the distinction from the separate lattice-only table are recorded in
[`research/radial_enrichment_20261114/sources`](../../../research/radial_enrichment_20261114/sources/).
This comparison does not assert an exhaustive worldwide priority search.

Run from a checkout with Python 3.11+ and `python-flint==0.9.0`:

```sh
python verifications/dim46/refinement_radial20261114/verify.py
```

Never use Python `-O`. No solver or network is needed for proof replay.
The entrypoint checks frozen input hashes and the current registry, then runs
all five primary and independent angular, cyclic and radial gates. Actual
completed replay status is in `full_replay_result.json`.

The angular polynomial is unchanged from `verifications/dim46/refinement_extended_pair`.
The new radial function is a nonnegative multiple of the original degree-111
48D seed plus a degree-63 correction in normalized dimension-46
Laguerre Fourier eigenfunctions. Pair coupling is proved through X=40, and
physical negativity only from X=40 onward. The exact M guard is
1000001/1000000. Read [PROOF.md](PROOF.md) for the global argument and
[independent/AUDIT.md](independent/AUDIT.md) for verification scope.

All preceding witnesses remain available. The discovery history and
deterministic exact-export checks are separate from these proof inputs.
