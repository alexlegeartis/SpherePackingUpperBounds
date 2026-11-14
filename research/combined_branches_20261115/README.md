# Combined branch results, 15 November 2026

The integration combines `research/radial-enrichment-20261114`,
`research/density-breakthroughs-20261113`, and the local lower-bound follow-up.
The first was already in main. Both other heads are ancestors of merge
`ca3fc0d`; the research branches were retained.

The principal new results are nine certified **upper** bounds. No new lower-bound
world record was established. Exact values and percentage comparisons are in
[the change record](../../CHANGES_20261115.md) and `REGISTERED_BOUNDS.json`.

## Certificate checks

Every result has a frozen package at
`verifications/dimN/refinement_combined20261115/`. `REGISTERED_REPLAY.json`
records complete registered wrapper runs, each executing all four proof paths.
`REGISTRATION_CONTROLS.json` records 63 controls on the actual wrapper prefix.
Current packages must match their own claimed bound exactly; historical
packages can accept a separately registered stronger successor. Neither
condition substitutes for the mathematical proof programs.

`EXACT_EXPORT_AUDIT.json` checks reconstruction of both rational polynomials,
the origin guard and objective, and the angular data from retained discovery
coordinates. The shared analytic theorem and proof are inherited from
`verifications/profile_enrichment_20261114/PROOF.md`.

```bash
OPENBLAS_NUM_THREADS=1 python research/combined_branches_20261115/replay_registered.py 23 26 27 28 29 30 31 33 34 --jobs 2
python research/combined_branches_20261115/test_registration_gates.py 23 26 27 28 29 30 31 33 34
python research/combined_branches_20261115/audit_exports.py 23 26 27 28 29 30 31 33 34
```

`BASELINE_REPLAY.json` records the successful clean merged-baseline aggregate:
89 proof invocations, with unchanged mathematical source files. Its 212 fresh
generated reports are retained in `baseline_replay_outputs.tar.xz`; the
complete output is in `baseline_upper_replay.log`. The nine new registered
packages were replayed separately, so this is not represented as a single
aggregate invocation at the final revised tree.

The compact incoming research checks are
`research/density_breakthroughs_20261113/verify_exact.py`,
`research/novel_routes_20261113/verify_checkpoint.py`, and
`research/lower_bound_search_20261114/broad_dimensions_followup/verify_all.py`.
The higher-dimensional checker rechecks the saved four-shard Ne64 enumeration
outputs; it does not repeat the 4.88-billion-node traversal. Parent minima remain
published inputs, exactly as documented on that branch.

## Combined lower-bound searches

The exact dual-cluster method from our lower-bound follow-up was applied to
the Ne64 and Gamma72 parents from the density branch. For physical cluster
diameter at most one it checks minimum-basis normal sections in 61D, 62D, 68D
and 69D. All four searches match the existing branch constructions; no graph
is left node-capped. The normal Gram denominator is three in the Ne64 cases,
and the score includes the parent covolume. `high_parent_cluster_results.json`
records the exact rational comparisons.

For the fixed Gamma72 root sections of ranks 6, 7 and 8, the exact maximum
cluster size at diameter 1/2 is nine. These checks cover every eligible dual
point and confirm the branch's nine-point clusters at this diameter. They are
not global packing bounds or bounds on other section shapes/diameters.

An order-two Ne64 hole of squared depth at least nine would give a 65D
construction of center density `3^16`, improving the branch's
`3^17/sqrt(10)`. The floating discovery screen tested 100,000 parity classes,
promoting 98,018 after a 1,000-node cap. The 100 most promising numerical
upper-distance candidates were refined to a 2,000,000-node cap. Exact witnesses
reject twelve; eighty-eight remain inconclusive. No distance lower certificate
or new packing bound follows from those caps. `ne64_order2_screen.npz` retains
the explicit residues and witnesses, and `verify_lower_checks.py` audits them.

```bash
OPENBLAS_NUM_THREADS=1 python research/combined_branches_20261115/verify_lower_checks.py
OPENBLAS_NUM_THREADS=1 python research/combined_branches_20261115/search_high_parent_clusters.py
```

## Discovery provenance

The selected angular and radial candidates remain in the
`research/ambitious_profiles_20261114/integrated_nN_*` folders. `extend_profiles.py`
calls the merged numerical tools, performs the full discovery gates, freezes
successful witnesses and runs the four direct proof programs. Those numerical
tools alone cannot register a result; `register_bounds.py` requires successful
complete frozen replays and matching source hashes.

Continuation histories are in `discovery_history.tar.xz` and
`discovery23_history.tar.xz`. Both archives were read back member by member and
SHA256-checked before the loose historical copies were removed. Their manifests
record all members and retained latest files. To restore, extract from the
repository root; the archived paths are relative to that root. Some retained
raw candidates say `certified:false`; the density claim comes from the complete
frozen proof packages, not from the solver status.

The earlier general self-orthogonal minimum-12 construction remains open.
