# Research branch reconciliation — 16 November 2026

The three branch heads were compared with main at
[`5fe8f91`](https://github.com/alexlegeartis/SpherePackingUpperBounds/commit/5fe8f91b15bdb1952f53bc9f277cb2ddf3b81482).
The density and radial branches were already ancestors of main. Only the
non-minimum-basis section branch had unmerged content, consisting of its
research README and search script. Its full history is incorporated by a
merge; the removed exploratory Actions workflow is not restored.

| Branch | Pinned head | Disposition |
| --- | --- | --- |
| `research/density-breakthroughs-20261113` | [`1d073d1`](https://github.com/alexlegeartis/SpherePackingUpperBounds/commit/1d073d1cffdef2c3f62533f6cec8595e3abb5c64) | Already merged by `0fb0b9b`; retain certificates and their later leaderboard successors. |
| `research/radial-enrichment-20261114` | [`8a080d6`](https://github.com/alexlegeartis/SpherePackingUpperBounds/commit/8a080d6a00909391a91ecbac097e4004c404c40d) | Already incorporated before `1871245`; retain the three historical certificates and the stronger profile entries. |
| `research/nonminimal-sections-20261116` | [`86f5a02`](https://github.com/alexlegeartis/SpherePackingUpperBounds/commit/86f5a02d8a3cc241c7443151c55627d85214628e) | Merge the unsuccessful search under `research/nonminimal_sections_20261116/`; no new density claim. |

The pinned commits preserve branch identity after merged branch heads are
retired. The [research index](../README.md) explains promotion and retention.

## Leaderboard reconciliation

All 43 upper-bound entries on each already merged branch are retained or
improved in the current registry. No historical value replaces a stronger
current one. In particular:

| Dimension | Radial branch's certified upper | Current registered upper | Current certificate |
| --- | --- | --- | --- |
| 45 | 0.000001287489137 | 0.0000012656297919 | [45D profile](../../verifications/dim45/refinement_profile20261114/) |
| 46 | 0.0000008756454737 | 0.00000087124499791 | [46D profile](../../verifications/dim46/refinement_profile20261114/) |
| 47 | 0.0000005972403006 | 0.00000059376953797 | [47D profile](../../verifications/dim47/refinement_profile20261114/) |

The density branch's certified constructions are already represented in the
lower-bound leaderboard, including later stronger successors in its dimensions.
Its exact local-optimality and rejection results remain research, with the
same inherited parent-minimum assumptions. The nonminimal branch adds no
lower or upper bound. Accordingly this reconciliation changes no bound values.

## Fresh checks

The [audit report](audit.json) records pinned heads, registry comparisons,
leaderboard coverage and the checks run during integration:

- `comparison.py`: exact percentages, ordering and generated leaderboard consistency.
- Radial export reconstruction and registration controls, including bad package,
  weaker bound, changed hashes and disabled-assertion rejection.
- All three historical radial wrappers: complete primary and independent proofs.
- Density research: all eight compact exact checks.
- Higher-dimensional checkpoint: all 18 compact checks. This checks the archived
  complete four-shard Ne64 outputs; it does not repeat their 4.88-billion-node traversal.
- Nonminimal search: all 50,000 fixed-seed trials. The three shortlisted forms
  again have exact minima 4, 2 and 2. [Fresh replay log](nonminimal-replay.log).

The old radial registration test treated a smaller successor bound as invalid.
Historical wrappers intentionally accept registered stronger successors. The
test now uses a bound twice the frozen certificate's bound to exercise rejection;
package/hash checks and mathematical proof programs are unchanged.

No full `run_all.py` aggregate or fresh literature-priority survey was performed
for this organizational merge. The current proof inputs and numerical leaderboard
entries are unchanged; the checks above state the actual fresh replay scope.
