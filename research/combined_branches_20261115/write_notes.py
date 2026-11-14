"""Generate the dated integration notes from the actual registered results."""

from pathlib import Path
from decimal import Decimal, localcontext
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
r = json.loads((HERE / "REGISTERED_BOUNDS.json").read_text())
rows = sorted(r["bounds"], key=lambda x: x["dimension"])
dims = [x["dimension"] for x in rows]
assert dims == [23, 26, 27, 28, 29, 30, 31, 33, 34]


def percent(row):
    with localcontext() as ctx:
        ctx.prec = 60
        return f"{100*(1-Decimal(row['strict_upper'])/Decimal(row['previous_strict_upper'])):.2f}%"


table = "| Dimension | Previous main upper bound | New strict upper bound | Further decrease |\n|---|---:|---:|---:|\n"
for row in rows:
    table += f"| {row['dimension']} | {row['previous_strict_upper']} | **<{row['strict_upper']}** | {percent(row)} |\n"
changes = f"""# Nine further certified upper bounds — 15 November 2026

The two requested research branches and the lower-bound follow-up are merged
into `main`. The radial branch was already an ancestor of the preceding main;
the density branch and lower-bound work were merged in `0fb0b9b` and `ca3fc0d`.

Applying the merged logarithmic-profile tools to previously untested dimensions
produced **nine new certified upper-bound refinements in 23D, 26–31D and 33–34D**.
All values below are packing fractions, not center densities. The percentage
column compares with the preceding main, not with a separately replayed paper
or an exhaustive publication-priority survey.

{table}
Each frozen certificate has `a=-96,b=120`, minimum squared distance `s=n/8`,
a newly optimized degree-19 angular polynomial and degree-95 Laguerre
enrichment of the retained degree-111 radial seed. The framework and analytic
proof are the [merged profile family](verifications/profile_enrichment_20261114/PROOF.md).
The improvements come from new witnesses in dimensions omitted by the preceding
campaign. No external neighbor-count bounds enter these nine certificates.

All nine packages passed the complete primary angular/cyclic, primary radial,
independent angular/cyclic and independent radial proof programs, first as
frozen witnesses and then through their registered wrappers. The implementations
share FLINT/Arb; this is not a separate arithmetic backend. The source-hash and
registry controls also reject corrupted witnesses and falsely tighter current
registry values. Historical packages remain replayable under a stronger
registered successor without changing their mathematical inputs.

Exact exports were reconstructed from the retained binary64 discovery
coordinates. Complete continuation checkpoints are retained in verified
archives; the selected latest witnesses remain directly readable.

No additional lower-bound world record is asserted. The merged density branch
contains a certified 65D construction and higher-dimensional periodic packings,
with construction-specific comparisons and unestablished record status. The
combined lower-bound search and its caps are documented in the
[integration research notes](research/combined_branches_20261115/README.md).

To replay the nine registered certificates:

```bash
OPENBLAS_NUM_THREADS=1 python research/combined_branches_20261115/replay_registered.py {' '.join(map(str,dims))} --jobs 2
python comparison.py
```

The complete merged-baseline replay passed all 89 proof invocations and the
shared algebra/rounding checks. The nine new registered packages were then
verified separately. `BASELINE_REPLAY.json` records that exact scope; the
updated `python run_all.py` also includes the nine new packages.
"""
(ROOT / "CHANGES_20261115.md").write_text(changes)
research = f"""# Combined branch results, 15 November 2026

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
OPENBLAS_NUM_THREADS=1 python research/combined_branches_20261115/replay_registered.py {' '.join(map(str,dims))} --jobs 2
python research/combined_branches_20261115/test_registration_gates.py {' '.join(map(str,dims))}
python research/combined_branches_20261115/audit_exports.py {' '.join(map(str,dims))}
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
"""
(HERE / "README.md").write_text(research)
readme = ROOT / "README.md"
text = readme.read_text()
start = "<!-- combined-20261115:start -->"
end = "<!-- combined-20261115:end -->"
note = f"""{start}
> [!IMPORTANT]
> **Nine further upper-bound refinements — 15 November 2026.** The merged profile tools improve **23D, 26D–31D and 33D–34D**, with further reductions up to **2.62% in 29D** and **2.60% in 30D**. All nine registered packages pass complete primary and independent proof programs. [Exact bounds and before/after table](CHANGES_20261115.md) · [Merged research and replay scope](research/combined_branches_20261115/README.md).
{end}

"""
if start in text:
    a = text.index(start)
    b = text.index(end, a) + len(end)
    text = text[:a] + note.rstrip() + text[b:]
else:
    text = text.replace("> [!IMPORTANT]", note + "> [!IMPORTANT]", 1)
marker = "Additional merged constructions include an explicit"
if marker not in text:
    anchor = "## Registered results"
    paragraph = "Additional merged constructions include an explicit **65D minimum-12 lattice** and periodic packings in 60–62D and 66–71D. Their comparison scope and unestablished world-record status are documented in [the certified higher-dimensional research](research/novel_routes_20261113/README.md).\n\n"
    text = text.replace(anchor, paragraph + anchor, 1)
readme.write_text(text)
print("Wrote release and research notes for", len(rows), "certified upper refinements")
