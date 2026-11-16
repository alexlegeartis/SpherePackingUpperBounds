# Research index

The current lower-bound leaderboard is the [packing constructions table](../README.md#packing-constructions-lower-bounds).
The current upper-bound leaderboard is [RESULTS.md](../RESULTS.md), generated
from [RESULTS.json](../RESULTS.json). Each entry links to its certificate,
normalization, dependencies and comparison sources.

A successful experiment enters a leaderboard only after a complete certificate
replay and a documented comparison with the registered result. An improvement
over a named comparison is labeled with that comparison; a best-known or SOTA
claim additionally needs a dated literature review. Numerical objectives,
capped searches and unembedded abstract Gram matrices remain research.

Discovery histories stay here even when a result succeeds. Proof inputs already
linked from a leaderboard can also remain here: moving them would break their
replay and provenance paths. Historical certificates remain available after a
stronger result replaces their leaderboard entry.

## Reconciled branches

| Workstream | Outcome and current location |
| --- | --- |
| [Density breakthroughs, 13 November](density_breakthroughs_20261113/README.md) | Already merged. Exact finite certificates, rejected candidates and local-optimality results; the 36D lower bound is certified in [constructions/dim36](../constructions/dim36/). |
| [Higher-dimensional constructions](novel_routes_20261113/README.md) | Already merged with the density branch. Certified lower bounds are linked from the lower-bound leaderboard; later improvements take precedence. |
| [Radial enrichment, 14 November](radial_enrichment_20261114/README.md) | Already merged. Historical certified improvements in 45D–47D, now superseded by the profile certificates selected in the registry. Discovery failures and proof exports remain available. |
| [Non-minimum-basis sections, 16 November](nonminimal_sections_20261116/README.md) | Unsuccessful randomized screen of 50,000 abstract rank-five forms. All three shortlisted forms fail exact minimum checks. No packing or leaderboard entry. |

The [branch reconciliation record](branch_integration_20261116/README.md)
pins all three branch heads and records fresh replay scope and leaderboard checks.

## Related merged investigations

| Workstream | Status |
| --- | --- |
| [Lower-bound search and follow-ups](lower_bound_search_20261114/README.txt) | Exact obstructions and unsuccessful low-dimensional searches, with explicit limits on what was excluded. |
| [Combined branch tools](combined_branches_20261115/README.md) | Certified upper-bound improvements plus lower-bound searches. The registry identifies the current certificates. |
| [Adaptive layers](adaptive_layers_20261115/README.md) | Certified constructions and finite research screens; see the lower-bound leaderboard for current values. |
| [48D–96D audit](range48_96_20261115/README.md) | Exact constructions, source comparisons and dimensions where no improvement was obtained. |
| [CM-plane construction and frontier](cm58_20261116/FRONTIER.md) | Certified advances and rejected or unresolved targets. |
| [Correlated normals and Craig/code routes](simplex_normals_20261116/PROOF.md) | Certified 86D–89D improvements and a separately scoped unsuccessful code search. |
| [Nine remaining lower-bound dimensions](remaining_gaps_20261116/PROOF.md) | Exact Ne64 local-optimality certificate; exhaustive rejection of a bounded 64D form family; a 92D code symmetry obstruction. No new packing record. |
| [Tilted 55D lamination and source correction](tilted55_20261116/PROOF.md) | New certified minimum10 lattice; explicit norm6 counterexample withdraws the former55D claim. Improves the repaired lamination; no world-record claim. |

Use each package's own reproduction commands. `run_all.py` replays the registered
upper-bound proof chains; lower-bound and exploratory checks have separate runners.
