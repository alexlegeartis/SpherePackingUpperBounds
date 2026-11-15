# Wider graph A: exact boundary diagnosis and a separate fallback

## 1. The narrow-A candidate has an exact violation

The first A/B/D global run certified 108 continuous branches and returned
`UNKNOWN` on exterior branch `(3,5)` after 6,304 nodes at depth 66. The
coordinator reported no violating corner. A targeted exact diagnostic
subsequently identified a genuine uncovered positive region.

Take the strict radial point

\[
r=\frac{1007580001}{10^9}=1.007580001,\qquad
s=\frac{13648031841661479}{10^{16}}=1.3648031841661479,
\qquad t=-\frac{751}{1000}.
\]

The exterior coordinate used by the verifier is defined exactly by

\[
z=\frac{2rs(1+t)}{(r+s)^2-49/16}.
\]

This is a valid exterior point with both radial coordinates strictly
inside their assigned pieces. The exact direct evaluator, including every
stored base energy and the complete pair-count gauge, proves

\[
\Gamma_0(r,s,t)>
\frac{179787}{10^{12}}>0.
\]

Its numerical display is approximately \(1.7978722486845845\times10^{-7}\).
Narrow A gives no decrement because \(t<-3/4\). Graph B gives none
because \(s<3/2\), and D gives none because both vertices are not outer.
Therefore this is an exact rejection of the narrow-A/B/D candidate,
independent of the earlier resource-limited `UNKNOWN` result.

The following values use the same strict \(r,s\); all are approximations
to rational enclosures in the retained exact diagnostic.

| Angular inner product | Complete base cyclic value | Narrow-A decrement | Wider-A decrement |
|---:|---:|---:|---:|
| \(-0.755\) | \(-8.7247104280\times10^{-8}\) | 0 | \(3\times10^{-7}\) |
| \(-0.751\) | \(+1.7978722487\times10^{-7}\) | 0 | \(3\times10^{-7}\) |
| \(-0.750001\) | \(+2.1808942118\times10^{-7}\) | 0 | \(3\times10^{-7}\) |
| \(-0.75\) | \(+2.1812236220\times10^{-7}\) | \(3\times10^{-7}\) | \(3\times10^{-7}\) |
| \(-0.7468\) | \(+2.7057613510\times10^{-7}\) | \(3\times10^{-7}\) | \(3\times10^{-7}\) |

At the cutoff itself the original graph switches on and the corrected
value is negative. Immediately below it, the positive base expression is
uncovered. The 19-point angular diagnostic covers selected values from
\(-0.8\) through \(-0.7\); every one is strictly negative after the
wider correction.

A separate 42-point radial diagnostic at \(t=-3/4\) found positive base
values at sampled middle radii 1.36, 1.3648031841661479, and 1.37 for near
radii 1.007580001, 1.007581, and 1.007585. At near radius 1.0076, the
sampled point with middle radius 1.3648031841661479 is still positive,
approximately \(1.5234768\times10^{-9}\). All sampled middle radii are
negative at near radii 1.0078 and 1.008. Every radial diagnostic point is
negative after the wider correction. These are point observations, not
a uniform assertion about the intervening radial or angular intervals.

Evidence and reproduction:

- [wide_a_angle_diagnostic.summary.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_a_angle_diagnostic.summary.json)
- [wide_a_angle_diagnostic.exact.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_a_angle_diagnostic.exact.json)
- [wide_a_radial_diagnostic.summary.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_a_radial_diagnostic.summary.json)
- [wide_a_radial_diagnostic.exact.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_a_radial_diagnostic.exact.json)
- [probe_wide_graph_boundary.py](sandbox:/workspace/scratch/93cb6ce085de/geometry/probe_wide_graph_boundary.py)

## 2. The wider-A local theorem

Retain the near class \([1,103/100]\) and the middle class
\((103/100,3/2]\). Widen the selected cross-class angular band to

\[
[-4/5,-1/2].
\]

The ordered edge count satisfies the exact inequality

\[
e_{A_{\rm wide}}\le862.
\]

The proof has the same matrix harmonic structure as the original A
theorem. All nonconstant harmonic coefficient matrices are explicit
nonnegative rational sums of rank-one matrices. Exact rational Bernstein
subdivision verifies every complete angular sign condition. The resulting
quadratic in the two class counts has its maximum over all nonnegative
integer pairs at \((21,60)\), with exact value

\[
\frac{13795458843363213}{16000000000000}
\approx862.2161777102008<864.
\]

The ordered count is even, giving 862. The count maximization is
unrestricted; it assumes no independent upper bound on either class
size. The minimum selected squared distance remains
\(30909/10000>49/16\), since the band's upper endpoint is unchanged.
Thus the wider graph still has zero full-triangle correction, zero
linear and diagonal terms, and an exterior decrement equal to its
weight when selected.

The new combined certificate retains graph B exactly, with capacity
1102. Both proofs freshly passed through the unchanged exact verifier.
The combined data are in
[exterior_wide_ab_graph_certificates.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/exterior_wide_ab_graph_certificates.json),
with proof report
[exterior_wide_ab_graph_verification.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/exterior_wide_ab_graph_verification.json).

At the same A weight \(3/10^7\), replacing capacity 820 by 862 costs

\[
\frac{(862-820)(3/10^7)}8
=\frac{63}{40000000}=0.000001575.
\]

Including unchanged B and D gives total graph density cost

\[
\frac{862(3/10^7)+1102(1/10^7)+6252(1/500000000)}8
=\frac{47663}{10^9}=0.000047663.
\]

For the unchanged robust base candidate, the corresponding objective
therefore satisfies the exact comparison

\[
J_{\rm base}+\frac{47663}{10^9}
<\frac{63010882009}{10^{11}}=0.63010882009.
\]

This is an objective budget and an exactly proved local repair. A global
bound at this value requires a complete successful continuous verifier.

## 3. Integration files and required coordinator changes

The separate helper
[wide_ab_graph_repair_bounds.py](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_ab_graph_repair_bounds.py)
has the same API as the original helper. Its three exact substitutions
are recorded in
[wide_ab_construction.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_ab_construction.json):

1. Change `CAPACITIES` from `(820,1102,9630)` to `(862,1102,9630)`.
2. Change A's angular lower endpoint from `-3/4` to `-4/5`.
3. Set the default weights to the current exterior-only values
   `(3/10000000,1/10000000,0)`.

The unused third slot is retained solely for compatibility; the current
coordinator requires its weight to be zero. All interval enclosures,
strict middle-class endpoint handling, and correction signs are otherwise
unchanged. Graph B and graph D sources are byte-identical to the earlier
ones.

The directory `geometry/wide_ab_graph_route` contains a ready-to-copy
bundle with all six canonical filenames expected by the integration:

- `exterior_graph_certificates.json`: the combined wider-A/B data.
- `graph_repair_bounds.py`: the compatible wider-A/B helper.
- `verify_exterior_graphs.py`: the unchanged exact matrix graph checker.
- `outer_exterior_graph_certificate.json`: unchanged D data.
- `outer_exterior_graph_bounds.py`: unchanged D interval helper.
- `verify_outer_exterior_graph.py`: unchanged D checker.

Fresh imports from that bundle were checked: D's proof utilities bind to
the freshly loaded graph verifier, and its interval utilities bind to the
fresh wider helper. All three local capacities were freshly reproved.
The result and source hashes are recorded in
[wide_ab_integration_ready.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_ab_integration_ready.json).

A proposed candidate metadata file is supplied as
[wide_ab_graph_energies.json](sandbox:/workspace/scratch/93cb6ce085de/geometry/wide_ab_graph_energies.json).
It pins the six bundled sources and the unchanged robust base data,
retains all three weights, and uses the conditional threshold
`63010882009/100000000000`.

To integrate in a separate coordinator copy, point `graph_dir` and every
corresponding source-contract path at the separate graph directory.
Update its hardcoded A/B capacity assertion to `[862,1102]`, its helper
capacity assertion to `(862,1102,9630)`, and A's required angular band to
`(-4/5,-1/2)`. Retain the canonical module names and the existing fresh
import checks so the unchanged three-graph cyclic helper resolves the
wider interval functions. Keep the summed correction upper bound,
separate origin accounting, complete pair gates, and every continuous
branch requirement intact. Bind the new coordinator and metadata in
the source contract and replay the complete verification.

No original certificate, helper, coordinator, or candidate source was
modified in preparing this variant.
