# Lower-bound research across dimensions, 15 November 2026

Two routes now have proofs: **compatible deep half-cosets** above published
72D/80D lattices, and **explicit antipodes** below L80. The former supplies
existence bounds in 73–78D and 81–85D; the latter supplies explicit finite
certificates in 76–79D. The best bounds from this experiment cover twelve
dimensions. The general self-orthogonal minimum-12 construction remains
open; these routes do not depend on resolving it.

**No new worldwide record is established here.** Every comparison below is
a specified construction derived from named published parents. It is not
a dimension-specific value printed by the source, or an assertion that no
previous packing is denser. In particular, possible higher-minimum sections
and sections of the denser Mordell–Weil 80D parent remain outside this
comparison. Sparse published tables do not establish priority.

The [15 November worldwide-status check](LITERATURE_STATUS.md) now identifies
**Chen's published 84D/85D guarantees** as direct comparisons: our bounds
exceed them by **899.675713%** and **364.438553%**, respectively. No later
stronger guarantee in those two dimensions was located in that search.
They are **record candidates**, with the remaining priority gaps stated in
the review. [Exact external-comparison replay](verify_external_comparisons.py).

## Best bounds in this experiment

All values here are **center densities** `delta`, not occupied-volume
fractions. The fraction is `delta*pi^(n/2)/Gamma(1+n/2)`. Exact rational
density squares, comparison formulas and proof inputs are in the linked JSON
reports. Ratios have the same value for center densities and fractions.

`N12` means [Nebe's Gamma72, minimum 8 and determinant 1](https://arxiv.org/pdf/1008.2862).
`BN98` means [Bachoc–Nebe's L80, minimum 8 and determinant 1](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html).
The comparisons in 73–83D are **derived orthogonal products** of the stated
parent with the listed classical root lattice, with both factors at minimum 8.
For a parent of dimension m its center-density square is
`2^m/(2^r*det(R_r))`, where `r=n-m`. The small Gram determinants are
`2,3,4,4,4,3,2` for `A1,A2,A3,D4,D5,E6,E7`.
The 84D/85D rows instead use the dimension-specific guarantees in
[Chen (2013), Theorem 2.4 and §3](https://doi.org/10.5802/jtnb.825).

<!-- results-table:start -->
| Dimension | Center density delta | Construction / proof | Comparison method / source | Ratio to that comparison |
| --- | ---: | --- | --- | ---: |
| 73 | 5.610922e+10 | [Compatible shifts (existence)](counting_results.json) | N12 Gamma72 x A1 | 1.63299x |
| 74 | 5.290028e+10 | [Compatible shifts (existence)](counting_results.json) | N12 Gamma72 x A2 | 2.66667x |
| 75 | 5.290028e+10 | [Compatible shifts (existence)](counting_results.json) | N12 Gamma72 x A3 | 4.35465x |
| 76 | 6.108398e+10 | [Compatible shifts (existence)](counting_results.json) | N12 Gamma72 x D4 | 7.11111x |
| 77 | 6.108398e+10 | [Compatible shifts (existence)](counting_results.json) | N12 Gamma72 x D5 | 10.05663x |
| 78 | 1.047318e+11 | [Explicit antipode](l80_results.json) | N12 Gamma72 x E6 | 21.11784x |
| 79 | 2.991785e+11 | [Explicit antipode](l80_results.json) | N12 Gamma72 x E7 | 69.65793x |
| 81 | 9.834331e+11 | [Compatible shifts (existence)](counting80_results.json) | BN98 L80 x A1 | 1.78885x |
| 82 | 1.015685e+12 | [Compatible shifts (existence)](counting80_results.json) | BN98 L80 x A2 | 3.20000x |
| 83 | 1.112628e+12 | [Compatible shifts (existence)](counting80_results.json) | BN98 L80 x A3 | 5.72433x |
| 84 | 1.407375e+12 | [Compatible shifts (existence)](counting80_results.json) | [Chen 2013 generalized Craig lattice + binary code; theorem guarantee](LITERATURE_STATUS.md#direct-published-comparison-in-84d-and-85d) | 9.99676x |
| 85 | 1.453531e+12 | [Compatible shifts (existence)](counting80_results.json) | [Chen 2013 generalized Craig lattice + binary code; theorem guarantee](LITERATURE_STATUS.md#direct-published-comparison-in-84d-and-85d) | 4.64439x |
<!-- results-table:end -->

The existence bounds do **not** contain explicit horizontal representatives.
They prove that suitable representatives exist by a strict counting
inequality. The explicit L80 certificates contain all data needed to define
the periodic packings; their minimum bound inherits the published parent
minimum. Neither checker formally proves the external parent theorem.

## The promising idea: preserve compatibility while counting deep shifts

An isolated deep hole helps in one more dimension. Several holes help in
several more dimensions only if **all their differences** remain deep.
Here we enforce that compatibility through a small binary subspace of
`L/2L`. The theta series bounds the number of bad classes. Averaging over
the finite orthogonal group proves that a subspace of a chosen quadratic
type avoids every bad class.

The quadratic type matters. For Gamma72, singular classes can be made to
have half-coset minimum at least 5, and nonsingular ones at least 9/2.
Matching the vertical metric to those two values yields minimum-8 lattices.
For example:

* In 76D and 77D, `delta = 2^39/9 = 61,083,979,320.888...`.
* The mixed-depth 77D construction improves the **uniform-depth variant
  derived in this experiment** by about 27.32%; this is not an external
  record comparison.
* Above L80 the corresponding depths are 5 and 11/2. The 85D proof needs
  A5's quadratic type; D5 fails the averaging inequality even though its
  root-lattice determinant is smaller. In 85D the construction gives
  `delta = 2^47/(25*sqrt(15)) = 1,453,530,529,565.988...`.

The complete [counting and gluing proof](COUNTING_PROOF.md) includes the
exact bad-class counts, subspace models, group-averaging argument, vertical
metrics, layer cases and determinant calculation. Arithmetic is checked by
two standard-library scripts. The Delta series is calculated in two ways;
the norm-8 coefficients agree with the published kissing numbers.

## Explicit L80 antipodes

This uses the established [Conway–Sloane antipode construction](https://arxiv.org/abs/math/0207182),
not a new construction theorem. Integer rows B specify a primitive normal
sublattice of L80 with Gram `K=B G B^T`. A set S of integer dual labels
has squared diameter beta in metric `K^-1`. The packing has period Gram
determinant `det(K)` and guaranteed squared separation `8-beta`, giving

```
delta^2 = |S|^2 * ((8-beta)/4)^n / det(K).
```

| Dimension | Normal rank | det(K) | Cosets | beta | Guaranteed squared separation |
| --- | ---: | ---: | ---: | ---: | ---: |
| 79 | 1 | 10 | 2 | 1/10 | 79/10 |
| 78 | 2 | 48 | 3 | 1/6 | 47/6 |
| 77 | 3 | 300 | 4 | 1/5 | 39/5 |
| 76 | 4 | 2000 | 5 | 1/5 | 39/5 |

The projection formula follows directly: every difference of ambient
points has squared norm at least 8 and loses at most beta in the normal
space. Distinct labels give distinct translates; otherwise a nonzero
parent vector would lie entirely in the normal space with norm below 8.
Within one translate the period minimum is at least 8.

`verify80.py` re-parses the archived source, checks positive definiteness
by exact leading principal determinants and checks determinant 1. Smith
forms verify saturation of B. An independent integer HNF calculation gives
a full period basis N and integral lifts R satisfying
`N G B^T=0`, `R G B^T=I`, and `det(N G N^T)=det(K)`. A label s gives
the explicit shift `s R - s K^-1 B`. Every cluster distance and resulting
density square is checked rationally. Minimum attainment after projection
is not asserted. A separate SymPy LDL audit also verified the parent Gram
and each normal Gram; the main checker uses FLINT for large matrices.

For comparison with the earlier section method, exact small-dimensional
Hermite constants give ceilings on **bare-section guarantees retaining
parent minimum 8**. These four antipodes exceed those method ceilings by
8.84%, 31.99%, 39.41%, 36.70%, respectively. Such a ceiling is not a
ceiling for sections whose actual minimum is higher, and is not an external
density record. The 76D/77D antipodes are weaker than the counting
constructions above; they remain useful explicit alternatives.

## What the unsuccessful searches taught us

The compact [discovery reports](discovery_reports/) retain actual budgets
and results. Capped or empty searches prove no nonexistence statement.

* **57–58D:** 1,000 extensions of the actual 59D normal were evaluated,
  using finite dual-point enumeration and greedy clusters. Best floating
  scores were 24,628.17 in 57D and 83,479.26 in 58D, below the existing
  derived references 151,706.08 and 151,348.91. The references come from
  a [Schütt56 torus extension and a Gamma72 section](../novel_routes_20261113/audits/PROPAGATED_BENCHMARKS.md).
  There were 480 point-cap skips; this was not exhaustive.
* The old **676-ray Ne64 pool spans only rank 44**. Exact rank and all
  pair operations at half the minimum inner product were checked: these
  operations add zero rays. Repeating the same closure cannot fix the
  missing twenty dimensions.
* **66D:** six exact catalogue automorphisms enriched a Gamma72 minimum
  shell sample to 450,000 unsigned rays. A distributed search tried 112
  roots, at most 6,000 attempts per root, without embedding the earlier
  abstract normal. It remains a conditional candidate, **not a packing**.
* **M80 control:** the sampled pool found the same 79D, 78D and 76D
  antipode parameters. It did not find the rank-3 target in its finite
  pool. Only the audited L80 certificates are retained as proof data.

Large sampled Gamma72 pools and M80 control data are scratch artifacts;
they are not needed to replay any positive result. The L80 source, normal
certificates and complete arithmetic reports are committed. The initial
L80 seed/closure pool generator is included and reproduces the parameters.

## Next experiments, in priority order

1. **Find explicit compatible shifts in 73–77D and 81–85D.** The proof
   gives useful success-probability bounds for random finite isometries:
   at least 31.66% for the 77D model, and 9.02% for the 85D model.
   The hard part is certifying half-coset depths. Use the Hermitian/Leech
   description of Gamma72 to reduce the nearest-point work to lower
   dimensional pieces, and cache short counterexamples as constraints
   while sampling. A coordinate certificate would turn an existence
   result into a reproducible explicit packing.
2. **Search the right fiber for the conditional 66D normal.** Its first
   two norm-8 vectors x,y have inner product -3, so `u=x+y` has norm 10
   and the other four normal vectors must lie in `u^perp`. Generate the
   fiber `w^2=8, w.u=5` rather than another uniform shell sample. If
   `v^2=8`, `v.u=0`, `w.v=+/-4`, then `w-sign(w.v)*v` is another
   point of the same fiber. Differences can generate new short vectors
   in the required perpendicular lattice. An embedding and saturation
   certificate are still required before claiming the conditional gain.
3. **Build a full-rank Ne64 dual pool for 57–59D.** Use independently
   verified parent automorphisms or additional short-vector seeds that
   escape the rank-44 span, then repeat the multi-point cluster search.
   Retain the stronger classical 59D comparator and the named derived
   references in 57–58D. Larger repetitions inside the closed old pool
   have little justification.

## Reproduce

The existence checks require only Python's standard library:

```sh
python research/cross_dimension_lower_20261115/verify_counting.py
python research/cross_dimension_lower_20261115/verify80_counting.py
python research/cross_dimension_lower_20261115/verify_external_comparisons.py
```

The explicit auditor needs NumPy, SymPy and python-flint:

```sh
python research/cross_dimension_lower_20261115/verify80.py
```

Use Python with assertions enabled. Verification reads saved results and
rejects a stale report. `--write-result` regenerates a report after a
deliberate change. The published parent minimum theorems, theta modularity,
classical finite quadratic-space classification/transitivity, and root
lattice properties are stated mathematical inputs, not formally verified
theorems in this repository.

Optional L80 discovery replay, writing outside the repository:

```sh
python research/cross_dimension_lower_20261115/prepare80_pool.py --output /tmp/L80_pool.npz
python research/cross_dimension_lower_20261115/search80.py --pool /tmp/L80_pool.npz --output /tmp/L80_normals
```

The larger discovery scripts accept user-supplied sampled pools and depend
on NumPy/SymPy/SciPy; they are not part of certificate replay.

## Sources and literature scope

* **N12:** G. Nebe, [*An even unimodular 72-dimensional lattice of minimum
  8*](https://arxiv.org/abs/1008.2862), J. reine angew. Math. 673 (2012),
  237–247. Theorem 2.3 supplies the parent, and p. 2 gives its kissing
  number. The catalogue automorphisms are archived in the earlier
  [Gamma72 package](../novel_routes_20261113/higher_parent/).
* **BN98:** C. Bachoc and G. Nebe, [*Extremal lattices of minimum 8
  related to the Mathieu group M22*](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf),
  J. reine angew. Math. 494 (1998), 155–171.
  [L80 catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html),
  [M80 catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/M80.html).
  The normalized L80 source and SHA-256 hashes are retained.
* **Antipodes:** J. H. Conway and N. J. A. Sloane,
  [*The Antipode Construction for Sphere Packings*](https://arxiv.org/abs/math/0207182),
  Invent. Math. 123 (1996), 309–313. The rank-2/3/4 cluster shapes here
  reuse the repo's earlier Gamma72 experiments; their embeddings in L80
  are explicit. [Chen et al. (2025)](https://arxiv.org/abs/2505.02394)
  also use suboptimal normal sections; their stated dimensions differ.
* **Related counting work:** G. Nebe and C. Xing,
  [*A Gilbert–Varshamov Type Bound for Euclidean Packings*](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/GVlat.pdf),
  Math. Comp. 77 (2008), 2339–2344. It uses lattice-point counts and
  theta series to obtain GV bounds. The author's
  [Leech/Gamma72 calculation](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/xingleech72gross)
  concerns dimensions divisible by 720. This is relevant prior counting
  methodology; the present experiment does not claim to originate
  counting-based packing constructions.

This research supplies concrete bounds and experiments. Establishing any
new SOTA requires comparison with the strongest prior constructions in the
same dimension, beyond the product and section-method controls used here.
