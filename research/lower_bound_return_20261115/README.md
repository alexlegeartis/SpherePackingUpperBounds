# Lower-bound return, 15 November 2026

**A new verified 59D packing improves the repository's certified 59D benchmark
by 4.722057%.** Center density rises from `262144` to
`274522.589741149678…`; packing fraction is
`2.637602590720003…e-12`. This is a six-translate antipode over a rank-five
primitive normal in Ne64's dual. The explicit construction, two exact proof
implementations, attained-minimum witness, and full comparison audit are in
[constructions/dim59](../../constructions/dim59/README.md). Worldwide priority
remains unestablished.

## Successful route: balanced path normal

Earlier Ne64 work retained the best normals of ranks one through four and
several higher-rank root sections. Here the rank-five simplex search tests all
`2^10=1024` numerator Grams with diagonal 12 and off-diagonal entries -2 or -3.
The successful graph of -3 entries is a five-vertex path. The normal numerator
is `14 I - 2 J - Adj(P5)`; the actual dual Gram is its third.

The stored 676-vector dual minimum pool suffices to find an explicit primitive
embedding after 105,692 backtracking attempts. The finite pool and floating
scores only discover it. Exact construction proofs re-read the original Ne64
source and do not trust the pool or any search result.

```bash
python research/lower_bound_return_20261115/ne64_rank5.py --output /tmp/ne64_rank5.npz
python constructions/dim59/verify.py
python constructions/dim59/verify_independent.py
```

The 59D comparison is to the repository's certified Gamma72 rank-13 section,
not merely to its weaker former Ne64 bare section. The second checker
reconstructs that comparison period and its determinant independently.

## Other routes and their actual status

* **65D, half-cosets:** the preceding branch left 88 of 100 refined Ne64
  classes unresolved. Pruned nearest-first enumeration finds an exact
  numerator-norm-at-most-34 witness for every one of the 88. Together with the
  earlier twelve witnesses, all 100 archived candidates fail the depth-nine
  target. This is a finite candidate exclusion, not an exclusion of all
  `2^64` half-cosets or of general holes. `half_coset_rejections.npz` archives
  the complete exact upper-witness list.
* **65D, another dual type:** 10,000 classes of norm `4/3 mod 2` were sampled
  with seed 651115. Exactly 8,497 have archived witnesses of numerator norm
  at most 66. The other **1,503 remain unresolved**, having reached a
  discovery cap. If one had minimum `28/3`, it would improve the existing
  65D layer height, but no such hole is certified here.
* **66D, eight-point circuits:** continuous normal optimization and an
  integer rounding search produce a positive integer rank-six Gram of
  determinant 21875 and an eight-label cluster of squared diameter `2/5`.
  Its *conditional* density is `85503676.224…`, above the existing 66D
  value. **No embedding into Gamma72 was found. This is not a packing
  certificate or a new lower bound.** Three complete finite searches used
  two million attempts each in a 150,000-vector unsigned pool; a fourth was
  interrupted. These searches do not prove non-embedding.
* **65–69D circuit relaxations:** fixed-cluster numerical searches include
  seven-, eight- and nine-label alternatives and point deletions from the
  former root clusters. `circuit_relaxations.json` and the saved root
  relaxations are exploratory optimizer output, not exact certificates.
* **58–63D Gamma72 sections:** dual-cluster searches on the saved laminated
  normals do not improve the best existing values. Some larger-diameter
  clique problems hit a 300,000-node cap. No full-family exclusion is asserted.

The general self-orthogonal minimum-12 code construction from the earlier
follow-up remains open. This run's success is a different construction in a
different dimension.

## Replaying the screens

```bash
python research/lower_bound_return_20261115/verify_screens.py
g++ -O3 -std=c++17 -fno-fast-math -ffp-contract=off \
  research/lower_bound_return_20261115/hole_screen.cpp -o /tmp/hole_screen
```

`verify_screens.py` checks every saved upper witness by integer multiplication,
its residue class, and the rejection threshold. Rows marked unresolved have
no lower-distance interpretation. It also checks the abstract 66D Gram and
diameter, explicitly retaining its unembedded status.

`hole_screen` takes `gram input_classes output_jsonl node_cap modulus cutoff`
and an optional pruning-slack argument. Each input line is a class ID followed
by 64 residues in the reduced basis used by
`../novel_routes_20261113/higher_parent/ne64_primal_best.npz`.
The successful half-class replay uses cutoff 34, modulus 2, cap 10,000,000 and
pruning slack 2. The second dual-type screen uses cutoff 66, modulus 3, cap
1,000,000 and slack 2. Without a witness, an empty floating or pruned traversal
still does not certify a hole.

`normal_relaxation.py --rank 6 --circuits --output /tmp/circuit6.json`
replays the exploratory circuit optimization. `embed_normal.py` documents the
finite Gamma72 embedding search; its vector pool is generated on demand and
is not committed. Use a fresh output path for large discovery artifacts.

Exact construction replay needs NumPy, SymPy and python-flint. Numerical
discovery additionally needs SciPy and fpylll. Machine-readable scope and
counts are in `search_scope.json` and `screen_audit.json`.
