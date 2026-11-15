# Fresh reproduction of the recent 88D computation

**Completed:** all stages of the author's full computational proof were
executed successfully in a fresh output directory on 15 November 2026.
Total wall time was **1697.53 seconds**. This removes the previous gap
“parent minimum computation not replayed” for the parent used by our
86D/87D antipode certificates.

Source: Ryotaro Sakamoto,
[rank-88-extremal-lattice-computations](https://github.com/RyotaroSakamoto0113/rank-88-extremal-lattice-computations/tree/628ce2a7ac23acf8d4c295677e4310c9ec5b21d9),
commit `628ce2a7ac23acf8d4c295677e4310c9ec5b21d9`, associated with
[arXiv:2609.15633v2](https://arxiv.org/abs/2609.15633v2).
The checkout's tracked files were unmodified.

## What ran

The command was the upstream full pipeline, rather than its quick
certificate/hash check:

```sh
python scripts/reproduce_all.py --output /path/to/new/full88_run
```

The eight top-level stages were L, J0, groups, factor_claims, L0_input,
L_perfection, L0_perfection and formulas. All returned zero. L completed
its projective-minimum branches, exact certificates, rank-4 obstruction,
field and ideal computations, orbit checks, histograms and catalogue
binding. All **78 independent branch verification reports** passed.
The J0 run accounted for all 8,159 rank-4 targets, found no forbidden
embedding, and concluded trace minimum 8 using the earlier rank-1–3
bounds. The group, factor, construction and perfection checks also
passed.

The freshly rebuilt L and L0 lattice files equal the distributed source
files exactly. More specifically, the rebuilt **L0 Gram equals the G
stored in both of our 86D/87D antipode certificates entry for entry**.
This binds the completed minimum computation to the actual parent of
our packings; matching just a paper title or a determinant would not
establish that connection.

The replay used Python 3.12.14, NumPy, PARI/GP 2.17.4, GAP 4.14.0 with
SmallGrp 1.5.4 and GAPDoc 1.6.7, GNU C++ 14.2.0 and GMP 6.3 on Debian 13.
PARI and GAP were extracted into a workspace-local sysroot; small
wrappers set their library/package paths. The upstream clang++ calls
were forwarded to the installed GNU compiler. No source-code change
was needed.

## Evidence and scope

The [top-level result](replay88/verification.json),
[L final verification](replay88/L/FINAL_VERIFICATION.json),
[J0 summary](replay88/J0/output/rank4_audit/FULL_RUN_SUMMARY.json),
and other compact stage reports are archived here.
The [manifest](replay88/manifest.json) records the pinned source-file
hashes, hashes of fresh stage reports and both fresh lattice files,
and all 78 branch-report hashes. Raw outputs occupy about 1.8 GB in
the execution workspace and are not required to interpret the verdict.
To regenerate this compact archive after your own full run:

```sh
python research/priority_resolution_20261115/archive_replay88.py \
  --upstream /path/to/pinned/checkout --run /path/to/new/full88_run
```

This is a **fresh reproduction of the author's computational proof**.
It is not an independently designed minimum algorithm or a formal
verification of the analytic reduction. The source's mathematical
reductions and its implementation remain dependencies. Credit for the
parent construction and minimum proof belongs to Sakamoto.

The older [conditional algebra report](../cross_dimension_lower_20261115/claimed88_comparisons.json)
is retained as a historical standalone audit: that checker itself only
checks the antipode algebra and does not perform this full replay.
Combining it with this replay discharges its computational parent-minimum
gap. Nevertheless **both antipodes are weaker than classical BCS88
sections**, so neither is a record candidate. The classical comparison
does not require this new source at all.
