SELF-ORTHOGONAL MINIMUM-TWELVE SEARCH: CONTINUATION

No new certified density record. The original best fibers had 404 and 262
norm-nine bad trit words in 34D and 36D. The new search reduced the 36D count to
258, with no density change: the attained minimum remains nine. The general
construction is open; the registry and certified 36D packing were not modified.

EXACT EXCLUSIONS

Write K for the ternary [4,2,3] tetracode, and N for the self-orthogonal [9,4,3]
code with generator

    111000000
    000111000
    000000111
    120120120

N has enumerator 1+6*z^3+66*z^6+8*z^9. Its weight-three supports are three
disjoint triples. Its weight-six supports are the three unions of two triples
and the 27 choices of two positions from each triple. Its weight-nine support
is the full nine-coordinate set.

The original saved 34D fiber is exactly N+K+K, and both the original and the new
36D fiber are N+N, where + denotes direct sum. The replay checks these
decompositions from the full codes.
It excludes every coordinate placement of these types, including all 190
two-position and 1140 three-position shortenings of the parent outer code.
Every length-17 instance checks all 6806800 N placements; at most two are
internally valid, and all 35 complementary K+K partitions fail for each.
Every length-18 instance checks all 13613600 placements; none is internally
valid. Column signs do not change the compatibility conditions.

Further shortening gives conditional minimum-twelve targets in 30D and 32D:

    dimension    outer code       ternary dimension    conditional density
        30        [15,8,6]                7                 2/sqrt(3)
        32        [16,9,6]                8                 4

These exceed the pinned references 1 and 3^16/2^24. However, K+K+K plus the
length-three repetition code fails for every one of the 15504 five-position
shortenings, and K+K+K+K fails for all 4845 four-position shortenings. The tests
include weight-nine conditions and conditions involving the repetition block.
Neither target was achieved.

The 206 verified support-symmetry orbits cover 21679 outer codes in these
explicit families. The group has order 120 and is checked on EVERY outer
weight-six, weight-seven and weight-nine support. The parent dual minimum is
exactly eight, so these shortenings have the claimed rates. A symmetry need
only preserve those support families; full-code automorphism is not assumed.

A separate obstruction excludes any ternary Golay12 direct factor. Every one
of the parent's 125970 twelve-coordinate subsets contains an outer weight-seven
support. Every seven-subset of the ternary Golay [12,6,6] code contains exactly
one of its 132 weight-six supports. Together these give a norm-nine difference
for any Golay-factor placement, including in the shortened outer codes.

Replay all these exclusions with NumPy and a C++17 compiler:

    python verify_exclusions.py --output work/exclusions.json

Saved results: results/exclusions.json and results/exclusions_replay.log.
The replay includes positive and negative controls. Do not disable assertions.

CATALOG SEARCH

The Harada--Munemasa database supplies these maximal self-orthogonal catalogs:

    length       15    16    17    18    19    20
    dimension     7     8     8     8     9    10
    code types   12     7    23   160    56    24

All 282 matrices are checked for rank and self-orthogonality. Completeness is
a published classification input, not independently reproved here. Provenance
and source hashes are in inputs/catalog_sources.json.

Each type was screened by a coordinate-permutation tabu search at its
30D/32D/34D/36D/38D/40D target. No valid minimum-twelve placement was found.
These are heuristic searches, not exhaustive permutation checks.

    python verify_catalogs.py
    python catalog_search.py --dimension 34 --steps 10000

Within the published length-17 catalog, types 0, 1 and 6 are excluded for every
placement by the tetracode, N-component and Golay arguments. The other 20 types
remain unresolved. Catalog placement results did not improve the old best
34D/36D fibers. Audited best placements are under results/catalog*.

TARGETED NEIGHBORS

An integer Fourier transform over F3 scores every hyperplane of a current
fiber. The search chooses low-cost cuts and enumerates their self-orthogonal
extensions. Dimension-seven cuts have four maximal extensions in length 17
and ten in length 18. Dimension-six cuts have 40 and 280, respectively.
The latter search caches projective-coset costs and sums them over isotropic
planes. Sampled cached energies are checked against direct code enumeration.

There were 2801557 extension/completion evaluations, including repeats. The
34D violation count stayed at 404; the 36D count improved from 262 to 258.
Every completion of a selected cut is enumerated; selection of cuts is
heuristic. This is not a global exclusion or an improved packing density.

    python neighbor_search.py 17 --seconds 180
    python neighbor2_search.py 18 --seconds 180
    python audit.py results/neighbor2_18_best.txt --output work/audit.json

Best generators are written immediately on discovery. The optional fifth
header field is the norm-six penalty (1000 for new neighbor outputs; 10 for
four-field catalog/older outputs). Use the exact audit's attained minimum and
density for any claim.

BOOLEAN MODEL

Install requirements_sat.txt to use sat_model.py. It encodes a coordinate
bijection, all weight-six conditions, verified support-automorphism symmetry
constraints, and lazily added weight-nine conditions.

    python sat_model.py 17 0 --conflicts 100000 --proof

For the already excluded tetracode type 0, Glucose returned UNSAT and an
independent drat-trim check returned VERIFIED. Source/CNF/proof hashes are in
results/sat_known_case_verified.json. Raw proofs are regenerated by the command.
The direct exclusion replay does not require SAT or DRAT software.

The other 20 length-17 types all returned UNKNOWN at 50000-conflict budgets.
UNKNOWN is not UNSAT. No new exclusion or record follows from those runs.
Empty-target controls are marked explicitly and are not packing candidates.
A new SAT UNSAT result should be independently proof-checked before being used
as a certificate.

STATUS AND SOURCES

The two saved fiber types are dead ends for minimum twelve in this outer
family. Other self-orthogonal types and other outer/fiber constructions remain
open. No new sphere-packing lower-bound SOTA is claimed.

https://www.math.is.tohoku.ac.jp/~mharada/research/codes/mso3.htm
https://www.math.is.tohoku.ac.jp/~munemasa/research/codes/data/3/17.magma
https://www.math.is.tohoku.ac.jp/~munemasa/research/codes/data/3/18.magma
M. Harada and A. Munemasa, On the classification of weighing matrices and
self-orthogonal codes, Journal of Combinatorial Designs 20 (2012), 40--57.
https://arxiv.org/abs/1011.5382
