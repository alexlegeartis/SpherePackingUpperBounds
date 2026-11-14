LOWER-BOUND RESEARCH RUN — 14 NOVEMBER 2026

Outcome: no new sphere-packing density record was obtained. The registry and
the certified 36D construction are unchanged. This package contains exact
finite obstructions, finite coset screens and numerical screens. Exclusions of
specified finite families are not upper bounds on the sphere-packing problem.

The continuation is documented in FOLLOWUP.txt. It strengthens the common-fiber
obstruction to arbitrary nonlinear fibers at the original radius sqrt(2) for
the parent code and all 1350 shortenings, excludes all 3213 selected-coordinate fiber families with up to
four selected coordinates, and certifies local optimality of a cyclic-code
34D seed. The sections below describe the initial run; their narrower stated
scopes should be read with those subsequent results.

NORMALIZATION12.txt describes an independent, implemented construction target:
self-orthogonal ternary fibers at squared minimum twelve. The conditional
densities would improve records in 34D/36D/38D/40D, but no valid target was
found. Product tetracode fibers are excluded; the general search remains open.

A later continuation is in selforthogonal_followup/README.txt. It exactly
excludes the two saved best fiber types, adds shorter-dimension product and
Golay-factor obstructions, and supplies catalog, neighbor and SAT searches.
The general construction is still unresolved; no new record was obtained.

The broader-dimension continuation is in broad_dimensions_followup/README.txt.
It adds exact arbitrary-cluster antipode searches in 19--22D and 43--46D,
longer rank-two section searches, a verified 44D periodic stress obstruction,
and further catalog and 37D lamination screens. No new density record is
asserted. Its verify_all.py replays the short exact checks; --full includes
the longer rank-five searches.

Here delta always denotes center density: packing fraction divided by the
volume of the unit ball in the same dimension.

The first obstruction closes the common rank-three ternary-fiber proposal in
constructions/dim36/provenance/ct36_symmetry_break_20261112/README.md. It also
closes the corresponding 34D construction obtained by shortening one position.

EXACT RESULT 1: THE COMMON LINEAR TERNARY FIBER CANNOT WORK

Let C18 be the archived GF(4) [18,11,6] code, using the exact saved generator.
Let C17 be its shortening at coordinate 17 (zero-based), an [17,10,6] code.
Their weight-six/seven supports number 1757 and 1110, respectively.

There is no ternary matrix H with three rows, nonzero columns, and the property

    sum(epsilon_i H_i, i in S) != 0 for every epsilon_i in {1,-1}

for every weight-six/seven support S of either code. Rank at most three is
covered, including ranks one and two.

Why this matters geometrically. The archived mixed hexagonal construction has
local squared costs 0 for (0,0), 4 for (0,nonzero), 3 for (nonzero,0), and 1
when both symbols are nonzero; the period minimum is 12. For the product
C x ker(H), a zero column would give a norm-four difference. A signed zero sum
on a weight-six/seven support gives a difference of that same squared norm,
below eight. Conversely, in the absence of these obstructions, all other
nonzero differences have squared norm at least eight. Thus this is an exact
feasibility obstruction for that particular common linear fiber, not merely a
failure of a chosen sufficient test.

The tempting conditional center densities were:

  34D: 4^10 / 3^(23/2) = 3.417476..., versus 3^(33/2)/2^25 = 2.222034...
       Conditional gain 2^45/3^28 - 1 = 53.7994%.
  36D: 4^11 / 3^12 = 7.892322..., versus 2^24/(3^3*7^6) = 5.281629...
       Conditional gain 7^6/(4*3^9) - 1 = 49.4297%.

These targets are infeasible with the specified codes and a common check of
rank at most three. Neither target is a packing lower bound.

Proof by exhaustive finite search. Nonzero ternary columns modulo sign have
exactly 13 projective representatives. Independent column signs do not affect
the signed-sum condition. For each support, the search maintains the exact set
of reachable signed sums in F3^3. A complete support reaching zero is rejected.
If a partial support reaches all 27 states, every completion reaches zero.
When one column remains, all incompatible colors are removed from its domain.
All other branching choices are exhausted; domain and state changes are undone
on backtracking. The positive control is independently checked by directly
enumerating all signs.

The normalized search returns UNSAT after 23687 recursive calls for each
instance. A second build implements translations by explicit coordinate
arithmetic instead of bit permutations and reaches the same result. A further
run disables both column normalizations and searches all 13 choices for every
column, returning UNSAT after 1683150 calls for each instance. The 17-coordinate
obstruction already implies the 18-coordinate one by restriction to its first
17 coordinates. This explains the identical search counts.

This is a reproducible computer-assisted finite proof, not a formal proof in a
proof assistant or an audit by an independent person. It does not exclude a
different outer code, nonlinear ternary fibers, or fibers depending on the
outer word.

EXACT RESULT 2: EXPURGATING THESE HIGHER-RATE OUTER CODES ALSO FAILS

Even arbitrary nonlinear distance-eight subsets of C17 and C18 are too small
when used with the full ternary parity fiber.

  dimension  code size   forbidden degree  least eigenvalue  maximum subset  needed
      34       1048576       3330              -166              49789       75754
      36       4194304       5271              -241             183386      311875

Build the Cayley graph on the additive code, joining words whose difference
has weight below eight. Integer Walsh transforms enumerate every character
eigenvalue. For a d-regular graph with N vertices and least eigenvalue lambda,
Hoffman's bound gives alpha <= N*(-lambda)/(d-lambda). The exact rational
calculation is rounded down only at the final integer bound. A distance-eight
subset is independent in this graph. Its packing center density is
M/3^(n/2+1), where n is the outer length. The comparisons above are exact.

verify_codes.py checks injectivity, complete weight enumerators, spectral
moments, the least eigenvalue, support counts, exact record thresholds and
input/output hashes. No floating arithmetic is used in either finite proof.

NUMERICAL SCREENS — SEPARATE FROM THE EXACT RESULTS

1. Direct sections of the certified 36D lattice. The tested short dual vectors
   give 35D center density 2.47025577..., below the cited 2*sqrt(2). Coordinate
   subspaces of a reduced dual basis also gave weaker 32–34D sections. This was
   not an exhaustive search over all dual sublattices. The 35D contact matrix
   had numerical rank 630; the sampled metric linear program returned zero
   descent objective.

2. Recent 38D and 40D P48p complements. Returned norm-six vectors were checked
   with integer arithmetic. Floating shell searches returned 381658 and
   1092000 vectors, agreeing with their source certificates. The sampled metric
   programs returned zero descent objective. These LP outcomes were not turned
   into rational eutaxy certificates. No claim about all deformations or all
   ambient lattices follows. Incomplete early shells gave apparent descent
   directions that disappeared after more contacts were included.

3. Lamination into 37D. If delta36 = 16777216/3176523, a height h gives
   delta37 = delta36*sqrt(2)/h at squared minimum eight. Beating the cited
   delta37 = 4*sqrt(2) requires

       h^2 < delta36^2/16 = 1.7434753066907749...
       dist(x,L)^2 > 8 - delta36^2/16 = 6.256524693309225...

   The second inequality is necessary for the layer difference one. Layer
   difference two must also satisfy dist(2*x,L)^2 + 4*h^2 >= 8; for h^2 >= 8/9
   all layer differences at least three are automatically safe.

   A 100-shift CVP screen and twelve gradient/LP/cutting-plane hole searches
   reached squared-distance upper bounds up to 5.446526481972746..., below the
   required threshold. The exact reconstructed point has norm
   130655275/23988734. The zero lattice point is an exact upper-distance
   witness. Nearest-point completeness was not certified, and the search is
   not a covering-radius upper bound. At this point, even treating the
   numerical nearest distance as exact would yield only center density
   4.67430676..., below 4*sqrt(2). A follow-on joint 37D metric LP hit its time
   limit and supplies neither a valid improved candidate nor an obstruction.

Do not promote floating shell completion, CVP output, LP feasibility or the
exact norm of one nearest-point witness into a global minimum certificate.

REPRODUCTION

Python 3.11+, NumPy and a C++17 compiler are sufficient for the exact proofs.

    python -m pip install -r requirements.txt
    python verify_all.py

The replay checks a satisfiable positive control, both obstructions using two
translation implementations, and both obstructions with symmetry reduction
disabled. It uses a temporary directory and needs no network or solver.
Never run with Python assertions disabled. Saved output: results/exact_replay.log.

Optional numerical experiments require the separate experiments/requirements.txt:

    OPENBLAS_NUM_THREADS=1 python experiments/replay.py holes37
    OPENBLAS_NUM_THREADS=2 python experiments/replay.py section35
    OPENBLAS_NUM_THREADS=2 python experiments/replay.py p48_section38
    OPENBLAS_NUM_THREADS=2 python experiments/replay.py p48_bkz40_28

The floating enumerator aborts at its node cap. Rerunning an optional experiment
can overwrite its generated experiment files; historical results live separately
under results/. The large shells are regenerated, not shipped. The supplied
NPZ files contain exact integer Gram matrices and basis changes; their use by
floating discovery tools does not certify a new minimum.

NEXT CONSTRUCTION TARGET

Change the finite coset structure before optimizing the metric. The full
archived higher-rate outer codes with a common ternary fiber at the original
minimum-eight target are excluded by the continuation, even for nonlinear
fibers. A minimum-twelve target is explored in NORMALIZATION12.txt. The tested contact sets also
stalled the direct metric tests. Arbitrary simultaneous changes to outer words,
word-dependent fibers and metric are not excluded.

A concrete 34D target is five compatible cosets of a GF(4) [17,7,>=8] subcode.
They give M=81920 words, exceeding the threshold 75754 by about 8.14% in density.
Inside a dimension-ten parent, the quotient has only 64 cosets. Compatibility
can be computed by exact weight enumeration and the five-coset condition by
an exact clique search. Parents must change: the archived [17,10,6] code is
already excluded by the spectral bound above. This is a search target, not
an existence claim.

The supplied experiments/cyclic17_cosets.py implements this target. Starting
from the cyclic [17,8,8] code with polynomial generator

    (x+1) * (1 + 2*x + 2*x^3 + 2*x^4 + 2*x^5 + 2*x^7 + x^8),

over GF(4) with 2^2=3, it uses four codimension-one subcodes, all 17 single-
coordinate moves, all 408 normalized two-coordinate moves, and 256 sampled
three-coordinate moves. Each move extends the odd-like dimension-nine parent
to dimension ten. All 64 coset minima are enumerated exactly, and an exhaustive
five-clique search is run. All 2724 tested parent instances reject five cosets;
some sampled moves may repeat. This does not exclude other parents, other
subcodes, or unrestricted code switching. Saved results are under results/.
Reproduce with: python experiments/cyclic17_cosets.py (NumPy only).

A surviving next direction is to change the parent/glue code more substantially,
then apply a metric or periodic-position optimization. That recommendation is
a research direction supported by these failures, not a claim of a new method
or of likely success.

PROVENANCE AND REFERENCES

Repository inspected and baseline used:
https://github.com/alexlegeartis/SpherePackingUpperBounds/tree/72028bb0c5eaaaeede50cd418ca09339ee404285

Archived code/fiber proposal and mixed-cost provenance:
https://github.com/alexlegeartis/SpherePackingUpperBounds/blob/72028bb0c5eaaaeede50cd418ca09339ee404285/constructions/dim36/provenance/ct36_symmetry_break_20261112/README.md

Certified 36D baseline:
https://github.com/alexlegeartis/SpherePackingUpperBounds/blob/72028bb0c5eaaaeede50cd418ca09339ee404285/constructions/dim36/README.md

38–43D comparison/construction source:
Ivan Dorofeev, Xiaoming Sun, Chengu Wang, arXiv:2607.20359v4 (29 August 2026)
https://arxiv.org/html/2607.20359v4
Certificate repository pinned to a0bc7f51444e42c3cbae0d2a910039cd50013500:
https://github.com/wcgbg/sphere-packing/tree/a0bc7f51444e42c3cbae0d2a910039cd50013500

Direct access to Cohn's live density table returned HTTP 502. Comparisons use
the explicit pinned references above and the repository's cited 34D/35D/37D
values. No exhaustive literature-priority claim is made or needed for these
negative finite-family comparisons.
