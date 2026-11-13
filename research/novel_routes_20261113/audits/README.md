# Independent exact audits and propagated benchmarks

This directory contains compact certificates and independent exact verifiers for several constructions found during the November 2026 search. The ambient lattice minimum theorems are inherited from their cited primary sources. Exact checks here establish the resulting packing lower bounds; they do not establish publication priority or a complete current record table.

## What has passed

| Construction | Exact checker | Result |
|---|---|---|
| Ne64 antipodes in dimensions 60–63 | `audit_ne64_antipodes.py` | Source Gram identity, dual scaling, primitive normal, complete section period, all projected shifts, diameter and density formula passed. |
| Ne64 bare sections in dimensions 57–59 | `verify_ne64_sections.py` | Source identity, D5/E6/E7 normal Grams, saturation, complete section periods and inherited density passed. |
| Γ72 sections from laminated normals, dimensions 56–71 | `verify_gamma72_leech_sections.py` | Source identity, ambient determinant and positivity, all normal Grams and saturation passed. Complete explicit period bases were also checked for the 58D and 59D sections. |
| Finite propagation benchmark envelope, dimensions 57–71 | `propagated_benchmarks.py` | Exact rational squared-density lower bounds from products, torus coloring, averaged projection and the preceding Γ72 sections. |

Run all four checkers with Python, NumPy and SymPy. The propagation checker also uses mpmath to propose rational radii; the final inequalities are verified with exact rational arithmetic. No search pool, compiled enumerator, fpylll installation, or floating short-vector search is needed to replay these audits.

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/audit_ne64_antipodes.py
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/verify_ne64_sections.py
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/propagated_benchmarks.py
```

The verifiers read archived primary source pages and compact certificates in the sibling `dual64` and `higher_parent` directories. The Γ72 period files store exact integer entries as decimal strings to avoid fixed-width overflow; NumPy does not need pickle to load them.

## A nonunimodular antipode detail

For Ne64, let L have Gram G, minimum 12 and covolume 3^16. A primitive row matrix A in dual integer coordinates defines K = A G^−1 Aᵀ. Primitivity gives L* intersect span(K) = K, hence the projected lattice is K*. The orthogonal section has covolume 3^16 sqrt(det K), not a covolume obtained by treating the ambient lattice as unimodular.

For dual labels whose squared diameter is beta < 12, projected cosets are distinct: equality would give a nonzero ambient lattice vector lying entirely in the normal space with squared norm at most beta. Every projected difference has squared norm at least 12 − beta. These facts justify the antipode density formula used in the exact checker. Explicit integer right inverses, complete kernel bases and rational projected shifts are retained in the small `*_explicit_period.npz` files.

## Benchmark corrections

The initially found Ne64 bare sections in dimensions 57–59 are valid, but all are superseded by bounds in this audit:

- Dimension 57: a torus extension of the published Schütt 56D lattice gives a rigorously certified center density above 151,706.08. The Ne64 section gives only 57,530.18.
- Dimension 58: the Γ72 normal has determinant 3·2^22, giving center density 2^18/√3 = 151,348.90897…. This exceeds the Ne64 section's 115,060.35.
- Dimension 59: the Γ72 normal has determinant 2^23, giving center density 262,144. This exceeds the Ne64 section's 244,079.86.
- Dimension 63: the Ne64 antipode gives center density about 13.106 million, below the older lower bound printed in Vardy's 1999 paper, whose log₂ center density is at least 23.92….
- Dimension 65: the initial Γ72 antipode is below the separately proved Ne64 discriminant-count extension. The newly certified explicit Ne65 lattice has center density 3^17/√10 = 40,837,705.24854…, and is displayed as a candidate rather than a prior benchmark.

The candidates in dimensions 60–62 and 65–71 exceed the finite benchmark families audited here. This does not by itself establish global SOTA. `PROPAGATED_BENCHMARKS.md` and `propagated_benchmarks.json` state the complete scope and keep Vardy's rounded published values separate from exact constructions.

## Why the Γ72 sections exist

The first principal 24×24 block of the published Γ72 Gram has even entries, determinant 2^24 and minimum at least 8. Dividing its Gram by 2 yields an even unimodular 24D lattice of minimum at least 4: the Leech lattice. Its coordinate embedding in Γ72 is primitive.

A short discovery search selected primitive rows of this Leech block with the standard laminated determinants. The compact rank-1 through rank-16 matrices are retained as `gamma72_leech_section_rank*.npz`. Their correctness is verified directly from integer matrices; completeness of the discovery search is irrelevant. The rank-13 and rank-14 determinants of the scaled normal Gram are 2^23 and 3·2^22, respectively. Since Γ72 is unimodular, a primitive normal and its orthogonal section have equal Gram determinants. The inherited minimum 8 then gives the displayed lower bounds.

## Mordell–Weil caution

`MW64.html` archives the catalogue's MW64 Gram, of determinant 2^52 and published minimum 12. Vardy's printed 57–63D numbers resemble a simple exact determinant sequence, but that identification has not been established. The available MW64 basis did not immediately supply the required unit dual normals. The empty floating enumeration recorded in `mw64_unit_dual_enum.json` is a discovery observation, not a proof of absence. No inferred exact Mordell–Weil density formula is used in the audit.

## Primary sources

- Ne64 Gram and published lattice data: https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64b.html
- Γ72 Gram: https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html
- Gabriele Nebe, *An even unimodular 72-dimensional lattice of minimum 8*: https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/Ext72.pdf
- Current published parent density inputs, including exact Schütt 56D density: https://cohn.mit.edu/sphere-packing/
- Small-dimensional density inputs: https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html
- Alexander Vardy, *Density doubling, double-circulants, and new sphere packings*, Trans. AMS 351 (1999), 271–283, Section 4: https://citeseerx.ist.psu.edu/document?doi=c9085ab03f7ab8f5154d3c38fe2be7ab6316fac4&repid=rep1&type=pdf

Large search pools, raw enumeration grids and routine logs are ignored. All seven 32D quotient certificates remain retained in the sibling `transfer_neighbors` directory, together with their exact all-class verifier and the resulting 33D family exclusion theorem.
