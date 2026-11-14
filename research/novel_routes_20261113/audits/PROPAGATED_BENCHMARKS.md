# Propagated benchmark audit, dimensions57–71

Finite propagated lower-bound audit, NOT an exhaustive SOTA/priority determination. It does not bound all sections, laminations, nonlattice packings or constructions from other parents.

All densities below are center densities. Principal computed bounds have exact rational squared densities in `propagated_benchmarks.json`. The torus bounds use strict rational volume inequalities, not floating hole searches.

| Dimension | Strongest bound in audited finite families | Construction | Current candidate | Vardy1999 printed log₂δ |
|---|---:|---|---:|---:|
| 57 | 151,706.083320 | 56D torus offsets, 1D lattice modulo 2 | 57,530.175827 | 15.37.. |
| 58 | 151,348.908966 | Gamma72 laminated normal rank 14 | 115,060.351653 | 16.46.. |
| 59 | 262,144.000000 | Gamma72 laminated normal rank 13 | 244,079.864699 | 17.75.. |
| 60 | 524,288.000000 | Gamma72 laminated normal rank 12 | 773,597.146771 | 19.04.. |
| 61 | 1,048,576.000000 | Gamma72 laminated normal rank 11 | 2,043,183.683386 | 20.54.. |
| 62 | 2,421,582.543463 | Gamma72 laminated normal rank 10 | 5,188,971.717320 | 22.13.. |
| 63 | 6,552,970.479141 | 64D averaged ball-window projection | 13,105,940.958282 | 23.92.. |
| 64 | 43,046,721.000000 | published 64D input | 19,144,465.698823 | — |
| 65 | 37,279,553.935621 | Ne64 discriminant-count extension (coset minimum≥8) | 40,837,705.248540 | — |
| 66 | 77,490,641.390820 | Gamma72 laminated normal rank 6 | 82,897,968.185303 | — |
| 67 | 189,812,531.248503 | Gamma72 laminated normal rank 5 | 247,629,167.976493 | — |
| 68 | 536,870,912.000000 | Gamma72 laminated normal rank 4 | 812,141,667.469561 | — |
| 69 | 1,518,500,249.988025 | Gamma72 laminated normal rank 3 | 2,342,574,740.001863 | — |
| 70 | 4,959,401,049.012478 | Gamma72 laminated normal rank 2 | 7,120,856,007.299829 | — |
| 71 | 17,179,869,184.000000 | Gamma72 laminated normal rank 1 | 19,663,555,319.937195 | — |

## Consequences

- The present Ne64 bare sections in57,58,59 are below explicit propagated bounds. The57D comparison follows from the Schütt56 torus extension; the58D and59D comparisons follow from newly certified primitive Γ72 sections.
- The63D antipode is below the older printed Vardy1999 lower bound. Its density is therefore not a record candidate.
- The older Γ72 65D antipode is below the independent Ne64 discriminant-count extension. The candidate column now uses the separately certified explicit Ne65 lattice with center density3^17/√10; this new construction is not counted as a prior baseline.
- Candidates in60–62 and65–71 exceed these audited finite families. This statement does not establish current SOTA, because these families and the available published tables are not exhaustive.

## Families included

1. Direct products of the published56D and64D packings with listed small-dimensional packings, including the nonlattice10D,11D and13D examples.
2. Greedy separated offsets on the parent torus, combined with a small-dimensional lattice colored modulo2. The low-dimensional lattice minimum, not a nonlattice density, is used for this family.
3. The analogous A_k coloring with k+1 colors. Its same-color squared minimum ratio is4 for k=1,3 for k=2,and at least2 otherwise.
4. Averaged projection through a ball window from64D and72D.
5. Primitive Γ72 sections along explicitly embedded laminated normals in its √2Leech block.
6. Direct products of the resulting bounds with listed small-dimensional packings.
7. The separately proved Ne64 extension with a3-torsion dual coset of minimum at least8.

## Formulas and proof scope

For a lattice L normalized to minimum4 and center densityδ, greedy torus selection gives s offsets at pairwise distance at leastR whenever (s−1)δ V_n R^n<1. If the vertical min1 latticeK has an s-coloring with same-color squared minimum at leastλ, scale it by t with t²=max(4−R²,4/λ). The packing consisting of L+h_color over each vertical point has center density δ/[t^k covol(K)]. The certificate chooses rational R² and checks the union-volume inequality using π<355/113. Full derivation and an independent57/58 check are in `../code_routes/TORUS_LAYER_BENCHMARKS.md`.

An averaged k-dimensional ball-window projection of an N-dimensional packing gives δ_(N−k)≥δ_N V_k(k/N)^(k/2)((N−k)/N)^((N−k)/2). These values are evaluated below using the certified lower boundπ>333/106. The formulas are lower bounds, not upper limitations on what projection can achieve.

For a primitive normal GramK inside unimodular Γ72, the orthogonal section has determinant detK and inherited minimum8. Thus δ²=2^(72−rankK)/detK. The rank13 and14 determinants are2^23 and3·2^22, yielding δ59=262144 and δ58=2^18/√3. Source identity, Smith saturation and explicit complete section-period bases are checked in `verify_gamma72_leech_sections.py`.

The Vardy values are kept as the printed rounded log₂ lower bounds. No inferred exact Mordell–Weil determinant sequence is used. The catalogue MW64 Gram has determinant2^52, but the attempted identification of its sections with those printed values has not been established.

## Replay

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/verify_gamma72_leech_sections.py
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/audits/propagated_benchmarks.py
```
