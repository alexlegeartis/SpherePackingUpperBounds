# Periodic 93–96D packings above the classical MW92 lattice

Let Λ be Shioda's published rank-92 lattice with squared minimum 16 and
covolume `v=47^7`. Its center density is `delta92=2^92/47^7`.
This parent theorem is inherited from Schütt–Shioda, Theorem 13.35 at
`p=47`. No input has dimension above 96.

## Greedy offsets in the torus

The flat torus `T=R^92/Λ` has volume v. A torus ball of radius √d is the
image of a Euclidean ball, so its volume is at most

```
B(d) = pi^46 * d^46 / 46!.
```

If `(N-1)*B(d)<v`, choose N points successively. At the j-th choice the
union of j forbidden balls has volume at most `j*B(d)<v`, leaving a point
outside. Thus the N offsets have pairwise squared torus distance at
least d. This is a continuum existence argument, without a nearest-point
oracle or explicit coordinates. The strict inequality makes boundary
conventions immaterial.

## Color a small vertical lattice

Let R have dimension r, minimum 2 and Gram Q. Choose a sublattice P of
index N and minimum `2*lambda`. Color R by its N cosets modulo P.
Scale R to a lattice M with squared minimum `t=16-d`, using length
factor `sqrt(t/2)`. Suppose `lambda*t>=16`.

Assign the N torus offsets to the N colors, and put a copy of Λ above
each point of M, translated horizontally by its color's offset.
The packing is periodic with period `Λ ⊕ sqrt(t/2)*P` and N translates.
For differences of points:

* In the same layer, the squared distance is at least 16.
* In different layers of the same color, the vertical squared distance
  alone is at least `lambda*t>=16`.
* In different colors, horizontal and vertical squared distances are
  at least d and t respectively, summing to 16.

The parent attains 16 in every layer, so the packing's squared minimum
is exactly 16. The volume per point is
`v*N*covol(M)/N = v*covol(M)`; the N offsets are **not** an extra
density factor. With packing radius 2,

```
delta_(92+r)^2 = delta92^2 * (8/(16-d))^r / det(Q).
```

This is the same greedy-torus/colored-layer construction used earlier in
[the repo's torus-layer benchmarks](../novel_routes_20261113/code_routes/TORUS_LAYER_BENCHMARKS.md).
The improvement here comes from the MW92 input and the color periods.

## Exact periods and choices

Matrices act on row coefficients in the standard root basis.

| r | R | Index N | Period minimum ratio lambda | Chosen d | Center density |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | A1 | 2 | 4 | 10.292706 | 8.182565090782445e15 |
| 2 | A2 | 3 | 3 | 10.138774 | 7.702212346014765e15 |
| 3 | A3 | 6 | 3 | 9.938814 | 7.410425802173365e15 |
| 4 | D4 | 9 | 3 | 9.837781 | 8.236644396376130e15 |

For A1 use `C=[2]`. For A2 use rows `(1,-1),(1,2)`; the resulting Gram
is `3Q`. For A3 use rows `(6,0,0),(-1,1,0),(-2,0,1)`, the kernel of
`a+b+2c mod 6`. All vectors with norm 2 or 4 are excluded by direct
integer enumeration. This is exhaustive: the diagonal of `Q^-1` is
`(3/4,1,3/4)`, and `z_i^2 <= (Q^-1)_ii*zQz^T` confines every vector of
norm below 6 to the enumerated box. The vector `(-1,1,0)` has norm 6.

For D4 use the norm-3 quaternion action. In physical integer coordinates,
with D4 the even-coordinate-sum sublattice of Z4, the matrix is

```
A = [[1,-1,-1, 0],
     [1, 1, 0, 1],
     [1, 0, 1,-1],
     [0,-1, 1, 1]].
```

It satisfies `A*A^T=3I` and maps D4 into D4. Its coefficient matrix C in
the root basis is integral with determinant 9 and `CQC^T=3Q`, proving
the index and minimum exactly. All matrices are in the JSON report.

The [checker](verify_torus_extensions.py) derives an exact rational upper
bound for pi using Machin's formula, then certifies the strict forbidden
volume inequality and the same-color condition for every row. The chosen
d is the largest multiple of `10^-6` passing that particular rigorous
volume bound; this is not a claim of an optimal offset configuration.

## Published 96D comparison

BCS1982 equation (15), with P48q of center density `(3/2)^24`, `m=48`,
`n=2`, and `a=1`, gives

```
delta_BCS96 = (3/2)^48 * 2^24 = 3^48/2^24.
```

Its logarithm rounds to the Table 1 entry 52.078. Our 96D value is
`delta92 * 32/(6.162219)^2`, with exact rational squared density in
[the report](torus93_96_results.json). Squaring the two endpoint ratios
proves exactly that the percentage gain lies in

```
73.240722% < gain < 73.240723%.
```

The decimal occupied fractions are approximately `4.842118544896e-22`
for the new packing and `2.795023294827e-22` for BCS96.
This is a proved improvement over that named published construction;
worldwide priority has the literature scope stated in the README.
