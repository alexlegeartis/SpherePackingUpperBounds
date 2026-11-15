# Allow one exceptional coset, then adapt the vertical lattice

This gives stronger lattice **existence bounds in 77D and 85D**. Horizontal
representatives are not explicit. The parent minimum and theta-series
theorems are the same published inputs as the preceding counting proof. All norm
bounds below are squared Euclidean norms.

For an even unimodular parent L, put q(v)=v²/2 modulo 2 on L/2L. In either
possible quadratic type, each nonzero type has population at least
`P_n=2^(n-1)-2^(n/2-1)-1`. An abstract five-dimensional quadratic model
embeds into five hyperbolic planes. Under a uniformly random orthogonal
image of that model, each nonzero class is uniform in its type.

If A_m counts parent vectors of squared norm m, at most
`(sum A_m)/2` binary classes contain such vectors: v and -v have the
same class. This is an upper bound, so collisions only help. The exact
shell counts are reconstructed from the extremal theta series.

For a vertical Gram M and binary representatives b_i, the full lattice
has generator blocks `(L,0)` and `(b_i/2, vertical_i)`. A vector with
coefficient parity a has horizontal part in the corresponding half-coset.
Its minimum follows by adding the horizontal and vertical bounds.
The full Gram determinant is det(M), so for minimum 8 its center-density
square in n+r dimensions is `2^(n+r)/det(M)`.

## 77D: move the exception into a longer D5 parity class

Use a totally singular five-space in Gamma72/2Gamma72. For its 31 nonzero
classes define the integer cost

```
Z = number with half-coset norm <5 + number with half-coset norm <4.
E[Z] <= 31 * ((A8+A12+A16)/2 + (A8+A12)/2) / P72
     = 1.398398587398988… <2.
```

Some five-space therefore has Z<=1: every nonzero half-coset has norm
at least 4, and at most one has norm below 5. This argument uses an
integer-valued cost; it does not require every bad class to disappear.

Choose `M=(3/2)*Gram(D5)`. Its minimum is 3. The nonzero binary class
represented physically by `2e1` has minimum 6. Identify the exceptional
horizontal class, if any, with that longer vertical class. GL(5,2) can
send any nonzero vector to it; the verifier checks explicit linear
bijections for all 31 possibilities. The horizontal model is totally
singular, so this change preserves its quadratic form.

Ordinary nonzero classes give `5+3=8`; the exception gives at least
`4+6>=8`; nonzero all-even vertical vectors have norm at least 12.
Thus

```
det(M)=243/8,
delta77²=2^80/243.
```

The gain over the preceding 77D bound is exactly `sqrt(4/3)-1`, or
**15.470053837925…%**.

## 85D: absorb one shallow nonsingular class in a fitted metric

Use q(a)=a*Gram(D5)*a^T/2 modulo 2, with 11 nonzero singular and 20
nonsingular classes. In L80 use the integer cost

```
Z = 2 * number of singular half-cosets with norm <5
    + number of nonsingular half-cosets with norm <11/2
    + number of nonsingular half-cosets with norm <9/2.
E[Z] <= [22*(A8+A12+A16)/2
         +20*((A10+A14+A18)/2 +(A10+A14)/2)]/P80
     = 1.139074004054923… <2.
```

Consequently every singular class has norm at least 5, every nonsingular
class has norm at least 9/2, and at most one nonsingular class lies below
11/2. All 20 nonsingular D5 classes are represented by its 40 roots
modulo sign. Signed coordinate permutations act transitively on them.
The verifier checks an integral unimodular D5 isometry for each class,
so the exception can be assigned to coefficient parity `(1,0,0,0,0)`.

Use the rational vertical Gram

```
M = [[ 7/2, -7/4,  1/2, -1/2,    0],
     [-7/4,  5/2, -5/4,    0,    0],
     [ 1/2, -5/4,  5/2, -5/4, -5/4],
     [-1/2,    0, -5/4,  5/2,    0],
     [   0,    0, -5/4,    0,  5/2]].
```

Its singular nonzero classes have norm at least 3, its nonsingular
classes at least 5/2, and the exceptional class at least 7/2. All-even
nonzero vectors have norm at least 8. Adding the horizontal bounds gives
minimum 8 in every case. These are exhaustive inequalities: positivity
and the exact inverse give coefficient bounds
`|z_i|² <= 8*(M^-1)_ii` for any potential violator, and every vector
in the resulting box is checked with exact integer arithmetic.

```
det(M)=3875/256,
delta85²=2^93/3875.
```

The exact ratio squared to our preceding bound is `75/62`, giving
**9.985336266014…%**. The metric was found by numerical determinant
descent; its proof uses only the rational matrix and exhaustive bounds.

[Exact verifier](verify_adaptive.py) · [Results](adaptive_results.json).
