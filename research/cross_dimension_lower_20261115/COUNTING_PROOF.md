# Compatible deep half-cosets above Gamma72

This proves existence of lattices in dimensions 73–78 and 81–85. It does **not** supply
explicit horizontal coset representatives, and establishes no world-record
or novelty claim. The numerical comparisons below are derived orthogonal
products, not published dimension-specific records.

The parent is Gabriele Nebe's even unimodular lattice `L = Gamma72`, of
minimum squared norm 8 and Gram determinant 1: [*An even unimodular
72-dimensional lattice of minimum 8*, arXiv:1008.2862v3, Theorem 2.3](https://arxiv.org/pdf/1008.2862).
The published minimum theorem is an input, not re-enumerated here.

## 1. Count bad half-cosets exactly

Put `V = L/2L`. This is a 72-dimensional binary quadratic space with
`q([w]) = (w,w)/2 mod 2` and polar form `(w,v) mod 2`. Evenness makes q
well-defined, and unimodularity makes the polar form nondegenerate.
The half-coset associated to `[w]` is `w/2 + L`.

The standard theta-series modularity theorem for even unimodular lattices
gives a weight-36 modular form. The four-dimensional space of these forms
has basis `E4^9, E4^6 Delta, E4^3 Delta^2, Delta^3`. The constant term 1
and the absence of norms 2, 4, 6 determine

```
Theta_L = E4^9 - 2160 E4^6 Delta + 965520 E4^3 Delta^2 - 27302400 Delta^3.
```

Here the coefficient of `q^j` counts vectors of squared norm `2j`.
`verify_counting.py` calculates E4 from divisor sums, Delta from its eta
product, and independently checks Delta through `(E4^3-E6^2)/1728`.
It solves the triangular vanishing equations and obtains:

| Squared norm | Number of vectors |
| --- | ---: |
| 8 | 6,218,175,600 |
| 10 | 15,281,788,354,560 |
| 12 | 9,026,867,482,214,400 |
| 14 | 1,989,179,450,818,560,000 |
| 16 | 213,006,159,759,990,870,000 |

The norm-8 count agrees with Nebe's independently computed kissing number
(paper p. 2). Theta modularity and the classical modular-form basis are
mathematical inputs; agreement of arithmetic formulas does not independently
prove those theorems.

For a nonzero singular class (`q=0`), exclude norms 8, 12, 16; every
remaining representative has squared norm at least 20, hence its half-coset
has minimum at least **5**. For a nonsingular class (`q=1`), exclude norms
10 and 14; its remaining norms are at least 18, hence half-coset minimum
at least **9/2**. Smaller positive norms vanish by the parent minimum.

Antipodes `w,-w` belong to the same binary class. Therefore upper bounds
for the two bad-class counts are

```
B0 = (a8+a12+a16)/2 = 106507593316845630000,
B1 = (a10+a14)/2    =    994597366303457280.
```

No listed vector belongs to the zero class: `2L` has minimum 32. Multiple
antipodal pairs in one class only make these upper bounds conservative.

The classification of nondegenerate quadratic spaces over F2 gives
`2^71 +/- 2^35 - 1` nonzero singular vectors and `2^71 -/+ 2^35`
nonsingular vectors. We use the common **lower** bound

```
D = 2^71 - 2^35 - 1 = 2361183241400462868479
```

for both populations, so the proof does not assume an Arf invariant.

## 2. Make the shifts compatible using orthogonal-group averaging

Consider the following subspaces W of V, specified by their restricted q.
The two ambient quadratic types both have Witt index at least 35, so the
models below embed using six hyperbolic planes. In hyperbolic coordinates
write `q(sum x_i e_i + y_i f_i) = sum x_i y_i`.

* Ranks 1–4: W is spanned by `e1,...,er`, and q vanishes identically.
* Rank 5: take `e1,...,e4,e5+f5`, giving `q(a)=a5`.
* Rank 6: take `e1,...,e4,e5+f5,e5+e6+f6`, giving
  `q(a)=a5+a6+a5*a6`. The last two coordinates form an anisotropic plane.

`verify_counting.py` enumerates the binary combinations of these small
models and checks independence and their restricted quadratic forms.
They contain the following numbers `(z,o)` of nonzero singular and
nonsingular vectors:

| Rank r | z | o | Upper bound `(z B0 + o B1)/D` |
| --- | ---: | ---: | ---: |
| 1 | 1 | 0 | 0.045108 |
| 2 | 3 | 0 | 0.135324 |
| 3 | 7 | 0 | 0.315755 |
| 4 | 15 | 0 | 0.676616 |
| 5 | 15 | 16 | 0.683356 |
| 6 | 15 | 48 | 0.696835 |

Apply a uniformly random element of the finite orthogonal group of V to
one fixed model. This group is transitive on the nonzero singular vectors
and separately on the nonsingular vectors. Linearity of expectation bounds
the number of bad vectors in the image W by `(z B0 + o B1)/D < 1`.
The number is a nonnegative integer, so **some image W contains no bad
nonzero class**. The script checks every strict inequality as a rational
number. This is averaging over isometries of the finite quadratic space;
they need not lift to automorphisms of L.

Thus all nonzero combinations of a basis of W give compatible half-cosets,
with depth 5 for singular combinations and 9/2 for nonsingular ones.
Neither an exhaustive CVP computation nor explicit representatives are
needed for this existence conclusion.

## 3. Match the vertical metric to the two depths

Let `h_i = w_i/2` be representatives of a basis of this good W. For an
r-dimensional vertical lattice M with basis `b_i`, form the lattice

```
Lambda = < (L,0), (h_1,b_1),...,(h_r,b_r) >.
```

For r at most 4, take M to be `sqrt(3/2)` times A1, A2, A3, D4,
respectively (roots of squared norm 2). Its minimum is 3.

For r=5, start with `sqrt(3/2) D4` in four orthogonal coordinates and add
the layer generator `(sqrt(3/2) c, sqrt(2))`, where
`c=(1/2,1/2,1/2,1/2)` is a D4 deep hole. The value q is the parity of
the fifth coefficient. If that parity is odd, the D4 component has squared
norm at least 3/2 and the last coordinate contributes at least 2, giving
**7/2**. If it is even, any nonzero vector has squared norm at least **3**:
zero last coordinate leaves a D4 vector, and a nonzero even last coordinate
already contributes at least 8. Here `2c` lies in D4.

For r=6 use the three nonzero glue classes of `D4*/D4`. Take
`c1=(1,0,0,0)` and `c2=(1/2,1/2,1/2,1/2)` for two of them, and attach
two vertical generators whose projected Gram is `[[2,-1],[-1,2]]`.
All three nonzero D4 glue classes have minimum 1. The parity of the last
two coefficients is zero exactly when q is zero. A nonsingular vector
has squared norm at least `3/2+2=7/2`. A nonzero singular vector has
norm at least 3; a nonzero even pair in the last two coordinates alone
contributes at least 8. These facts also show `min(M)=3`.

For completeness, `D4={z in Z^4: sum(z_i) even}` has minimum 2. Its odd
integer glue class has minimum 1, and its half-integer glue classes have
minimum 1 because all four coordinates have absolute value at least 1/2.
This proves the glue-class bounds without a numerical nearest-point search.
The A2 projection has norm `2(a*a-a*b+b*b)`, positive and at least 2
on every nonzero integer pair. The script constructs the rational Grams
and checks positive definiteness and their determinants exactly.

For any vector of Lambda with vertical coefficient vector a:

* `a=0`: a nonzero vector belongs to L and has squared norm at least 8.
* `a mod 2` is nonzero: the good half-coset contributes at least 5 or
  9/2 horizontally; M contributes at least 3 or 7/2, respectively. The
  total is at least 8.
* `a mod 2=0`, `a!=0`: its vertical component is a nonzero vector of 2M,
  with squared norm at least 12.

The embedded parent attains 8, so **min(Lambda)=8**. A block triangular
basis gives `covol(Lambda)=covol(M)`, since the parent is unimodular.
Consequently the **Gram** determinants in dimensions 73–78 are

```
3, 27/4, 27/2, 81/4, 81/2, 243/4,
```

and the center-density squares are `2^(72+r)/det(M)`. These lattices
are not asserted to be integral or unimodular.

## Comparison scope and remaining research

The comparison in `counting_results.json` is the construction
`Gamma72 x R_r`, with the small root lattice scaled to minimum 8. Its
center-density square is `2^72/(2^r det(R_r))`. These are derived
constructions using Nebe's parent and classical root lattices. No claim
is made that a source prints these 73–78D values, that they were the best
previous packings, or that this gluing method is historically new.

The mixed-depth refinement in rank 5 improves the uniform-depth-9/2
version, whose determinant is `7^5/2^8`, by the density factor
`sqrt(7^5/(2^7*81)) = 1.273...`. This is a comparison between variants
derived in this experiment, not an external record.

The next substantial step is constructive: find a good W and certify its
15/31/63 nonzero half-cosets efficiently. The counting proof by itself
does not provide their coordinates. Improving depth 5 for 31 singular
classes also requires a stronger argument: the present union bound for
a totally singular rank-5 space is `31 B0/D > 1`, so it proves nothing
in that case.

Replay without network or scientific dependencies:

```sh
python research/cross_dimension_lower_20261115/verify_counting.py
```

## 4. The 80D parent gives deeper nonsingular half-cosets

The same argument applies to the [Bachoc–Nebe L80](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html)
or M80, both even unimodular of minimum 8. The parent minimum is an
external theorem: C. Bachoc and G. Nebe, *Extremal lattices of minimum 8
related to the Mathieu group M22*, J. reine angew. Math. **494** (1998),
155–171. This is the same theorem inherited by the explicit antipodes.

Weight-40 theta modularity, in basis
`E4^10, E4^7 Delta, E4^4 Delta^2, E4 Delta^3`, now gives coefficients

```
1, -2400, 1360800, -103488000.
```

The norm-8 count is 1,250,172,000, agreeing with the L80 catalogue. There
are sufficiently few short vectors to make **singular depth 5** and
**nonsingular depth 11/2** compatible in suitable ranks up to 5. Thus B0
still counts norms 8,12,16, whereas B1 now counts norms 10,14,18. The
common population lower bound becomes `2^79-2^39-1`. Exact counts and
expectations are saved in `counting80_results.json`.

Choose the restricted quadratic form of W to be the one on `R/2R`,
`q_R(a) = a G_R a^T/2 mod 2`, for the following root lattices R. Their
actual forms embed in five hyperbolic planes: write the quadratic
polynomial as `a^T U a`, with U upper triangular, and embed a as
`(a,U a)` in a hyperbolic space with form `x.y`. This also handles
degenerate restricted forms. The code checks these embeddings on every
binary vector, not just the counts of singular and nonsingular vectors.

| Dimension | R | Nonzero singular z | Nonsingular o | Expected bad-class bound |
| --- | --- | ---: | ---: | ---: |
| 81 | A1 | 0 | 1 | 0.056324 |
| 82 | A2 | 0 | 3 | 0.168972 |
| 83 | A3 | 1 | 6 | 0.338513 |
| 84 | D4 | 3 | 12 | 0.677595 |
| 85 | **A5** | 15 | 16 | 0.909728 |

The root lattice **D5 does not work in the last row**: its actual restricted
form has 11 nonzero singular and 20 nonsingular vectors, giving expectation
upper bound `1.1327... > 1`. Equal dimensions do not make these finite
quadratic forms interchangeable. This failure is why A5 is used.

Set `M = sqrt(5/4) R`. For a nonsingular coefficient class the unscaled
root-lattice norm is at least 2, giving vertical norm at least 5/2;
the horizontal depth 11/2 supplies a total of 8. For a nonzero singular
coefficient class the unscaled norm is divisible by 4 and hence at least
4, giving vertical norm at least 5. With horizontal depth 5 the total
is at least 10. An all-even nonzero coefficient vector lies in 2M, of
minimum 10. The embedded parent still attains 8. The block determinant
argument is unchanged.

This proves minimum-8 lattices in dimensions 81–85 of Gram determinants

```
5/2, 75/16, 125/16, 625/64, 9375/512.
```

The explicit comparisons use `L80 x A1/A2/A3/D4/D5` rescaled to common
minimum 8. In particular the 85D comparator uses the denser D5 product,
even though our proof needs A5 vertically. These are derived products,
not a claim about published records. No horizontal representatives have
been found by this proof.

```sh
python research/cross_dimension_lower_20261115/verify80_counting.py
```
