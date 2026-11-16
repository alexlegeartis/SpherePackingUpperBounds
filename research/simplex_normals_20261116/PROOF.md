# Correlated dual normals improve 86D–89D

These four explicit periodic packings strengthen already registered
best-known lower bounds. All numbers below are **center densities**.
The main README uses occupied fractions instead.

| Dimension | New center density | Gain over commit 52d3fea | Named external comparison | Gain over that comparison |
| --- | ---: | ---: | --- | ---: |
| 86 | 4,208,042,413,315.359956… | 19.124698% | BCS88 perpendicular to sqrt(2) A2 | 65.722287% |
| 87 | 16,512,868,694,273.547878… | 35.346663% | BCS88 perpendicular to sqrt(2) A1 | 165.489721% |
| 88 | 61,981,923,829,060.643589… | 22.853263% | Bos–Conway–Sloane's 88D construction | 252.326446% |
| 89 | 227,360,271,990,343.399570… | 11.466525% | MW92 section from the full Rogers theorem | >218.140269% |

The 86D, 87D and 89D comparison numbers are **derived guarantees**,
not numerical records printed by those authors. They are the external
comparisons already established in the [preceding source audit](../mw92_sections_20261115/PROOF.md#external-comparisons-including-the-full-rogers-consequence).
For 89D the gain uses the upper endpoint of the full Rogers-derived
interval, not a weaker approximation. No stronger published guarantee
was located in the carried-forward review and the additional targeted
searches in [sources.json](sources.json).

## Why the normals matter

The preceding search primarily favored small normal determinants before
finding translate clusters. This continuation searches for correlated
dual normals whose projection labels form a nearly regular simplex.
That reduces the maximum distance removed by projection.

For example, the 86D squared projection loss falls from
`4743947052913/6874337093334 ≈0.69009` to
`35821491265610/100789271925211 ≈0.35541`. Seven translates now beat the
earlier eleven-translate construction by 19.12%. Shortness of individual
normal vectors alone does not determine the resulting density.

Discovery generated a 155,664-ray orbit of a dual vector of norm 244/47,
then searched signed subsets with negative correlations. Every new normal
has norm 244/47. A finite floating search selected the tuples, but neither
its completeness nor its stopping rules enter the proof. The supplied
[search.py](search.py) regenerates this discovery from an integer seed
and the previously verified unitary isometries.

## Source-bound parent and complete packing certificate

The parent is Shioda's rank-92 Mordell–Weil lattice, minimum 16 and
determinant 47^14. These are inherited published inputs from
Schütt–Shioda, *Mordell-Weil Lattices*, Theorem 13.35, p.357, at p=47.
The [existing reconstruction](../mw92_sections_20261115/PROOF.md) supplies
180 polynomial points on `y²=x³+1+t^48` over F_(47²)(t), all 16,110
pairwise heights, and a full integer spanning matrix. The new auditor
repeats every curve, height, determinant and span check.

Let n be the target, k=92−n, and B be its integer k-by-92 normal-label
matrix in [sections.json](sections.json). Put `K=B G^-1 B^T`. The physical
dual normals have coefficient rows `B G^-1`. Use the k+1 labels

```
S={0,e1,...,ek}.
```

An exact unimodular HNF transformation T satisfies `T B^T=[I_k;0]`.
Thus the first k rows A give integer lifts of every label, while the
remaining rows N are the complete integer kernel. The section period
Gram is `Q=N G N^T`, and `det(Q)=47^14 det(K)`.

For each s in S the projected shift in parent coordinates is

```
h_s=s A-s K^-1 B G^-1.
```

The packing is the union of `Z^n N+h_s`. Its maximum squared projection
loss is the exactly computed rational number

```
beta=max_(s,t in S) (s-t) K^-1 (s-t)^T.
```

Every preprojection difference is a parent vector. Distinct labels give
a nonzero vector, so every nonzero packing distance has squared length
at least `16-beta`. The same-label distances inherit 16. Since beta<16,
the translates are distinct. Therefore

```
delta_n²=(k+1)²*((16-beta)/4)^n / (47^14 det(K)).
```

All labels, integer lifts, full kernels, determinants and projected
shifts are checked exactly in [verify.py](verify.py). No minimum
enumeration or unclassified shortest-vector shell is used. This proof
does not claim that the guaranteed minimum is attained.

## A separate attack on the nine missing dimensions

[verify_craig.py](verify_craig.py) checks the generalized Craig/binary-code
route against the named published guarantees in 48D, 54D–56D, 63D–64D,
72D, 80D and 92D. This supplies **no new packing SOTA**.

Chen (2013), Theorems 2.3–2.4, gives the guaranteed center-density square

```
2^(2k-n) * m^n / (ell^(2m-2)*(n+1)),
```

using an even binary [n+1,k,≥8m] code and an odd prime ell≥n+1. The
shorter-code version in Theorem 2.4 is included by parity extension.
The smallest admissible prime maximizes the formula, so checking it
covers all larger primes. The verifier checks all m in the published
ranges and also the adjacent boundary value. If 8m>n+1, a nonzero code
is impossible and only dimension zero needs consideration.

Hamming and Griesmer bounds eliminate most parameters. Twelve saved
rational Delsarte polynomials handle the remaining cases. For normalized
binary Krawtchouk polynomials `P_j(w)=K_j(w)/binomial(N,j)`, each certificate
has the form

```
f(w)=1+sum_(j>=1) a_j P_j(w), a_j>=0,
f(w)<=0 at every permitted nonzero even distance.
```

The character-sum identity makes every P_j positive definite on the
binary cube. Summing f over codeword pairs gives `M²<=M f(0)`, hence
`M<=f(0)=1+sum a_j`. If r coordinates were shortened first, a linear
code loses at most r dimensions, giving `2^k<=2^r f(0)`. The checker
verifies every coefficient, grid sign and resulting density comparison
with rational arithmetic. No numerical LP status is trusted.

The resulting **ceilings on this guarantee formula** lie below the
comparison guarantees in eight of the nine target dimensions. They are
not upper bounds on actual Craig-lattice densities; a separate stronger
minimum proof or a different construction could improve the guarantee.

After optimizing ell, the only surviving 92D parameter case is
`m=3, ell=97` with an
**even binary [93,43,24] code**. Such a code would improve MW92 by
24.9189116…% through Chen's formula. **No such code is constructed here.**
The consulted [Grassl table](https://www.codetables.de/BKLC/BKLC.php?q=2&n=93&k=43)
provides a distance-18 construction and an upper bound of 24. Our exact
code LP also permits dimension 43, so it leaves this case open.

A 220-start bounded metric search in 85D retained the preceding metric;
a finite three-vector orbit sieve in MW92 found no vector shorter than
244/47. These finite failures are not nonexistence or minimum proofs.

## Reproduce

```sh
python research/simplex_normals_20261116/verify.py
python research/simplex_normals_20261116/verify.py --export NEW_DIRECTORY
python research/simplex_normals_20261116/verify_craig.py

# Optional complete discovery replay; keep BLAS thread counts modest.
OPENBLAS_NUM_THREADS=2 python research/simplex_normals_20261116/search.py --output NEW_SEARCH_DIRECTORY

# Optional exact code-LP reconstruction, independent of a floating LP solver.
python research/simplex_normals_20261116/produce_code_bound.py --length 93 --distance 24 --output NEW_CODE_CERTIFICATE.json
```

The packing auditor needs NumPy and python-flint. The Craig certificate
checker uses only the standard library. Both expose their exact results
in [results.json](results.json) and [craig_results.json](craig_results.json).
