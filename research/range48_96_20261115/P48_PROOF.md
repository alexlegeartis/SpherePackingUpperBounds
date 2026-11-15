# Compatible half-cosets above a classical 48D extremal lattice

**Earlier checkpoint:** the [subsequent deeper-coset proof](../covering_and_deep_cosets_20261115/DEEPER_PROOF.md)
gives stronger explicit 49D–52D lattices. The counting proof and
certificates below remain valid for their original, weaker bounds.

This extends the repo's existing binary quadratic-space counting method
to a different parent and dimensions **49–54**. It gives useful bounds in
49–52D; the 53D/54D variants are below the stronger classical comparisons
in the range audit. The existence proof below is now accompanied by
**explicit, fully enumerated shifts in 49–52D**. The weaker 53D/54D
variants remain existential.

## Parent and short-vector counts

Let L be a published even unimodular rank-48 lattice of minimum 6, such
as Leech–Sloane's P48q. This classical minimum theorem is inherited.
BCS1982 §5 gives the P48q input; the catalogue also lists P48p, P48q,
P48n and P48m. The conclusion uses no result above dimension 96.

Its theta series is the weight-24 modular form

```
Theta_L = E4^6 - 1440 E4^3 Delta + 125280 Delta^2.
```

The coefficients are determined by the constant term 1 and vanishing at
norms 2 and 4, using the standard basis of weight-24 modular forms.
The first relevant shell populations are

| Squared norm | Number of vectors |
| --- | ---: |
| 6 | 52,416,000 |
| 8 | 39,007,332,000 |
| 10 | 6,609,020,221,440 |
| 12 | 437,824,977,408,000 |

The checker independently expands Delta using both its eta product and
`(E4^3-E6^2)/1728`; it solves the two vanishing equations and checks the
norm-6 coefficient against the published kissing number.

## Good binary subspaces

On `V=L/2L`, use the nondegenerate quadratic form
`q(x)=||x||^2/2 mod 2`. For a nonzero singular class, exclude a norm-8
representative. Its minimum representative norm is then at least 12,
so its half-coset has squared depth at least **3**.
For a nonsingular class, exclude norms 6 and 10. Its half-coset then has
squared depth at least **7/2**.

Counting opposite vectors only once gives bad-class upper bounds

```
B0 = 19,503,666,000,
B1 = 3,304,536,318,720.
```

For either possible Arf invariant of V, the number of classes of each
type (excluding zero for the singular type) is at least
`D=2^47-2^23-1`. Finite orthogonal-group averaging gives the sufficient
condition `(z*B0+o*B1)/D < 1`, where z,o count the nonzero singular and
nonsingular vectors of the desired subspace.

Take the root-lattice quadratic form `q_R(z)=zQz^T/2 mod 2` for
`R=A1,A2,A3,D4,D5,E6`. Each embeds in r hyperbolic planes, available
in either 48D quadratic type. Explicitly, with hyperbolic bases e_i,f_i,
use `w_i=e_i+f_i+sum_(j<i) (Q_ij mod 2) f_j`. The e-coordinates prove
independence, and the diagonal and polar pairings agree with q_R.
The checker verifies this identity on every binary coefficient vector.

| r | z | o | Expected bad classes, upper bound (approx.) |
| --- | ---: | ---: | ---: |
| 1 | 0 | 1 | 0.023481 |
| 2 | 0 | 3 | 0.070441 |
| 3 | 1 | 6 | 0.141020 |
| 4 | 3 | 12 | 0.282178 |
| 5 | 11 | 20 | 0.471128 |
| 6 | 27 | 36 | 0.849027 |

All six exact rational expectations are below 1, so a good subspace
exists. Ranks 7/8 fail this particular sufficient inequality; this is
not a nonexistence result.

## Glue the vertical metric

Set `M=sqrt(5/4)*R`, of minimum 5/2, with Gram `G=(5/4)Q`.
Choose representatives b_i in L for a basis of a good subspace. Generate
a lattice in `R^48 ⊕ R^r` by L in the first factor and the rows
`(b_i/2, v_i)`, where v_i are the basis vectors of M.

For an integer layer vector z:

* If z is zero, the parent minimum is 6.
* If z is nonzero and even, its vertical squared norm is at least
  `4*(5/2)=10`.
* If q_R(z)=1, the horizontal norm is at least 7/2 and the vertical
  norm is at least 5/2, summing to 6.
* If z is nonzero modulo 2 and q_R(z)=0, the horizontal norm is at
  least 3. The positive even integer `zQz^T` is divisible by 4, so the
  vertical norm is at least 5, summing to at least 8.

Thus the resulting lattice has minimum **exactly 6**, attained in L.
Its Gram determinant is `det(G)` by the block triangular basis.
Consequently

```
delta_(48+r)^2 = (3/2)^48 * (6/5)^r / det(Q).
```

The new useful center densities are

| Dimension | Guaranteed center density (approx.) |
| --- | ---: |
| 49 | 13,039.6472366809 |
| 50 | 11,663.0150495343 |
| 51 | 11,064.5075824050 |
| 52 | 12,120.5607811403 |

They exceed the derived product/section comparisons identified in
[TABLE.md](TABLE.md). No dimension-specific previous worldwide record
was identified for these four dimensions, so they are not labelled
confirmed SOTAs. The 53D/54D variants are below the Leterrier54 input
and its 53D Blichfeldt section, and are excluded from the improvements.

Replay `python research/range48_96_20261115/verify_range.py`.
The exact expectations, Grams, density squares and unsuccessful ranks
are in [results.json](results.json), under `p48_counting`.

## Explicit 49–52D lattices

The [certificate](explicit_p48/certificate.npz) contains the archived
P48q Gram G, a unimodular basis change T, the reduced Gram `Q=TGT^T`,
four original and reduced binary residue rows, and the vertical D4 root
Gram. Its leading principal blocks are A1, A2, A3 and D4, so the same
four compatible shifts give all four constructions by taking prefixes.
The original catalogue Gram is independently parsed from the included
[source snapshot](explicit_p48/P48q.html).

For each of the **15 nonzero binary combinations** of the four shifts,
a complete congruence-class enumeration finds **no vector of squared
norm at most 10**. The checker verifies the quadratic type of every
combination. Thus each nonsingular half-coset has minimum at least 7/2,
and each nonzero singular half-coset has minimum at least 3, exactly as
required above. No unenumerated sampled class is used in the certificate.

For r=1,2,3,4, let B be the first r original residue rows and let R be
the corresponding leading root Gram. An explicit full Gram is

```
[[ G,            G B^T/2                       ],
 [ B G/2,        B G B^T/4 + (5/4) R           ]].
```

The auditor reconstructs it and verifies its determinant independently
of the counting formula. The preceding layer argument proves minimum
exactly 6. Original residue rows, vertical Grams, determinants and exact
density squares are also supplied in human-readable
[explicit_results.json](explicit_results.json).

The enumerator uses outward-rounded IEEE-754 binary64 bounds, exact
dyadic enclosures of a rational LDL decomposition, and exact integer
norm checks at leaves. Compile with `-fno-fast-math -ffp-contract=off`.
Its safety guards abort on unsupported ranges. This is a complete
enumeration under those documented arithmetic/compiler assumptions,
not a formal proof assistant certificate. It adapts the repository's
earlier Ne64 coset enumerator to dimension 48, modulus 2 and cutoff 10.
A deliberately short class is also enumerated and **rejected with an
exact short-vector witness**, providing a negative control.

Quick audit of exact inputs and the archived completed reports:

```sh
python research/range48_96_20261115/verify_explicit.py
```

Fresh complete replay, compiling GNU C++ and checking all fifteen classes
plus the negative control (NumPy, SymPy and python-flint required):

```sh
python research/range48_96_20261115/verify_explicit.py --output /tmp/p48_replay
```

The output directory must be new. The archived
[manifest](explicit_p48/replay/manifest.json) records source, data,
input and binary hashes and the compiler flags. The quick audit does
not rerun the enumerations; the `--output` command does. Both retain
the published minimum-6 theorem of P48q as an external input.
