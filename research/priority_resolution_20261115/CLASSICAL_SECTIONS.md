# Stronger classical comparisons from the BCS88 construction

These are **derived classical constructions**, not dimension-specific
records printed in BCS1982. They establish that our previous 86D/87D
results cannot be records, independently of the recent 88D preprint.

## The 88D parent

Apply Bos–Conway–Sloane's Theorems 1–2 and equation (15) to E8, with
`m=8`, `n=11`, `a=3`, and nested codes over GF(16) with parameters

```
[11,10,2], [11,8,4], [11,4,8].
```

Use ordinary Reed–Solomon evaluation at eleven distinct field elements;
the spaces of polynomials of degrees below 10, 8 and 4 are nested.
Their distances meet the BCS requirements 2, 4 and 8. The
[finite-field checker](verify_classical_sections.py) uses the irreducible
polynomial `x^4+x+1`, checks all 506 maximal minors, and checks constant
words, which make every coordinate projection surjective.

For E8 the parent minimum is 2, covolume 1 and center density 1/16.
Equation (15) therefore gives an 88D lattice of minimum 2 with

```
delta88 = (1/16)^11 * 2^(4*(10+8+4)) = 2^44.
```

Scale all lengths by 2. The resulting lattice Λ88 has minimum 8 and
covolume 1; its center density is unchanged. Integrality is not needed
for this classical lattice to be a valid sphere packing.

## Dual coordinate-block normals

Let D be the E8 norm-doubler and `T=D^-1`; thus T scales lengths by
`1/sqrt(2)`. The coordinate projection of the unscaled BCS lattice onto
any E8 block is exactly `T^3 E8`. Surjectivity of each code projection
supplies every quotient digit, while the base E8 supplies the integral
part. After scaling by 2 this becomes `2T^3 E8`, a rotated copy of
`E8/sqrt(2)`.

For any lattice L and coordinate subspace U,
`L* intersect U = (projection_U L)*`. Thus the dual intersection of
Λ88 with one 8D coordinate block is a rotated copy of `sqrt(2)*E8`.
Primitive root sublattices A1, A2, A3 and D4 inside E8 provide dual
normal Grams `K=2Q`, with determinants 4, 12, 32 and 64.

The [matrix reconstruction](construct_bcs88.py) makes this argument
concrete: it constructs the GF(16) code lifts, the complete rank-88 HNF
basis, the rational Gram of determinant 1, and each primitive dual
normal. HNF verifies saturation and computes the full section basis,
whose Gram determinant equals det(K). The saved Gram numerator and
denominator are in [bcs88_data.npz](bcs88_data.npz).

For a primitive dual normal, the perpendicular section has covolume
`covol(Λ88)*sqrt(det(K)) = sqrt(det(K))`. Its minimum is at least 8.
It attains 8 in an untouched coordinate block, which contains a copy of
`2E8`. Therefore its exact center density is

```
delta_(88-k) = 2^((88-k)/2)/sqrt(det(K)).
```

| Dimension | Normal | Classical center density | Consequence for our result |
| --- | --- | ---: | --- |
| 87 | sqrt(2) A1 | `2^(85/2)` | Above our 88D-parent antipode |
| 86 | sqrt(2) A2 | `2^42/sqrt(3)` | Above both our L80 counting bound and our 88D-parent antipode |
| 85 | sqrt(2) A3 | `2^40` | Our bound is higher; ratio squared `16384/9375` |
| 84 | sqrt(2) D4 | `2^39` | Our bound is higher; ratio exactly `64/25=2.56` |

The [exact report](classical_bcs88_sections.json) compares the complete
rational density squares. The stronger external gains are **156% in
84D** and **32.197831550546…% in 85D**. Chen's 2013 theorem gives weaker
numbers in both dimensions, so it is retained as historical attribution
instead of the leading README comparison.
