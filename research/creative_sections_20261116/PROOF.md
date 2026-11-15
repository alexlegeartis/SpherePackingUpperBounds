# Stronger explicit bounds in 57D, 78D and 79D

These three certified constructions improve the repository's previous
bounds. They strengthen three already registered dimensions; the nine
remaining dimensions are still open here. All densities below are
**center densities**, not occupied fractions.

| Dimension | New center density | Gain over commit f021719 | Named external comparison | Gain over that comparison |
| --- | ---: | ---: | --- | ---: |
| 57 | 165,830.432487143… | 1.678302% | Shioda MW56 × minimum-matched A1 | 69.030851% |
| 78 | 113,833,877,612.492613… | 8.690819% | Bachoc–Nebe L80 section perpendicular to 2 A2 | 43.457189% |
| 79 | 331,702,761,803.955754… | 10.871173% | Bachoc–Nebe L80 section perpendicular to a minimum vector | 20.672762% |

The comparison numbers are **derived consequences** of the named
published constructions. They are not same-dimensional numerical records
printed by those authors. The main README continues to use these external
comparisons; the incremental column above separately measures our progress.

## 57D: retain one certified half-coset layer

The preceding [58D proof](../cm58_20261116/PROOF.md) supplies the explicit
MW56 Gram Q, determinant 29^8, minimum 10, and the binary row c with

```
min ||c/2+MW56||² >=13/2.
```

This inequality has two complete enumeration reports with exact LDL
enclosures and the arithmetic assumptions stated there. Retaining only
the first vertical generator gives the explicit lattice

```
L57 = < (MW56,0), (c/2,sqrt(7/2)) >.
```

Odd vertical coefficients have horizontal norm at least 13/2 and vertical
norm at least 7/2. Nonzero even coefficients have vertical norm at least
14. Zero coefficients inherit the parent minimum 10. A parent basis
vector attains 10. Thus the minimum is exactly 10, and

```
det(L57)=(7/2)*29^8,
delta57²=(5/2)^57 / ((7/2)*29^8).
```

The full rational Gram is [gram57.json](gram57.json). The external product
MW56×A1 has half the MW56 center density, so the exact density ratio is
`sqrt(20/7)`. This is a new extracted consequence of our certified coset;
it requires no new shortest-vector enumeration.

## 78D and 79D: antipodes from the reconstructed MW80

Schütt–Shioda, *Mordell-Weil Lattices*, Theorem 13.35, p.357, at p=41,
gives the parent of rank 80, minimum 14 and determinant 41^12. The
[source-bound reconstruction](../cm58_20261116/FRONTIER.md) checks 160
unitary polynomial points on `y²=x³+1+t^42`, every curve equation, all
12,720 pairwise heights by direct polynomial gcd, and a full integer span.
The parent facts are inherited published inputs.

The integer normal-label matrices B and complete label clusters S are in
[sections80.json](sections80.json). The 79D normal Gram is `[234/41]`,
with labels 0 and 1 and squared projected diameter 41/234. In 78D,

```
K = B G^-1 B^T = [[188/41,-88/41],[-88/41,178/41]],
S = {(0,0),(1,0),(0,1)},
beta=779/2572.
```

The target in the normal search was the density after projection, rather
than the shortest possible dual normal. In 79D the two-translate formula
is maximized at normal norm `80/14=40/7`; 234/41 is close to that target.
In 78D, correlated pairs make the three projected distances nearly equal.
The search used finite BKZ reductions and unitary-isometry orbits. None
of its stopping rules or completeness is needed for the proof.

For each certificate, exact unimodular HNF gives
`T B^T=[I;0]`. The first rows of T are integer lifts A, and the remaining
rows N form the complete section lattice. Its Gram is `N G N^T`, with
determinant `det(G)det(K)`. For each s in S the projected shift is

```
h_s=s A-s K^-1 B G^-1.
```

The packing is the union of these translates of the section lattice.
Every nonzero difference comes from a parent vector of norm at least 14,
with projected squared length at most beta removed. Hence the proved
minimum is at least `14-beta`. The saved parent norm-14 vector in each
certificate has a label difference attaining beta, proving equality
for both minima. Exactly,

```
delta_n²=|S|²*((14-beta)/4)^n / (41^12*det(K)).
```

The auditor checks all kernels, determinants and shifts rationally. Its
`--export NEW_DIR` option emits full integer period bases, period Grams
and rational shifts in period coordinates.

The named external parent is Bachoc–Nebe's published even unimodular
80D lattice of minimum 8. Its section normals have determinant 48 in
rank two and 8 in rank one. The corresponding inherited-minimum section
guarantees have center-density squares `2^74/3` and `2^76`. The existing
[exact L80 audit](../cross_dimension_lower_20261115/l80_results.json)
binds these numbers to actual normals. They do not assert an upper bound
on the unknown actual minimum of every section.

## Failed creative routes, with counterexamples

The 92D deformation is now [rejected by an exact vector](../cm58_20261116/FRONTIER.md#92d-exact-counterexample-rejects-the-deformation).
The vector has parent norm 16 and deformed norm 29/2. Its order-three
orbit also completes a proof of local lattice optimality of the parent.
The earlier conditional gain was never registered as a lower bound.

In 56D, three F4-linear checks define an index-64 sublattice that deletes
all 2,436 previously supplied shortest rays. Raising its minimum from 10
to 12 would give a factor `(6/5)^28/64>2.57`. However, an additional
norm-10 vector survives. The [negative certificate](rejected_code56.json)
contains the three binary masks and that exact integer witness. The
auditor regenerates the known shell, checks its deletion, and verifies
the surviving vector. This candidate is excluded from the density table.

## Run the checks

```sh
python research/creative_sections_20261116/verify.py
python research/creative_sections_20261116/verify.py --export NEW_DIRECTORY
python research/cm58_20261116/verify_frontier.py
python research/range48_96_20261115/verify_range.py
```

NumPy, python-flint and SymPy are required. Exact gains and their certified
rational intervals are in [results.json](results.json). Source locators,
review scope and search queries are in [sources.json](sources.json).
