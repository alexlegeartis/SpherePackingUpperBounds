# A further certified 55D improvement, and new attacks on 54D and 56D

The new 55D lattice has minimum **exactly 10** and center density
**91,333.867143595418…**. This is **0.064703986…%** above our preceding
91,274.808704… construction, and **0.966206016…%** above the repaired
Leterrier lamination. This is a small improvement, **not a breakthrough or
a new SOTA claim**. The 54D and 56D bounds remain unchanged.

The complete rational Gram is [gram55.json](gram55.json). Run:

```sh
python research/cross_coset55_20261116/verify.py
```

The auditor requires NumPy and python-flint. It reconstructs the published
54D and 56D parents, proves the new layer bounds by rational arithmetic,
and checks the negative certificates below. No floating-point search or
shortest-vector completeness assumption is used in these proofs.

## 55D: interpolate between two nontrivial cosets

Use Leterrier's rank-54 lattice L with minimum 10 and determinant 3^25.
The non-narrow point Q satisfies

```
3Q in L,    Q not in L,    ||Q||² = D = 28/3,
minimum(Q+L) = minimum(2Q+L) = D.
```

These are the same published inputs as in the
[preceding construction](../tilted55_20261116/PROOF.md), independently of
the faulty 55D lamination proposition. The verifier repeats their explicit
source-to-Gram binding.

The previous tilt used a minimum vector of L, of norm 10. Instead use Q
itself, of smaller norm 28/3. Although Q is not in L, adding it takes one
nontrivial coset to the other, whose minimum is equally large. Set

```
alpha = 1/21,
s = alpha Q,                 ||s||² = 4/189,
h² = 10/9 - ||s||² = 206/189,
L55 = < (L,0), ((22/21)Q,h) >.
```

For w in Q+L, both w and w+Q have squared norm at least D. Hence

```
||w+alpha Q||²
  = (1-alpha)||w||² + alpha||w+Q||²
      - alpha(1-alpha)||Q||²
  >= D - D alpha(1-alpha).
```

The layer checks are complete:

* Layer zero inherits minimum 10.
* In layers ±1, the displayed inequality plus h² gives exactly 10.
* In layers ±2, the horizontal part is in 2Q+2alpha Q+L. The reverse
  triangle inequality gives at least `D(1-2alpha)²`; adding `4h²` gives
  **12**.
* In layers ±3, the horizontal coset is 3s+L. Since `9||s||²=4/21<10/4`,
  zero is a nearest parent point: for every nonzero v in L,
  `||v+3s|| >= sqrt(10)-||3s|| >= ||3s||`. The resulting squared norm is
  `9(||s||²+h²)=10`.
* For every absolute layer index at least 4, the vertical contribution
  alone is at least `16h²=3296/189>10`.

The coefficient vector `(-3Q,3)` in the displayed basis attains norm 10.
The positive Schur complement is h², so

```
det(L55) = 3^25 * 206/189,
delta55² = (5/2)^55 / (3^25 * 206/189).
```

The exact squared density ratios are **1547/1545** over our previous
tilted lattice and **105/103** over the repaired untilted lamination.
The full values and certified gain intervals are in [results.json](results.json).

## External comparison

The valid source inputs are Leterrier (2022),
[Corollary 1.4 and the Mordell–Weil height bounds](https://doi.org/10.1007/s40993-022-00321-0),
and his [2024 thesis](https://infoscience.epfl.ch/record/307330/files/EPFL_TH10417.pdf),
Theorem 3.2.7 and the discussion preceding Proposition 3.2.22.

Proposition 3.2.22 prints center density **116,783.149646…**, which remains
larger than our new number. Its displayed construction has an explicit
norm-6 vector and therefore does not prove its stated minimum 10; the
[prior correction](../tilted55_20261116/PROOF.md#external-source-correction)
documents this precisely. We have not disproved the existence of another
packing attaining the printed density. Here the squared ratio to the
printed claim is **63/103**, so **we make no 55D world-record claim**.

## A tempting larger 55D deformation fails

A full Gram deformation, constrained by the supplied 54D roots and 2,187
closest points, left the horizontal Gram unchanged and suggested another
0.567% density gain. BKZ with block size 32 did not find a vector below 10.
Neither observation is a minimum certificate.

We generated additional closest points algebraically. Over F_(3^6), take
`b^9=b` and `c³+c=b⁴`. Then

```
x(t) = -b³ t^7+c,
y(t) = t^14+b t^7+b²
```

satisfy `y²=x³+x+t^28`. Translating the base variable gives more points.
Their lattice coordinates are checked against the reconstructed Gram.
One such point supplies a vector of squared norm
**19951291/2000000 = 9.9756455** in the proposed 55D Gram, rejecting it.
Even using that value as an upper bound on its true minimum gives a
density upper bound below our preceding 55D lattice; this comparison is
also checked exactly.

[rejected_deformation55.json](rejected_deformation55.json) supplies the
complete new border of that Gram and the exact short vector. The auditor
reconstructs and rejects the proposal directly. Thus the reported positive
bound above does not rely on the incomplete closest-point list.

## 54D: an index-nine overlattice target

Adding Q already gives an index-three overlattice. We tried adding a
second independent order-three dual class R, giving

```
M54 = L + Z Q + Z R,    [M54:L]=9.
```

If its minimum were 28/3, its density ratio to Leterrier54 would be
`9*(14/15)^27 = 1.397125689…`: a **39.71%** improvement. This is a
conditional target, not an achieved density.

The 24 initial candidates all failed. Each has an explicit vector of
squared norm at most 8. This is a finite screen, not an impossibility
result for other choices of R.

A further **300** adaptive trials, with BKZ block size 32, also failed.
Their exact matrices and witnesses are in [adaptive54.npz](adaptive54.npz).
Each has a supplied vector of squared norm at most **22/3**. The auditor
checks all 300 extension indices, witness memberships, norms and density
upper bounds. These are short-vector upper bounds, not assertions that
the actual minima equal the values returned by reduction. The initial
and adaptive screens used different reduction block sizes.

## 56D: an index-29 overlattice target

The [published MW56 parent and its exact reconstruction](../cm58_20261116/PROOF.md#published-parent-and-exact-gram-binding)
have minimum 10 and determinant 29^8. Adjoining
a nonzero dual class gives an index-29 overlattice. A proved minimum of
**258/29** would already yield density ratio

```
29*(129/145)^28 = 1.098008807… .
```

This would improve MW56 by **9.80%**. Again, it is only a target.
All 24 sampled candidates have exact short-vector counterexamples; the
largest of their supplied minimum upper bounds is **200/29**.

The 48 initial candidates and their witnesses are in
[overlattice_screen.json](overlattice_screen.json). For each one the auditor
checks the dual condition, the exact extension index, witness membership,
its norm, and the rational squared density upper bound
`index²*(witness_norm/10)^dimension < 1`.

The finite discovery screen can be repeated with NumPy, python-flint and
fpylll; its output is separate from the proof:

```sh
python research/cross_coset55_20261116/search_overlattices.py --output /tmp/overlattice_replay.json
python research/cross_coset55_20261116/search_adaptive54.py --output /tmp/adaptive54_replay.json
```
