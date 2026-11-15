# A six-translate sphere packing in dimension 59

This construction has packing fraction

\[
\Delta_{59}=2.63760259072000297875\ldots\times10^{-12},
\]

with center density **274,522.589741149678…**, or
`log2(delta)=18.066565…`.

An external historical reference is **Elkies's section E59, reported by
Alexander Vardy in *Density doubling, double-circulants, and new sphere
packings*, Trans. AMS 351 (1999), 271–283, §4**
([paper](https://doi.org/10.1090/S0002-9947-99-02169-8)). Its printed lower
bound is `log2(delta(E59)) >= 17.75…`. Our packing exceeds this reported
historical guarantee. **An exact percentage over the actual Elkies packing,
and a comparison with the strongest subsequent published 59D result, remain
unestablished.**

Separately, the repository's certified Gamma72 59D section has center density
`262144`. The new construction improves that **internal** reference by
**4.722057243785735…%**. Both the new packing and this exact internal ratio
pass two implementations. The internal reference is stronger than the floor
obtained from Vardy's printed prefix, so it remains useful as a progress check.

The idea is a balanced rank-five normal in the dual of Ne64, followed by a
six-point antipode construction. Five dual minimum-length vectors have squared
length 4. In their numerator Gram, consecutive vertices of a five-vertex path
have inner product -3; the other pairs have inner product -2. This balances the
six projected labels more effectively than the previously tested root normals.

## Exact parameters

| Quantity | Value |
|---|---|
| Ambient | Published Ne64, dimension 64, minimum 12, Gram determinant `3^32` |
| Primitive normal | `K = (14 I_5 - 2 J_5 - Adj(P_5))/3` |
| Normal determinant | `33800/81` |
| Dual labels | `0, e_1, e_2, e_3, e_4, e_5` |
| Maximum squared normal distance | `45/104` |
| Dimension / translates | 59 / 6 |
| Period Gram determinant | `33800 * 3^28 = 773235584977681800` |
| Exact packing squared minimum | `1203/104` |
| Shift denominator in the saved period basis | 1560 |

The exact center density is

\[
\delta_{59}
=\frac{6}{3^{16}\sqrt{33800/81}}
  \left(\frac{1203}{416}\right)^{59/2},\qquad
\delta_{59}^2=\frac{3^{33}401^{59}}{2^{296}5^2 13^{61}}.
\]

Multiplication by `pi^(59/2)/Gamma(61/2)` gives the packing fraction.
`result.json` records the exact rational squared density and comparison ratio;
the displayed decimals are evaluations of these exact expressions.

## Why this defines a packing

Let `L` be Ne64 with the archived Gram `G`. Rows of the integer matrix `A` in
`data/normal_coefficients.txt` describe normal vectors in the **dual basis** of
`L`. Their Gram is `K = A G^-1 A^T`. A unit minor and Smith invariant factors
all equal to 1 prove that these rows generate the whole primitive normal
`L* intersect U`, where `U` is their real span.

An integer primal row `x` has projected dual label `x A^T`. For each of the six
labels `s`, the saved integer lift has that label. Its projection to `U^perp` is

\[
y_s=x_s-sK^{-1}AG^{-1}.
\]

The saved 59-row integer matrix `N` generates the entire kernel of `A`, and
`Q=N G N^T` is the period Gram. Its free-coordinate block is an identity matrix;
together with the unit minor, this proves completeness of the kernel basis.
The saved rational rows `c_s` satisfy `c_s N=y_s`. Thus an entirely explicit
description is the union of the six translates `Z^59 + c_s`, with metric `Q`.

For any difference between two centers, the corresponding ambient difference
has squared norm at least 12. Its normal part has squared norm at most
`beta=45/104`. Consequently its perpendicular part has squared norm at least
`12-beta=1203/104`. Differences within one translate have squared norm at least
12. The translates are distinct: equality would give a nonzero ambient vector
lying wholly in the normal space, of squared norm at most `beta<12`.

An explicit norm-12 ambient vector in `data/contact.json` projects to a
difference of squared norm exactly `1203/104`. The checker verifies its label,
membership in the appropriate translate difference, and exact norm. The
minimum is therefore attained, conditional only on the published parent
minimum theorem.

Balls of radius `sqrt(1203/104)/2` about these centers give the stated density.

## Replay

From the repository root, with Python 3, NumPy, SymPy and python-flint installed:

```bash
python constructions/dim59/verify.py
python constructions/dim59/verify_independent.py
```

`verify.py` re-parses the archived Ne64 catalogue Gram, checks exact rational
LDL positivity, the normal embedding, primitive section basis, every translate,
all cluster distances, the attained minimum, and the density comparison.
It does not import the search or generator.

`verify_independent.py` uses FLINT integer Hermite forms to obtain a different
complete period basis. It checks a unimodular change to the saved period,
uses leading principal determinants for positivity, and verifies the same
rational packing data. It also reconstructs the **comparison** period from the
published Gamma72 Gram and the archived rank-13 normal, proving its determinant
`2^23` and center density `2^18`. It does not import the first implementation.

The saved explicit period data can be regenerated from the normal embedding:

```bash
python constructions/dim59/generate.py
python constructions/dim59/verify.py --write-result
```

Run these verifiers without Python's `-O` flag. No optimizer, floating minimum
enumeration, search cap, or completeness assumption about a sampled vector pool
is part of the proof.

## Sources and comparison scope

The inherited mathematical input is the published minimum 12 of
[Nebe's Ne64](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html),
from the Nebe–Sloane Catalogue of Lattices. The original catalogue HTML is
[archived in this repository](../../research/novel_routes_20261113/dual64/CQ64b.html).
Source and data hashes are recorded by both verifiers.

The independent comparison also inherits the published minimum 8 of
[Nebe's Gamma72](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/Ext72.pdf).
The prior explicit rank-13 section and its full period audit are in
[the propagated benchmark package](../../research/novel_routes_20261113/audits/PROPAGATED_BENCHMARKS.md).
This is the strongest 59D value in that finite audited collection, not an
exhaustive catalogue of all published packings.

The named historical same-dimensional comparison is **Elkies's section E59,
reported by Alexander Vardy in *Density doubling, double-circulants, and new
sphere packings*, Trans. AMS 351 (1999), 271–283, Section 4**
([DOI and paper](https://doi.org/10.1090/S0002-9947-99-02169-8)). The table there
prints `log2(delta(E59)) >= 17.75…`. Our value has `log2(delta)=18.066565…`.
This exceeds the displayed historical benchmark; the printed lower bound with
ellipses does not supply an exact density for the actual Elkies lattice or
establish the best subsequent published 59D packing. We do not calculate a
percentage improvement over that lattice from the truncated entry.

Current public tables, the older Elkies section values discussed by Vardy, and
recent antipode work were consulted. The commonly used public table lists only
selected dimensions above 48, so its omission of 59D cannot establish priority.
No global optimality or global SOTA claim follows from this package.

[Discovery code and other dimensions](../../research/lower_bound_return_20261115/README.md).
