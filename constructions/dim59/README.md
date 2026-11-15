# A six-translate sphere packing in dimension 59

This construction has packing fraction

\[
\Delta_{59}=2.63760259072000297875\ldots\times10^{-12},
\]

with center density **274,522.589741149678…**, or
`log2(delta)=18.066565…`.

The strongest comparison currently verified here uses **Nebe's published
Gamma72 structure and the classical κ13 section of the Leech lattice**.
They give a 59D packing of exact center density `2^22/(9*sqrt(3))`, or
**269,064.727051458…**. Our improvement over this construction is
**2.028457148397469…%**. This 59D number is **derived here from published
inputs**, not quoted as a result of Nebe's paper. See
[the derivation and source survey](#literature-derived-comparison).

Separately, the repository's certified Gamma72 59D section has center density
`262144`. The new construction improves that **internal** reference by
**4.722057243785735…%**. Both the new packing and this exact internal ratio
pass two implementations. The κ13 comparison above is stronger than this
internal reference and the historical prefix reported for Elkies's E59.

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
`result.json` records the exact rational squared density and the earlier
internal comparison ratio. `literature_section_comparison.json` records the
stronger κ13 comparison used in the README table. Displayed decimals evaluate
these exact expressions.

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
python constructions/dim59/verify_external_comparison.py
python constructions/dim59/verify_literature_section.py
```

`verify.py` re-parses the archived Ne64 catalogue Gram, checks exact rational
LDL positivity, the normal embedding, primitive section basis, every translate,
all cluster distances, the attained minimum, and the density comparison.
It does not import the search or generator.

`verify_independent.py` uses FLINT integer Hermite forms to obtain a different
complete period basis. It checks a unimodular change to the saved period,
uses leading principal determinants for positivity, and verifies the same
rational packing data. It also reconstructs the **earlier internal comparison** period from the
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

## Literature-derived comparison

**This is a mathematical consequence of published constructions, not an
externally reported 59D record.** It is stronger than the two comparison
values used earlier, so the README now uses it for the numerical gain.

[Nebe, *An even unimodular 72-dimensional lattice of minimum 8*](https://arxiv.org/pdf/1008.2862),
Theorem 2.3 and §4.2, supplies a primitive section
`K1 = sqrt(2)*Leech` of the unimodular Gamma72 lattice. Its 48D perpendicular
`K2` has minimum 8. The classical
[κ13 catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA13.html)
supplies a 13×24 basis, Gram determinant **972**, and cites Conway–Sloane,
*Sphere Packings, Lattices and Groups*, Chapter 6.

The verifier reads that basis `B` and the
[published Leech coordinate basis](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/Leech.html)
`L`. It verifies `B B^T/8 = Gram(kappa13)` and that `C=B L^-1` is integral
with all 13 Smith invariant factors equal to 1. Thus κ13 is a **primitive**
section of Leech, independently of any normal found by our discovery code.
The source HTML and integer embedding are retained in
[comparison_sources](comparison_sources/) and
[literature_section_comparison.json](literature_section_comparison.json).

Primitivity composes, so `sqrt(2)*kappa13` is a primitive normal in Gamma72.
Its determinant is `2^13*972`. In a unimodular lattice a primitive normal and
its full orthogonal section have equal Gram determinants. The 59D section
contains `K2`, so its minimum is **exactly 8**, using Nebe's published result.
Consequently

\[
\delta_{\mathrm{classical}}^2
=\frac{2^{59}}{2^{13}\cdot972}
=\frac{2^{44}}{3^5},\qquad
\delta_{\mathrm{classical}}=\frac{2^{22}}{9\sqrt3}.
\]

Its occupied fraction is `2.585163653786987618…e-12`, compared with our
`2.637602590720002979…e-12`. Comparing exact rational squared densities gives
**2.028457148397469…%**. Replay with
`python constructions/dim59/verify_literature_section.py`.
The Γ72 structure/minimum theorems are external mathematical inputs; the
embedding, determinant, primitivity and comparison arithmetic are checked
exactly here. No isometry search or capped enumeration is trusted.

### What other works actually report

The source checks on 15 November 2026 distinguish a reported 59D number
from constructions in other dimensions and their consequences:

| Work | Reported result relevant to this comparison |
|---|---|
| [Vardy (1999), §4, p. 278](https://doi.org/10.1090/S0002-9947-99-02169-8) | His own `J59` has exact center density **131072**. He reports the denser Elkies section with a lower-bound logarithm beginning **17.75…**, attributed to personal correspondence. |
| [Nebe (2012), Theorem 2.3 and §4.2](https://arxiv.org/pdf/1008.2862) | Reports Gamma72 of minimum 8 and its scaled Leech section. **269064.727051… in 59D is our derivation above**, not a number printed in her paper. |
| [Chen–Hu–Li–Wang–Wu (2025), Table 1](https://arxiv.org/html/2505.02394v1) | Reports density improvements in **19, 20, 21, 23, 44, 45, 47D**; no 59D row. |
| [Dorofeev–Sun–Wang (2026), Table 1](https://arxiv.org/html/2607.20359v4) | Reports density improvements in **38–43D**; no 59D row. |
| [Cohn's table](https://cohn.mit.edu/sphere-packing/) | Covers **1–48, 56, 64, 72D**, with no 59D entry. Direct retrieval returned HTTP 502; the indexed copy was inspected and reported a crawl age of four months. |
| [Nebe–Sloane density table](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html) | The live table has selected dimensions above 48, including **54, 56, 64, 72, 80, 128**, and no 59D row. |
| [Felpix-Studios repository](https://github.com/Felpix-Studios/sphere-packing) | Its inspected README reports a **4D upper bound**, not a 59D lower bound. |

The absence of a 59D entry in these sources does not prove no stronger result
exists. The present survey has **not identified the strongest previously
reported 59D packing**. The meaningful quantified comparison currently
available is the **2.028457% gain over the derived classical construction**;
the older roughly 24% reported-guarantee comparison is retained below as
historical context.

## Verified external comparison

The full [Vardy paper](https://doi.org/10.1090/S0002-9947-99-02169-8),
including §4 on p. 278 and reference 10 on p. 283, was checked directly.
An [archived copy of the AMS PDF](https://web.archive.org/web/20231125005643id_/https://www.ams.org/journals/tran/1999-351-01/S0002-9947-99-02169-8/S0002-9947-99-02169-8.pdf)
is accessible; its SHA-256 and the source observations are retained in
[the source registry](../../research/novel_routes_20261113/audits/benchmark_sources.json).

Let `B` denote the **guarantee reported for E59**, rather than the unknown
actual lattice density. Reading `17.75..` as the printed decimal prefix gives

\[
2^{71/4}\le B<2^{444/25}.
\]

For our exact squared center density `D = 3^33*401^59/(2^296*5^2*13^61)`,
the percentage `g=100*(sqrt(D)/B-1)` therefore satisfies

\[
\boxed{23.675981\%<g<24.536216\%.}
\]

This interval is proved by rational integer-power comparisons; its endpoints
are rounded outward. The factors `V_59` cancel, so it is also the percentage
for occupied fractions. The high-precision endpoint evaluations are
`23.675981094714…%` and `24.536215572032…%`; these decimals are for display.

The checker independently establishes

```text
D^25 > 2^888 * (1 + 23.675981/100)^50
D^4  < 2^142 * (1 + 24.536216/100)^8
```

using only Python's exact integers and `Fraction`. It also checks the source
metadata and the saved packing density formula. It does not import discovery
code. Replay with `python constructions/dim59/verify_external_comparison.py`;
[external_comparison.json](external_comparison.json) records the result.

**What remains unresolved:** Vardy states a lower bound on the section density.
The bound `B<2^17.76` is about the **reported guarantee**, not about the actual
E59 lattice. Thus the interval above cannot certify an improvement over that
lattice. The paper gives no exact E59 Gram, determinant, or minimum, and its
section-specific reference is personal correspondence. A formula matching
the displayed digits is insufficient to recover those missing parameters.

The same page supplies an exact external control, Vardy's `J59`, with center
density `2^(59-42)=131072`. Our exact gain over it is `109.444114487571…%`.
It is weaker than the reported Elkies guarantee, so it is not the comparator
in the README table.

Further source checks included the full sphere-packing chapter of the
[Schütt–Shioda author draft dated 16 March 2017](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf),
§13.4, pp. 356–360, and
[Chen–Hu–Li–Wang–Wu (2025), Table 1](https://arxiv.org/html/2505.02394v1#S1.T1).
Neither supplies the missing exact E59 comparison. The draft does confirm
the separately used 56D formula in Theorem 13.35, equation (13.24), p. 357;
it is not the 2019 final edition. These checks do not establish a current
59D literature record.

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

Current public tables, the older Elkies section values discussed by Vardy, and
recent antipode work were consulted. The commonly used public table lists only
selected dimensions above 48, so its omission of 59D cannot establish priority.
No global optimality or global SOTA claim follows from this package.

[Discovery code and other dimensions](../../research/lower_bound_return_20261115/README.md).
