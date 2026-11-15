# Explicit MW92 antipodes in dimensions 86 through 91

The new objects are finite unions of explicitly specified translates of an
integer-kernel lattice. They are sphere packings with the stated minimum
**guarantees**; equality for their actual minima is not assumed.

## Reconstructing the published parent

Schütt–Shioda's *Mordell-Weil Lattices*, author draft of 16 March 2017,
Theorem 13.35 (p.357), applied at p=47, gives the Mordell-Weil lattice of

```
E: y² = x³ + 1 + t^48 over F_(47²)(t),
rank = 92, minimum = 16, Gram determinant = 47^14.
```

These parent facts are inherited published inputs. We do not claim a new
construction or minimum theorem for the 92D lattice. Section 13.4.2,
p.359, supplies explicit polynomial points using unitary matrices.

Use F_(47²)=F_47[i], i²=-1. The integer a+47b encodes a+bi. Every matrix
U=[[a,b],[c,d]] in [unitary_sections.json](unitary_sections.json) satisfies
U U^*=I. Its point is

```
x_U(t)=-(at+b)^16, y_U(t)=(ct+d)^24.
```

The exact finite-field checker verifies each unitary identity and each
curve equation. Every such section is disjoint from the zero section:
its coordinates have degrees at most 16 and 24, including in the smooth
chart at infinity. The singular fibres are irreducible of type II;
unitarity prevents x_U=y_U=0 at their cusps. Consequently the height
formula, Theorem 6.23, Eq.(6.15), p.123, is

```
<P_U,P_V> = 8 - (P_U.P_V)  (distinct non-opposite sections),
<P_U,P_U> = 16, <P_U,-P_U> = -16.
```

At a smooth point the intersection multiplicity of two graphs is the
minimum of the orders of their x and y differences. Thus the finite
intersection is deg gcd(x_U-x_V,y_U-y_V). At infinity it is the minimum
of 16-deg(x_U-x_V) and 24-deg(y_U-y_V), omitting identically zero
differences. The verifier computes all **16,110 pairwise heights** by
this direct polynomial-GCD method. Discovery used a separate method:
apply the relative unitary matrix, then count intersections among the
16 linear factors of the x difference. Both computations agreed.

The certificate contains an integral symmetric positive-definite 92D
Gram G and integral coordinates Y of all 180 polynomial points. Exact
checks establish

```
Y G Y^T = the directly computed height Gram,
row_module_Z(Y) = Z^92,
det(G) = 47^14.
```

The row-module check is HNF(Y)=[I_92;0]. Hence G describes the integer
span of these actual points. Its determinant agrees with the published
full Mordell-Weil lattice, so the index is one. This identifies the full
parent, rather than a finite-index sublattice. The published minimum 16
therefore applies to this precise explicit Gram.

## Complete periods and projected shifts

For each target n, let k=92-n and let B be the k by 92 integral normal
**label matrix** in [sections.json](sections.json). The physical dual
normal coefficient rows are Z=B G^-1; put K=B G^-1 B^T.

An exact unimodular HNF transformation T satisfies

```
T B^T = [I_k;0].
```

Its first k rows give integer lifts L with L B^T=I. Its remaining n rows
N form the complete integer kernel. Thus the period Gram is N G N^T,
with determinant det(G)det(K). This is a complete section, not an
arbitrarily chosen sublattice of that section.

For every distinct integer row s in the cluster S, the projected shift
in parent coefficient coordinates is

```
h_s = s L - s K^-1 B G^-1.
```

It lies in the kernel real space. The packing is the union of
Z^n N + h_s, with row-vector conventions. For two labels s,t, the
squared length removed by projection is

```
(s-t) K^-1 (s-t)^T <= beta.
```

Every preprojection difference is a parent lattice vector. For distinct
labels it is nonzero. Nonzero differences within one translate inherit
minimum 16. Consequently every distance squared in the union is at
least 16-beta. Since beta<16, the translates cannot coincide modulo the
period. With M=|S|, the resulting center-density square is

```
delta_n² = M² * ((16-beta)/4)^n / det(N G N^T).
```

All labels, shifts, diameters, determinants and density comparisons are
checked exactly. The discovery clique search may be incomplete; its
optimality is unnecessary. The saved finite cluster is the certificate.

`verify.py --export NEW_DIRECTORY` emits the complete period bases,
period Grams and rational shifts. No search, numerical optimizer or
floating-point enumeration is used by this proof auditor.

## External comparisons, including the full Rogers consequence

For 86D–88D the comparisons are the BCS1982 construction and its
primitive sqrt(2)A2/A1 sections. Their exact center-density squares are
2^84/3, 2^85 and 2^88. See the [existing exact reconstruction](../priority_resolution_20261115/CLASSICAL_SECTIONS.md).

For 89D–91D we also beat the section guarantee obtained from the **full
Rogers regular-simplex theorem**, beyond our earlier second-moment
evaluation. Let R_d be the theorem's numerical density upper bound.
For simplex-uniform T, S=sum(T_i²),

```
R_d = (d+1) 2^(-d/2) E[(1+S)^(-d/2)].
```

Jensen implies

```
R_d² >= (d+1)² * ((d+1)/(2(d+3)))^d.
```

Therefore the numerical Hermite-bound formula g_d=(4^d R_d²/V_d²)^(1/d)
is bounded **from below** by an exactly certified rational c_d. This is
a bound on the formula g_d, not a lower bound on the actual Hermite
constant. Iterating

```
D_(d-1)² = (D_d²)^((d-1)/d) / g_d
```

with c_d and outward integer ceilings produces an **upper endpoint** for
the full Rogers-derived section guarantee. The earlier second-moment
calculation supplies a lower endpoint. All cube/power/root comparisons
and the Machin enclosure for pi are rational. The reported external gain
uses the upper endpoint; it cannot benefit from weakening the Rogers
comparison. That upper endpoint is never inserted as an achieved packing
in the range table.

The Rogers theorem is stated and attributed to Rogers (1958) in
Fejes Tóth (2022), Section 2.2. [Source URLs, locators and PDF hash](sources.json).
