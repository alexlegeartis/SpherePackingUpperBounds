# Sharper conservative section guarantees in 53D, 55D and 89–91D

The same regular-simplex geometry also sharpens the Hermite constants
used in the previous Blichfeldt section controls. These are **derived
classical guarantees**, not explicit normals or dimension-specific
records printed in the parent papers.

**55D source correction:** Leterrier's 2024 thesis, Proposition 3.2.22,
gives the stronger guarantee **116,783.149646…**. The old 55D control
below remains a valid calculation but is superseded in the current table.
[Source and exact comparison](REMAINING_DIMENSIONS.md#the-55d-correction-from-a-previously-omitted-source).

Rogers's packing simplex theorem is stated alongside the covering
theorem in [Fejes Tóth2022, §2.2](https://arxiv.org/abs/2202.11358),
with attribution to Rogers1958. Let T be uniform on the simplex
`T_i>=0, sum T_i=1`, and put `S=sum T_i^2`. The same cone integral
as in [COVERING_PROOF.md](COVERING_PROOF.md), now using balls tangent
along the simplex edges, gives the packing-density upper bound

```
Delta_n <= (n+1)*2^(-n/2) * E[(1+S)^(-n/2)].
```

The exact Dirichlet moments and support are

```
S >= 1/n,
mu = E[S] = 2/(n+1),
Var(S) = 4*(n-1)/((n+1)^2*(n+2)*(n+3)).
```

For `f(s)=(1+s)^(-n/2)`, the second derivative is positive and
decreasing. Taylor's theorem therefore gives
`E[f(S)] <= f(mu)+f''(1/n)*Var(S)/2`. Explicitly,

```
Delta_n <= U_n := (n+1)*2^(-n/2) * [
  ((n+1)/(n+3))^(n/2)
  + n*(n-1)/(2*(n+1)^2*(n+3)) * (n/(n+1))^(n/2+2)
].
```

This is finite in every dimension used here, without an asymptotic
error term. It implies `gamma_n^n <= 4^n*U_n^2/V_n^2` for the Hermite
constant. The checker encloses every odd-dimensional square root above
by an exact rational, bounds V_n below using Machin's pi enclosure,
and proves the chosen rational gamma bound by an integer-power check.
It also proves that every bound used is strictly stronger than the
**unrounded** Blichfeldt expression; this is not just a rounding gain.

For a lattice of center density delta_n, a primitive shortest dual
vector supplies a perpendicular section with

```
delta_(n-1)^2 >= (delta_n^2)^((n-1)/n) / gamma_n.
```

Recursive rational squared-density floors are certified by raising
both sides to the n-th power. The method is the same as the earlier
section auditor, with the sharper proved Hermite bound substituted.

| Dimension | New conservative center-density guarantee (approx.) | Published parent |
| --- | ---: | --- |
| 53 | 17,624.661755620… | Leterrier54, Corollary1.4 |
| 55 | 55,694.295039977… | Shioda MW56 |
| 89 | 71,446,504,435,357.756… | Shioda MW92 |
| 90 | 361,283,064,288,202.736… | Shioda MW92 |
| 91 | 1,861,659,817,306,416.211… | Shioda MW92 |

The parents are the same source-attributed lattice theorems used in
the preceding audit. No minimum above 96D, explicit dual embedding, or
worldwide-record improvement is inferred. The numbers improve our
previous conservative section controls and are included in the full
range table in that capacity.

Replay the standard-library [checker](verify_sections.py). Exact
gamma bounds and density squares are in [sections_results.json](sections_results.json).
