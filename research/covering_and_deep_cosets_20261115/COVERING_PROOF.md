# Use sphere-covering bounds to obtain deeper torus offsets

The earlier greedy argument only used that balls of total volume less
than a torus cannot cover it. A stronger classical obstruction comes
from **Coxeter–Few–Rogers's covering simplex bound**. This is a use of
a published covering theorem, not a claim to originate that theorem.

## A finite, explicit covering lower bound

Let theta_n be the infimum density of coverings of R^n by congruent
balls. Coxeter, Few and Rogers,
[*Covering space with equal spheres*](https://doi.org/10.1112/S0025579300002059),
Mathematika 6 (1959), 147–157, prove `theta_n >= D_n`, the regular
simplex covering density. The statement applies to arbitrary sphere
coverings, including periodic coverings with several centers per cell.
A precise statement is also in Fejes Tóth's
[2022 author survey, §2.2](https://arxiv.org/abs/2202.11358): D_n is the
density of n+1 unit balls centered at vertices of a regular simplex
of circumradius 1, measured within that simplex.

Here is a self-contained finite estimate for D_n. The regular simplex
has edge length `sqrt(2(n+1)/n)` and volume

```
W_n = sqrt(n+1) * ((n+1)/n)^(n/2) / n!.
```

At each vertex, the edge cone has Gram proportional to `I+J`, where J
is the all-ones matrix. Its solid-angle fraction is

```
Omega_n = sqrt(n+1)/(2*pi)^(n/2)
          * integral_(x_i>=0) exp(-(sum x_i^2+(sum x_i)^2)/2) dx.
```

The opposite facet is at distance `(n+1)/n > 1`, so a unit ball's entire
sector in this cone lies inside the simplex. Consequently
`D_n=(n+1)*V_n*Omega_n/W_n`.

In the integral put `x_i=s*T_i`, with T uniform on the simplex
`T_i>=0, sum T_i=1`. Integrating s exactly and cancelling constants gives

```
D_n = (n+1) * (n/(n+1))^(n/2)
      * E[(1+sum T_i^2)^(-n/2)].
```

The function `(1+u)^(-n/2)` is convex, and
`E[sum T_i^2]=2/(n+1)`. Jensen therefore proves the finite bound

```
theta_n >= D_n >= C_n := (n+1)*(n/(n+3))^(n/2).
```

This uses no asymptotic remainder. For n=1 it gives 1; for n=2 it gives
6/5, below the exact planar covering density `2*pi/(3*sqrt(3))`.
Its asymptotic value `n*exp(-3/2)` agrees with the classical theorem.
For the even-dimensional parents used here, C_n is an exact rational.

## Improved greedy selection

Normalize a parent lattice L to minimum at least 4 and covolume
`v=1/delta`, using the same convention as the preceding range audit.
Suppose j torus balls of squared radius d covered `R^n/L`. Their lifts
would give a periodic sphere covering of R^n of density

```
j * V_n * d^(n/2) / v = j * delta * V_n * d^(n/2).
```

If that expression is strictly below C_n, the covering theorem rules
out coverage. Thus N offsets can be chosen greedily whenever

```
(N-1)*delta*V_n*d^(n/2) < C_n.
```

No assumption that the torus balls are disjoint or injectively embedded
is made. The lifted covering argument accounts for all overlaps.

The [colored-layer proof](../priority_resolution_20261115/TORUS_PROOF.md)
then applies without a change to its distance cases or density formula.
For vertical root Gram Q, period index N and same-color minimum ratio
lambda, use `t=max(4-d,4/lambda)`. The result is a periodic packing with
squared separation at least 4 and center-density square

```
delta_new^2 = delta^2 * (2/t)^r / det(Q).
```

## Exact audit and attribution

The standard-library [checker](verify_covering.py) rechecks all the
preceding scan's 200 color/parent combinations, replacing 1 by C_n in
the forbidden-covering inequality. It chooses rational depths in units
of `10^-6`, bounds pi above by Machin's formula, and verifies

```
(N-1)^2 * delta^2 * V_n_upper^2 * d^n < C_n^2
```

with exact fractions. It separately checks the same-color distance
condition and every density square. The parent lattices and small
colorings remain the named published/audited inputs from the preceding
range package. No parent or target is above 96 dimensions.

The improvements over the checkpoint are new calculations here, not
numerical records quoted from Coxeter–Few–Rogers. In particular, the
old 81D number was not the strongest consequence available from the
published MW80 parent and classical covering geometry. The explicit
96D comparison is still to BCS1982's reported construction, whose exact
center density is `3^48/2^24`.

All new offsets in this covering-based construction remain existential.
Published covering and parent-minimum theorems are mathematical inputs;
this checker verifies the derived arithmetic, not formal proofs of those
external theorems.
