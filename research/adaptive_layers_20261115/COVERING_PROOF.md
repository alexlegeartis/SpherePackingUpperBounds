# Improve the finite covering estimate with an odd Taylor polynomial

The preceding covering proof derived

```
D_n = (n+1)*(n/(n+1))^(n/2) * E[f(S)],
f(s)=(1+s)^(-n/2), S=sum(T_i^2), T uniform on the n-coordinate simplex.
```

Coxeter–Few–Rogers's covering theorem gives `theta_n>=D_n`. This published
theorem is inherited. Jensen used only `E[f(S)]>=f(E[S])`. A globally
valid cubic Taylor polynomial gives a stronger finite estimate:

```
E[f(S)] >= f(mu) + f''(mu)*Var(S)/2 + f'''(mu)*E[(S-mu)^3]/6.
```

The fourth derivative of f is positive for s>=0. Taylor's remainder
therefore has the required sign on either side of mu; no concentration
or asymptotic approximation is used. The exact Dirichlet moments are

```
mu = 2/(n+1),
E[S^2] = [24n + 4n(n-1)]/[n(n+1)(n+2)(n+3)],
E[S^3] = [720n +144n(n-1)+8n(n-1)(n-2)]/[n(n+1)...(n+5)].
```

Put `a=n/2`, `v=Var(S)`, `c=E[(S-mu)^3]`. The new bound is

```
C'_n = (n+1)*(n/(n+3))^(n/2)
       *[1+a(a+1)*v/(2(1+mu)^2)-a(a+1)(a+2)*c/(6(1+mu)^3)].
```

For every even parent used here this is rational and strictly exceeds
the preceding Jensen estimate. The auditor rechecks the same 200 color
and parent combinations, all at dimension<=96, then verifies

```
(N-1)^2 * delta_parent^2 * V_n_upper^2 * d^n < (C'_n)^2.
```

The colored-layer separation cases and density formula are unchanged.
All new depths, period constraints and density comparisons are checked
with exact fractions and a Machin enclosure for pi. Horizontal offsets
remain existential. The seven incremental gains are **0.011–0.040%**;
they are smaller than the new deep-plane, explicit-section and adaptive
counting gains and are reported separately for that reason.

[Exact checker](verify_taylor_covering.py) · [Results](covering_results.json) ·
[Preceding covering proof](../covering_and_deep_cosets_20261115/COVERING_PROOF.md).
