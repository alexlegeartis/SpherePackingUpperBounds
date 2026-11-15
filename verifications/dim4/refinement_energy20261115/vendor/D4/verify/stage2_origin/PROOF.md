# Stage 2: the written argument

This is the mathematical composition that `run.py` checks. It is taken from
the private review of 2026-10-20 (section "Mathematical composition"). The
changes are: links point to the paper instead of private files, and notation
follows [README.md](README.md). Section 4 of [`../../paper.pdf`](../../paper.pdf)
gives the same argument in shorter form.

## Normalization

Use center separation one, ball volume pi^2/32 and scalar integral pi^2/4.
The local-energy theorem (Theorem 2.1 of the paper) then gives density at most
(f(0)+A-beta)/8 when the completed local energy is at least beta on every
admissible neighbourhood, including the empty one. Its finite-patch and
periodization argument covers arbitrary upper density.

## Projection floor

Let E0=(z+u)^T M(z+u) be the fixed degree-zero Gram energy. The five cumulative
counts have radii

```
501473/500000, 1005643/1000000, 50379/50000,
50767/50000, 102597/100000.
```

For their indicator matrix J, the exact solution MW=J gives D=J^T W>0
and m=J^T u. Completing the square proves

```
E0 >= Q(n) = (n+m)^T D^(-1)(n+m).
```

The checker reconstructs M from positive integer rays and verifies 245 range
equations, 49 origin equations, 25 inverse equations and exact positive LDL
pivots. Bernstein feature sums give the actual closed cumulative counts.
The first radial piece includes its left endpoint; later pieces are open on
the left and closed on the right. The stage 1 sign proof covers all opposite
one-sided limits.

## Affine floors

For M=sum c_i v_i v_i^T, each supplied rational ray-center vector a gives

```
E0 >= C+p.z,
p = 2 sum c_i a_i v_i,
C = sum c_i (2 a_i v_i.u-a_i^2).
```

The gap is sum c_i(v_i.(z+u)-a_i)^2. Fresh exact Bernstein subdivision proves
p.phi(r)>=ell_j on each near shell and >=0 on both outer pieces. This gives
the affine count floor L(n)=C+sum ell_j(n_j-n_(j-1)). All twelve planes and
all seven closed radial pieces pass. The energy floor is max(Q,L_1,...,L_12).

## Count domain

The unchanged three integer-count energies are added once. Fresh spherical
count certificates and affine-latitude implications give the caps
(24,24,24,25,26), monotonicity, and

```
(23-n1)(n5-n1)>=0, (23-n2)(n4-n2)>=0.
```

These define all 165,503 enumerated integer vectors. Realizability of these
vectors is not assumed.

## Spherical surplus

The two unchanged spherical energies combine as Bbar*N-H_gbar. Their common
rational polynomial is reconstructed exactly, with ambient basis U_l/(l+1),
Gegenbauer parameter one. Auxiliary anchored kernels use transverse Legendre
harmonics, which suit S2. For exact cardinality N, positivity is required
only after substituting the constant feature sum N-1. Integer Gram rays prove
that restriction. For the bounded-count theorem, the unrestricted anchor
matrix is PSD and its nonnegative count slack requires N<=26.

The repair T=q/q0-1 has zero constant harmonic, nonnegative higher harmonics,
T<=-1 on the entire off-diagonal interval and T(1)=H. With delta=epsilon*H+eta,
adding epsilon*sum T plus delta*((N-1)-n) at fixed N preserves positivity.
For all N<=M use delta*((M-1)-n), which is nonnegative. The pair polynomial
decreases by eta and the cyclic polynomial by at least 3epsilon; the origin
charge increases by (N-1)delta or (M-1)delta. The common g is unchanged.

Three local proofs suffice for the final credit:

| Scope | Triangle nodes | Surplus |
|---|---:|---|
| N=24, all radii <=R1 | 222761 | approximately 0.0034998436290743415 |
| Every 0<=N<=26, radius R5 | 38969 | N times approximately 0.00011381806090652776 |
| N=25, all radii <=R4 | 80469 | approximately 0.003951599198732102 |

Pair, repair and full closed Gram-domain checks pass for each. A proof at a
smaller radius R_j applies only if n_j=N. Indeed R_j^2<2 makes
(r^2+s^2-1)/(2rs) increasing in each radius on [1,R_j], giving cosine at most
1-1/(2R_j^2). The checked angular cutoff dominates that exact value.

The driver starts with zero spherical surplus and takes maxima over only these
three proofs. Earlier fixed-N and weaker uniform proofs are not needed;
leaving them out still reproduces the exact final credit, so none of them is
assumed.

## Minimization and bound

The complete finite minimization is

```
beta = min_n [max(Q(n), L_1(n),...,L_12(n))
              + the three original count energies
              + max(applicable spherical surpluses,0)].
```

It gives beta approximately 0.0036901023146053216, attained at
(24,24,24,24,24). The N=25 and N=26 minimum floors are approximately
0.003960201005545091 and 0.0037781191961084527. N=0 is included.
Every other original local energy remains nonnegative and unchanged.

Grouping the count vectors by (N, j), with N=n5 and j the first index with
n_j=N, gives the 124 group floors that stage 3 refines. Within a group the
applicable spherical surplus is a constant.

The raw objective includes all graph charges, exactly 4769/200000000 in density
units, and is approximately 0.63000029144357012698. Subtract beta/8 from that
raw objective. The old origin credit is replaced once, not added again.
The reduced rational objective agrees exactly with the stored final result,
and the next 10^-12 grid point is the strict ceiling 0.629539028655. Only
K(0,0) changes, so every physical origin coupling, nonorigin value, pair
expression, cyclic expression and graph overlap charge is unchanged from V7.
