# Ne64: an exact rejection of the global order-three multiplicity shortcut

This note concerns the 64-dimensional even 3-modular lattice Ne64, with
minimum 12 and Gram matrix `dual64/source.npz:G`. Write `H=3 G^{-1}` and
`D=L*/L`, so `D` has dimension 32 over F3.

## The proposed shortcut and its consequence

Suppose an automorphism tau of L satisfied

```
tau^2 + tau + I = 0,
G^{-1}(tau-I) integral.
```

The first identity makes tau fixed-point-free on the real vector space;
the second makes it fix every discriminant class. Every nonzero vector in
each class would therefore lie in a three-element orbit of vectors of the
same norm. Dividing the count of dual vectors of norms 14/3 and 20/3 by
three would prove that some class of norm residue 2/3 has minimum at least
26/3. Such a class would yield a 65-dimensional lattice with centre density
`3^17/sqrt(10)` by the three-layer construction.

**The required tau does not exist in the published automorphism subgroup.**
The calculation below does not exclude an additional automorphism outside
that subgroup, nor does it exclude a class of minimum 26/3 by another method.

## Exact matrices and convention

The primary catalogue page
<https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64b.html>
states that the first five generators supplied for
`L8,2.otimes.L_32,2` generate a subgroup of Aut(Ne64). Those matrices are
published at
<https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64a.html>.
The page is pinned as `CQ64a.html`, and its matrices are also available in
`cq64_generators.npz`, key `generators`.

Coordinates are rows, so every first-five generator g satisfies
`g G g^T = G`. The sixth generator does not preserve Ne64's G and is excluded.
For this convention, g fixes every class in D exactly when
`H(g-I) = 0 mod 3`. The reversed matrix product is not the right test.

The first two generators generate a group A of order 4896, the SL2(17)
factor. The next three generate a commuting group B of order 120, the
SL2(5) factor. Their intersection is exactly `{I,-I}`; together they generate
the published central product of order 293760.

`verify_ne64_discriminant_action.py` regenerates both groups directly from
the pinned HTML and tests all elements. The discriminant-action kernel of
each factor is trivial. In a central product of these factors, an element
of the full kernel commutes with each factor: its commutator with A lies in
the kernel intersected with A, and similarly for B. It must therefore lie
in the centre. The centre is `{I,-I}`, and `-I` does not fix D. Thus the
entire published subgroup acts faithfully on D.

## Partial fixed-subspace version also fails to give the proposed count

There are twenty elements tau of B satisfying `tau^2+tau+I=0`.
Each fixes a 16-dimensional subspace of D. No other nonidentity element of
B fixes a nonzero vector of D. In particular, the available fixed classes
number at most `20 * 3^16`, vastly fewer than the approximately `3^31`
classes required for the proposed global counting argument.

There is an additional exact structural obstruction for the desired
nonisotropic classes. Put `T=tau-I`. The quadratic identity gives
`T^2=-3 tau`, hence `T^2=0` on D. Since T has rank 16 on D, its image equals
its kernel. On the real space, the isometry identity gives
`||T x||^2=3 ||x||^2`. Consequently `T(D)` is totally isotropic for the
discriminant quadratic form. All tau-fixed classes are therefore isotropic.
The SL2(5) factor acts freely on each nonisotropic discriminant class orbit;
its order-three symmetries cannot supply multiplicity within the classes
whose target minimum is 26/3.

## Relation to the original construction

Nebe's primary paper, *Some cyclo-quaternionic lattices*, J. Algebra 199
(1998), 472–498, describes this construction on manuscript page 24:

<https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/cycloquat.pdf>.

It obtains the 3-modular form by a totally positive prime above 3 in the
real biquadratic field Q(sqrt(17),sqrt(5)), and proves modularity using a
Galois automorphism of the underlying unimodular lattice. This construction
does not itself assert an Eisenstein-unit automorphism acting trivially
on the discriminant group.

Local PDF and OCR working copies are ignored by version control. The exact
group claims above depend on the retained catalogue matrices and the
integer checker, not on OCR. The primary PDF can be downloaded again from
the displayed URL if needed.

## Replay and optional group caches

```
python research/novel_routes_20261113/code_routes/verify_ne64_discriminant_action.py
```

This regenerates both groups from `CQ64a.html` and repeats all assertions.
Add `--save-groups` to regenerate the ignored `sl2_17_group.npz` and
`sl2_5_group.npz` caches, each with key `group`. Those full-group caches are
not required by the checker. The small retained `cq64_generators.npz`
contains the six source generators; only its first five preserve Ne64.
