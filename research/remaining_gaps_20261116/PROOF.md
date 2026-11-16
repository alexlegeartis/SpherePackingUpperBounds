# Continuation on the nine dimensions without a new lower-bound record

**No new sphere-packing SOTA is claimed by this package.** The remaining
dimensions are 48, 54, 55, 56, 63, 64, 72, 80 and 92. Their named external
benchmarks remain those in the [range table](../range48_96_20261115/TABLE.md).
The record table has not been changed.

This continuation gives two exact obstructions and an exhaustive, bounded
search in a larger 64D arithmetic family. These results help avoid repeating
failed searches; they do not satisfy the objective of a new lower bound.

Run with Python, NumPy and python-flint:

```
python research/remaining_gaps_20261116/verify.py
```

The verifier reconstructs its Gram matrices from the archived primary
sources, checks every supplied integer vector, rebuilds the field algebra,
and repeats all finite search decisions exactly. [results.json](results.json)
contains the output. Numerical eigensolvers, lattice reduction, and convex
hulls were used in discovery but are absent from verification.

## 64D: the apparent deformation directions disappear

The parent is Nebe's even 3-modular lattice Ne64, minimum 12 and determinant
3^32, from the [Nebe–Sloane catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64b.html).
Its minimum and modularity are inherited published inputs. We independently
check the source Gram, determinant, positive definiteness and level-three
integrality conditions; we do not enumerate all 138,458,880 minimum vectors.

Several initially generated minimum-vector orbits span only 1,960 of the
2,080 symmetric quadratic-form directions. That is insufficient evidence
for an improving deformation. For minimum vectors u,v with inner product
6 or −6, the vector u−sign((u,v))v also has norm 12. Applying this identity
to the supplied orbits produces additional exact minimum vectors and fills
the missing directions.

[ne64_perfection.npz](ne64_perfection.npz) gives 2,080 integer vectors.
The verifier checks their norms are exactly 12 in the primary Gram and
checks that the matrix of products v_i v_j, i≤j, has rank 2,080 modulo
65,521. Thus its rank over Q is also 2,080: **Ne64 is perfect**. The
certificate makes no assertion that its vector list is a complete shell.

There is also a short, classical harmonic-theta proof of eutaxy. For any
harmonic quadratic polynomial P, the weighted theta series of Ne64 is a
holomorphic modular form of weight 64/2+2=34 on Gamma0(3). Its order at
infinity is at least 6 because the minimum is 12. At cusp zero, Poisson
summation gives the corresponding weighted theta series of the dual.
The dual minimum is 4 by 3-modularity, and the cusp width is 3, so its
order is also at least 6. The polynomial may change under the duality
isometry; it need not be a Fricke eigenvector. This does not change the
vanishing-order argument.

A nonzero form has total cusp order at most

```
34 [SL2(Z):Gamma0(3)] / 12 = 34/3 < 6+6.
```

Therefore every such weighted theta series vanishes. In particular the
minimum shell is a spherical 2-design. Its equal positive weights give
eutaxy. Perfection and eutaxy imply local optimality for lattice Gram
deformations by Voronoi's criterion. **This is a local conclusion, not a
global sphere-packing bound.**

## 64D: exact rejection of a distant integral-form search

Local optimality does not exclude distant lattices. We therefore searched
a larger arithmetic family using the underlying even unimodular lattice
[CQ64a](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64a.html).
Write its Gram as G0 and the Ne64 Gram as G. Consider

```
span_Q {G0, G, G G0^-1 G, G G0^-1 G G0^-1 G}.
```

These four forms are independent. To obtain *all* even integral forms in
this particular rational span, we vectorize a symmetric Gram with diagonal
entries divided by two, then saturate the resulting coefficient lattice
using HNF. This gives the four matrices B0,...,B3 in
[global64_forms.npz](global64_forms.npz). The verifier reconstructs this
basis, rather than trusting the saved matrices. It also checks invariance
under the five published Ne64 subgroup generators. We do not claim to
classify the entire space of invariant forms.

The search checks every

```
F = a B0 + b B1 + c B2 + d B3,   −32 ≤ a,b,c,d ≤ 32.
```

There are 17,850,625 tuples, giving 2,427,638 positive definite forms,
of which 2,237,409 have primitive coefficient vectors. Multiples have the
same packing density, so testing primitive vectors covers the entire box.
**Every positive definite form in this box has packing density at most
Ne64's.** This is an upper bound within this finite family, not merely a
failure to certify a better minimum.

Here is the exact finite proof. The commuting endomorphisms G0^-1 Bi span
a degree-four field. The verifier reconstructs their multiplication table,
checks an irreducible quartic minimal polynomial, and verifies that the
64D representation is sixteen copies of the regular representation.
For the field element beta corresponding to F,

```
det(F) = Norm(beta)^16.
```

The four elementary symmetric functions of its four real conjugates are
the integer polynomials in [norm_polynomials.json](norm_polynomials.json).
The verifier rebuilds these polynomials from principal minors of the
regular representation. All conjugates are real because the corresponding
operators are self-adjoint for the positive metric G0. They are all
positive exactly when all four elementary symmetric functions are positive:
in the reverse direction, the characteristic polynomial at a negative
argument has strictly positive coefficients. This gives an exact
positive-definiteness test for every tuple.

[scan64_witnesses.npz](scan64_witnesses.npz) supplies 630 nonzero integer
vectors. For each candidate, at least one of these vectors has F-norm mu
satisfying

```
mu^4 ≤ 2304 Norm(beta).
```

Since the actual lattice minimum is at most mu, its squared center density
is at most

```
(mu/4)^64 / Norm(beta)^16 ≤ 3^32,
```

which is Ne64's squared center density. [scan.py](scan.py) checks this for
every primitive positive tuple, using integer square roots twice to compare
fourth powers without overflow. Integer-arithmetic bounds are checked
before the NumPy operations. No shortest-vector completeness claim or
floating-point comparison enters this proof.

In discovery, two forms of determinant 2^32 appeared to have minimum 10
after weaker reduction. Stronger reduction found norm-8 vectors in both.
They are rejected by the finite certificate. A verified minimum 10 at that
determinant would beat Ne64; these two forms do not have it.

This box is specified in the supplied integral basis. The conclusion does
not exclude coefficients outside it, other arithmetic families, or
nonlattice packings.

## 92D: the binary code target cannot use a shift by three

The [previous exact Craig/code audit](../simplex_normals_20261116/PROOF.md#a-separate-attack-on-the-nine-missing-dimensions)
left an even binary [93,43,≥24] code as a conditional route to a
24.9189116% improvement over the MW92 guarantee. **No such code has been
constructed here.** Grassl's [binary n=93,k=43 entry](https://www.codetables.de/BKLC/BKLC.php?q=2&n=93&k=43)
reports distance between 18 and 24; its upper endpoint is not a construction.
Codes with these parameters over GF(3) or GF(4) do not supply the binary
input.

We can exclude ordinary cyclic constructions, and more generally every
even binary [93,43] code invariant under the shift by three coordinates
(index-three quasicyclic codes). This permutation has three cycles of
length 31. Over F2,

```
x^31 − 1 = one degree-1 factor times six degree-5 factors.
```

The verifier factors this polynomial exactly and checks ord_31(2)=5.
Since 31 is odd, the representation is semisimple. The fixed subspace
has dimension three; even parity restricts it to dimension at most two.
Every other simple summand has dimension five. Consequently any invariant
even code has

```
dimension = k0 + 5t,   0 ≤ k0 ≤ 2.
```

But 43 is 3 modulo 5, so this is impossible, independently of the minimum
distance. Cyclic invariance implies shift-by-three invariance and is also
excluded. A permutation consisting of 31 cycles of length three is a
different case and is not excluded. General binary codes remain open.

## What remains

There is no improved lower bound in any of the nine target dimensions in
this continuation. The [prior frontier proof](../cm58_20261116/FRONTIER.md)
and the [Craig guarantee-formula audit](../simplex_normals_20261116/craig_results.json)
remain applicable. Neither their obstructions nor this finite 64D search
establish optimality in any of these dimensions. A future record entry
still needs a packing, a certified minimum-distance guarantee and a strict
comparison with the named external benchmark.
