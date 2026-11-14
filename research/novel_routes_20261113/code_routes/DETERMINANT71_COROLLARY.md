# An even71-dimensional minimum-eight lattice has determinant at least eight

## Statement

Let L be a positive definite even integral lattice of rank71 with squared
minimum at least8. Then its Gram determinant satisfies

```
det(L) >= 8.
```

For squared minimum exactly8, its centre density is consequently at most
`2^34`. This is a bound for the stated integral, normalized class, not a
universal density upper bound for all71-dimensional lattices or packings.

## Proof

An odd-dimensional even Gram matrix is alternating modulo2 and therefore
has even determinant. If the positive determinant d were below8, then
d would be2,4, or6.

Use the discriminant norm convention `q(x+L)=(x,x) mod2Z`. Milgram's
signature formula says its normalized Gauss sum is `exp(7*pi*i/4)`, since
the signature of L is71. Classifying the tiny possible discriminant forms
is elementary:

* For a cyclic group of order d, a generator has q-value a/d with a coprime
  to d, taken modulo2d. For d=2,4,6, the only possibility with Brown
  invariant7 is respectively a=3,7,11. Thus the form is the opposite of
  the discriminant form of the one-dimensional lattice〈d〉.
* When d=4, the other group is `(Z/2Z)^2`. Its Gauss sum is a Gaussian
  integer divided by2, so if it has absolute value1 its Brown invariant
  is even. It cannot have invariant7. The checker also enumerates all
  nondegenerate forms on this group directly.

An anti-isometry between the discriminant forms of L and〈d〉therefore
exists. Glue their direct sum along its graph. The result is an even
unimodular positive definite lattice N of rank72. If w denotes the
generator of〈d〉, then w is primitive in N, its norm is d, and

```
N intersect w-perp = L.
```

The last equality follows directly from graph gluing: a vector with zero
〈d〉component must have zero discriminant class on both sides.

But N contains the nonzero vector w of norm at most6. The theorem in
`MODULAR_OBSTRUCTION.md` says that for any nonzero normal of norm at most10
in such a nonextremal even unimodular72-dimensional lattice, its
orthogonal complement contains a nonzero vector of norm at most6.
Applying it to w contradicts min(L)>=8. Hence d>=8.

## Inputs and exact check

The finite Gauss sums and the order-four classification are verified by
`verify_71_discriminants.py` using exact SymPy expressions. The modular
identity used in the final step is independently checked in
`verify_modular_obstruction.py`.

Standard discriminant-form facts are as in V. V. Nikulin, *Integral
symmetric bilinear forms and some of their applications*, Mathematics of
the USSR-Izvestiya14(1980),103–167, especially the discriminant-form and
overlattice machinery in Sections1.3–1.6:
<https://www.mathnet.ru/eng/im1677>.
Milgram's signature formula is given in J. Milnor and D. Husemoller,
*Symmetric Bilinear Forms*, Springer1973, Appendix4. The underlying
normalized Gauss-sum statement is also displayed with proof in
R. Miranda and D. R. Morrison, *Embeddings of Integral Quadratic Forms*
(2009 manuscript), ChapterIII, Section5, Theorem5.1:
<https://web.math.ucsb.edu/~drm/manuscripts/eiqf.pdf>.

No historical originality claim is made for this corollary.
