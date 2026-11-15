# BCS96 is a strict local maximum among lattice deformations

This explains an unsuccessful research direction and provides an exact
finite certificate. It does not assert global optimality. The denser
[periodic packing](TORUS_PROOF.md) is compatible with this local result.

## Reconstruction and published inputs

Use the physical coordinate description of P48q in §5 of BCS1982.
Integer coordinates in this document are divided by sqrt(12) to obtain
physical lengths. Generators are the 47 cyclic shifts of twice the
indicator of `{0} union QR47`, the vector `(-5,1,...,1)`, and the
vectors `6(e_i +/- e_j)`. The exact HNF has determinant `12^24`, giving
an even unimodular parent Λ48. Its **published minimum 6 is inherited**,
not re-enumerated in this audit.

Let N be the symmetric 48-by-48 integer matrix indexed by infinity and
the residues modulo 47, with

```
N_inf,inf = 5,
N_inf,i = N_i,inf = -1,
N_i,j = -5                 if i+j = 0 mod 47,
        Legendre(i+j,47)   otherwise.
```

Then `N^2=72I`. The map `D=N/6` satisfies `D^2=2I` and maps Λ48 into
itself. The auditor checks this in the complete HNF basis; it also checks
that every row N_i lies in Λ48. Define

```
Gamma96 = {(x,y) in (D^-1 Λ48)^2 : x-y in Λ48}.
```

BCS's theorem gives minimum 6. The block basis `(b,0)` and
`(D^-1 b,D^-1 b)` has covolume `2^-24`, hence Gram determinant `2^-48`
and center density `3^48/2^24`. The exact Gram determinant is checked.

## A finite set of minimum contacts

In the integer-coordinate convention above, take the 2,256 root rays
`v=6(e_i +/- e_j)`, `i<j`. Form contacts

1. `(v,0)` and `(0,v)`;
2. `(u,u)` and `(u,-u)`, with `u=vN/12=D^-1 v`;
3. For each `i<j`, put `u+=(N_i+N_j)/2`, `u-=(N_i-N_j)/2` and use
   `(u+,u-)`, `(u+,-u-)`, `(u-,u+)`, `(u-,-u+)`.

These are all in Gamma96. For the second family `2u=Dv` belongs to
Λ48. For the third family, `u+ +/- u-` is a row of N in Λ48, and
`D u+ = 6(e_i+e_j)`, `D u- = 6(e_i-e_j)` belong to Λ48.
Every contact has integer-coordinate norm 72, hence physical norm 6.
The auditor checks that there are **13,536 distinct unsigned rays**.

## Perfection

Let a symmetric form H vanish on every contact, and partition it as
`[[X,A],[A^T,Y]]`. The first family implies X=Y=0: the D48 root outer
products span all symmetric 48-by-48 matrices. Indeed both signs at
each pair force the off-diagonal entries to vanish and give
`X_ii+X_jj=0` for every pair, forcing the diagonal to vanish too.

The second family now gives `u^T A u=0` for all transformed roots u.
Since N is invertible, their outer products also span the symmetric
matrices, so A is skew-symmetric. The third family gives

```
u+^T A u- = -N_i A N_j^T/2 = 0  for every i<j.
```

The rows of N are a basis. A skew form vanishing on all their pairs is
zero. Thus H=0, proving perfection: the contact outer products span
`Sym96`, of dimension **4656**. This is a structural spanning proof,
rather than a numerical rank calculation.

## Strict eutaxy and local consequence

For the integer contact matrix S, exact integer multiplication gives
`S^T S=10152 I96`. After division by sqrt(12), the physical contacts w
therefore satisfy

```
sum_w w w^T = 846 I96.
```

Every weight `1/846` is positive. Perfection and strict eutaxy imply
that Gamma96 is extreme by Voronoi's criterion: a strict local maximum
of lattice packing density, up to rotations and scaling. Consequently
a sufficiently small lattice Gram deformation cannot beat this
construction. The conclusion concerns this lattice, not all lattices
or all periodic packings.

Replay [audit_bcs96.py](audit_bcs96.py); the compact result is
[bcs96_audit.json](bcs96_audit.json).
