# Contact geometry of the existing 36D packing

**No new density record.** This folder now contains an exact local-optimality
certificate for the 36D lattice already committed in `constructions/dim36/`.
The certificate explains why continuing infinitesimal lattice deformations
at that endpoint cannot yield another improvement.

## Exact positive identity and perfection

Let Q be the committed integer reduced Gram matrix, with physical Gram Q/20
and minimum 160 in Q units. `contacts36.npz` contains94,050 antipodal integer
vectors z. Every norm `z Q z^T=160` is rechecked exactly. No completeness claim
for this contact list is needed.

For D equal to the old contact-preserving direction in reduced coordinates,
the values `z D z^T` and weights in the new eutaxy identity are:

| Slope | Number of antipodal contacts | Positive weight numerator |
|---:|---:|---:|
| -8 | 288 | 1326 |
| -4 | 5472 | 2581 |
| -2 | 59904 | 698 |
| 0 | 28386 | 1326 |

With common denominator 417590880, the exact identity is

```
(sum_z numerator(z) z^T z) Q = 417590880 I.
```

The 666-dimensional space of symmetric36-by 36 forms is spanned by these
contact outer products. `perfection36.json` selects666 contacts whose
constraint matrix has rank 666 modulo the prime 65521, which proves full rank
over the rationals.

For completeness, here is the local argument. Normalize nearby forms to keep
minimum 160, and write Q+E. Every displayed contact satisfies `z E z^T>=0`.
The positive identity gives `tr(Q^-1 E)>=0`, with equality only if every
contact value vanishes. Full rank then forces E=0. On the unit slice of this
closed constraint cone the trace has a strictly positive minimum. It follows
that `log det(Q+E)-log det(Q)=tr(Q^-1 E)+O(||E||^2)>0` for every sufficiently
small admissible nonzero E. Thus Q is a strict local density maximum among
lattices, modulo scale and isometry, conditional on its previously committed
minimum certificate. This is not a global optimality theorem.

Run:

```bash
python research/density_breakthroughs_20261113/contact_walk/verify36_local.py
```

## Other finite screens

`shell.cpp` is a floating Cholesky discovery enumerator with exact integer leaf
norms and an explicit node cap. Its complete-looking output does not constitute
an exact enumeration proof. The original 36D run visited126,699,355 nodes and
found the 94,050 contact representatives above without reaching its cap.

`simplex43.py` screened1,048,576 abstract rank-five Gram matrices with diagonal6
and off-diagonal entries from `{-3,-2,-1,0}`, using one dual-basis simplex and
explicit small-vector rejection checks. The best score reproduces the known
43D value 226.517274889. This finite search is neither a classification of
normal sublattices nor a global cluster optimum.

`periodic43.py` sets up all 946 metric and 215 relative-shift variables for the
published six-coset43D packing. Early sampled LP descents violated additional
known contacts. The successful numerical nonnegative-eutaxy check is saved
in `../sections/periodic43_fista_report.json` and its compact witness. It is
explicitly numerical, unlike the 36D identity and modular rank certificate.
