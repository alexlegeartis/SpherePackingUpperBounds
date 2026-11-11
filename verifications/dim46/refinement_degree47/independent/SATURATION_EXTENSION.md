# Contact duals also bound free scalar amplitudes

This lemma extends the existing degree-31 saturation proof without changing
any of its numerical inputs or weights. It does not extend the radial span or
remove physical negativity beyond squared radius 7.

Let Q have nonnegative normalized harmonic coefficients, constant coefficient
1 and Q(1)=B>0. Suppose the nonorigin scalar kernel has common radial feature
h, normalized by h(s)=1 at contact, and free amplitudes

```
K(x,x)=D h(x)^2,
K(0,x)=-A h(x),
K(0,0)=C.
```

The degree-zero Gram is

```
[ C     -A  ]
[ -A    D/B ].
```

PSD therefore gives C>=0, D>=0 and CD>=B A^2. Here A may be any real number.
The contact pair condition is f(s)<=2A-D.

If D>0, the identity

```
A^2-D(2A-D)=(A-D)^2>=0
```

gives 2A-D<=A^2/D<=C/B. If D=0, PSD forces A=0, so f(s)<=0<=C/B.
Thus in every case the necessary contact row is exactly

```
f(s)-C/B<=0.
```

This is the row used by the embedded radial duals: multiplying by exp(6*pi)
turns it into P(6*pi)-exp(6*pi) C/B<=0. Every selected pair row is at contact;
all other selected rows are the unchanged physical and Fourier sign rows.
The actual objective is f(0)+C, also exactly the dual objective. Consequently
the same inverse-defined dual weights bound the full free-amplitude 2-by-2
scalar cone, not merely the balanced choice C=BM, A=D=M.

If an angular proof gives B>=L_B, then C>=0 implies C/B<=C/L_B. Therefore every
such kernel satisfies the contact row at L_B and the same radial lower proof
composes with the angular bound. The zero-D case satisfies the row directly.

For the all-rate result, the angular certificate uses only necessary cap rows
and the radial certificate uses only contact, so neither depends on the
positive Gaussian rate. All finite angular harmonic degrees are covered by
the replayed integral-moment tail proof.

Restrictions that remain: squared separation 6, the fixed degree-111 source
plus normalized Laguerre degrees 0 through 31 with nonnegative source scaling,
positive Fourier normalization, nonnegative origin cost, and physical
negativity for squared radii at least 7. Nonconstant harmonic-dependent radial
matrices, multiple radial features, and larger correction spans are excluded.

This file is an additive theorem explanation in the new audit package. The
previous frozen saturation package was not modified.
