This uses a cyclic generalization, not the literal pointwise-sign hypothesis of Theorem 1.1. The full finite-sum implication is in REVIEW.md.

# Solid-harmonic cyclic certificate in dimension 18

The primary full exact/interval replay proves Delta_18 < 0.013284728.
The density enclosure is centered at
0.01328472704022940068952086716325760124601533840031853772476992896674363378012,
with radius below 1.97e-78. An independently implemented audit passed; this is not external peer review.

Use precisely the definitions, proof, periodic cyclic summation, boundary cases
and validation methods in [the 17D framework proof](../../dim17/PROOF.md), with n=18, s=r^2=9/4, T=R^2=26999999999999997/10^16 (the exact stored rational, slightly below 2.7), and
the exact matrices and polynomials in certificate.json. There are five radial
hats, angular degree16 and Fourier/physical degree31. The solid-harmonic
recurrence uses n=18 and scale s=9/4; its normalized spherical kernel is the
18-dimensional Gegenbauer kernel. The objective is

pi^9/Gamma(10)*(s/4)^9 * [P(0)+A_0[0,0]]/H(0).

No neighbor-count theorem is used. The compact triangle proof covers20 radial
cell triples; the far proof covers10 radial cell pairs, with no triangle-domain
exclusions needed. The full coupled radial inequality, all PSD blocks, exact
Fourier identity and both global sign conditions pass. Details are in producer.log.

Against the previously registered outward0.013408931 this is approximately
0.92627% smaller. Against the pinned published0.0134853404450862 it is
approximately1.48763% smaller. These are improvements of upper bounds, not
claims of optimal packings or percentage changes in the unknown optimum.
No worldwide-priority claim is made.

Run without Python -O:

python verify.py

Exact functions are the same as witness18/dim18.json; only descriptive metadata
and the outward upper threshold differ.
