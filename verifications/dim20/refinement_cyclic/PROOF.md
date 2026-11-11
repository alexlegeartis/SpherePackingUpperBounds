This uses a cyclic generalization, not the literal pointwise-sign hypothesis of Theorem 1.1. The full finite-sum implication is in REVIEW.md.

# Solid-harmonic cyclic certificate in dimension 20

The primary complete exact/interval replay proves Delta_20 < 0.006977702.
The density enclosure is centered at 0.00697770162755537123555216263623502973315245883428172144901608211093138178049.
An independently implemented audit passed; this is not external peer review.

Use exactly the definitions, explicit periodic cyclic proof, support-boundary
cases and verification methods in [the 17D framework proof](../../dim17/PROOF.md), replacing dimension by n=20,
s=r^2 by the exact rational 2.5, T=R^2 by the exact rational
3.0, and all polynomial/matrix data by certificate.json.
There are five radial hats, angular degree 16 and physical/Fourier degree 31.
The solid Gegenbauer recurrence uses this n and s. The density normalization is

pi^(n/2)/Gamma(n/2+1)*(s/4)^(n/2) * [P(0)+A_0[0,0]]/H(0).

All PSD blocks, exact Fourier identity, global Fourier and tail signs, full
coupled radial inequality, 20 compact triangle cells and 10 far-edge cells pass.
No triangle-domain exclusions were needed. See full.log. No external neighbor
count theorem is used. The result is unrestricted, not lattice-only, and does
not claim optimality, worldwide priority or external peer review.

Run without Python -O:

python verify.py

The exact functions are frozen and identical to witness20/dim20.json; only
metadata and the outward decimal threshold differ.
