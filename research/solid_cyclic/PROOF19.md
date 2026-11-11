# Solid-harmonic cyclic certificate in dimension 19

The primary complete exact/interval replay proves Delta_19 < 0.009557491.
The density enclosure is centered at 0.00955749071723961053311040089623003976576655784428912332212152422727976367995.
An independently implemented audit passed; this is not external peer review.

Use exactly the definitions, explicit periodic cyclic proof, support-boundary
cases and verification methods in PROOF17.md, replacing dimension by n=19,
s=r^2 by the exact rational 2.375, T=R^2 by the exact rational
2.85, and all polynomial/matrix data by certificate19.json.
There are five radial hats, angular degree 16 and physical/Fourier degree 31.
The solid Gegenbauer recurrence uses this n and s. The density normalization is

pi^(n/2)/Gamma(n/2+1)*(s/4)^(n/2) * [P(0)+A_0[0,0]]/H(0).

All PSD blocks, exact Fourier identity, global Fourier and tail signs, full
coupled radial inequality, 20 compact triangle cells and 10 far-edge cells pass.
No triangle-domain exclusions were needed. See full19.log. No external neighbor
count theorem is used. The result is unrestricted, not lattice-only, and does
not claim optimality, worldwide priority or external peer review.

Run without Python -O:

python verify.py certificate19.json

The exact functions are frozen and identical to witness19/dim19.json; only
metadata and the outward decimal threshold differ.
