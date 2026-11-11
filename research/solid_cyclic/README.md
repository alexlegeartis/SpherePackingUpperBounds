# Solid-harmonic cyclic packing-density certificates

This package extends the prior annular search structurally: solid spherical
harmonics turn the full three-sided cyclic condition into polynomials in squared
edge lengths. Exact domain certification then replaces the stronger condition
of individual negativity at every compact triangle vertex.

All four certificates passed primary and independently implemented audits:

- 17D: density <0.018335552.
- 18D: density <0.013284728.
- 19D: density <0.009557491.
- 20D: density <0.006977702.

17D adds a previously unregistered dimension; the others tighten registered
bounds. SUMMARY.json has pinned comparisons. No optimality, external peer-review
or worldwide-priority claim is made. The existing mathematical framework is
Cohn-de Laat-Salmon, https://arxiv.org/abs/2206.15373. PROOF17.md supplies the
complete periodic cyclic summation, and each dimension has an exact parameter
note. In 18D the stored upper squared radius is 2.6999999999999997, not 27/10.

## Genuine use of cyclic relaxation

At squared sides (187/80,187/80,17/8), the exact 17D witness has positive individual
unnormalized kernel value 1.50660346262744226547e-6, but its cyclic sum is negative,
approximately -2.43693368692587691746e-8. Both signs were independently checked over
the rationals. This witness is infeasible for the previous pointwise-negative
annular cone. No claim is made that any finite numerical program was optimized
exactly or globally over its parameters.

## Reproduce

Install verification dependencies in requirements.txt, then run:

python run_certificates.py

Never use Python -O. Only NumPy and python-flint are required for verification.
Search also uses requirements-search.txt. An example discovery/replay sequence:

python search.py --dimension 17 --output candidate17.json
python export.py candidate17.json exported
python verify.py exported/dim17.json

Exact rational witnesses define the functions. Optimizer status is never proof.
The first 17D candidate was rejected by a far-domain positive value; its witness
and logs are retained. Discovery was repaired by checking every real angular
stationary point on each radial grid fiber, then the final full-domain checks
passed. Noncompetitive/rejected 15D,16D,21D,22D,23D experiments are not bounds.

search_balanced.py tests a diagonal congruence of the same solid basis for
conditioning at larger supports. Its outputs remain unverified unless separately
certified. The .venv directory is an execution dependency cache and must not be
included in a source archive or publication.
