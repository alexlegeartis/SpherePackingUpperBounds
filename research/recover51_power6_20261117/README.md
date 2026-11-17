# Fresh 51D compact-kernel reconstruction

Research only. No new sphere-packing upper bound is established here.

The unpublished power-six checkpoints did not survive workspace recovery. This experiment instead starts from the published 51D power-eight certificate at commit `34f59648d034e92a52768c95af96221a8f589ac9`. It does not claim to recover the lost numerical candidate.

`screen.py` tests changing the common radial power while holding all harmonic matrices fixed. Powers 4–7 fail sampled cyclic geometry; at power 6 an optimized violation is approximately 645 (before normalization). Simply editing the exponents is invalid. Powers 8 and above pass this numerical screen only.

`refit.py` constructs a fresh semidefinite program for thirty 3-by-3 harmonic blocks and their origin extension, retaining the published Fourier-positive scalar polynomial. It minimizes the origin cost under sampled cyclic, support-crossing Bernstein-control, and radial residual constraints. Existing annular allowances serve only as sampled targets. Every output is marked numerical-only; exact PSD repair, full geometry, annular/tail proofs, normalization, and independent replay remain necessary.

```sh
python -m pip install numpy scipy sympy python-flint cvxpy
python bootstrap.py
OPENBLAS_NUM_THREADS=1 python screen.py
OPENBLAS_NUM_THREADS=1 python refit.py 6
```

Do not put candidate values in the certified tables.

## Guarded reconstruction

`refit_guarded.py` uses a shifted-Legendre radial basis and positive PSD margins. The saved power-six iteration 1 passed its sampled/optimized triangle screen and had positive matrix eigenvalues, but an exact support-crossing angular audit found a positive interior extremum (about 0.09738 before normalization). It is rejected and is not a bound.

`cap_extrema.py` isolates all derivative roots of the six exact rational angular controls. `refit_capguarded.py` incorporates these critical points as additional SDP cuts and records resumable cut sets. `prepare_producer.py` rounds a numerical candidate to exact rationals, checks all block principal minors, sets the origin Schur complement exactly, and prepares new producer checks without copying historical PASS receipts. Candidate geometry and radial inequalities still require fresh proofs and independent audits.

A failed solver now stops rather than reusing stale variable values. Earlier saved candidates retain their explicit numerical-only status. No lost historical power-six result is claimed recovered.

## Exact edge rejection and broader reconstruction

The next cap-guarded candidate passed all six exact angular controls but failed full cyclic geometry. An independent direct nine-entry contraction certifies a positive value about7.19105388238 at squared sides(1,1,67/52). See `independent_p6/edge_counterexample.py` and its exact frozen input; run that script to reproduce the rejection. Its conditional scalar value1.3598002e-7 is not a bound.

`refit_edgeguarded.py` adds all12 edge extrema to the cutting-plane loop; `separate_candidate.py` also searches all6 faces and multiple full-cube starts. These remain numerical screens, not proof substitutes.

`refit_joint.py` additionally optimizes a72-term normalized Laguerre perturbation of the Gaussian scalar function, together with the PSD compact kernel. Fourier normalization, sampled Fourier positivity, physical tail signs and annular residual targets are imposed. Any candidate still requires exact positive-factor repair and independent full-domain certification. This joint route addresses the fixed-scalar limitation; no result is claimed here.
