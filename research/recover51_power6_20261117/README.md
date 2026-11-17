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
