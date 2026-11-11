# Dimension 43 Gaussian cyclic certificate

Strict density upper bound: **0.000002629078004**. All primary exact/interval checks pass. A separate independent audit passed. This is compared with the pinned published upper 0.0000027169074727, not an exhaustive worldwide priority claim.

Within this repository run `python verify.py`. For a standalone copy run `python verify_primary.py` and `python independent/verify_all.py` with `python-flint==0.9.0`. No optimizer, network, external theorem-count input or neighboring research directory is required. The repository wrapper checks the registration and runs both proof paths. The primary entrypoint verifies the angular files, full angular and cyclic-domain proofs, complete Fourier and physical root isolation, radial coupling, exact Fourier identity, origin PSD and rational outward rounding. Do not run with Python optimization flags.

`PROOF43.md` gives the full cyclic mechanism and analytic tail/scale arguments. The certificate is not an instance of the literal pointwise triangle-sign hypothesis of Theorem 1.1; the summed cyclic form is used.

Discovery and rejected candidates are in the separate dimension-midrange research package. They are not proof inputs.
