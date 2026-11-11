# Improved six-dimensional unrestricted density certificate

This frozen computer-assisted certificate proves density strictly below
**0.409135399** for arbitrary sphere packings in R^6, using the standard
Fourier/PSD-kernel facts and the published neighbor-count bounds in PROOF.md.

The pinned paper-reported comparator is **0.4102905** in
[arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1).
The strict decimal decreases that upper by approximately **0.281532%** and
closes approximately **3.093223%** of its additive gap to the E6 lower density
pi^3/(48 sqrt(3)). The previous repository result was 0.410189987.
These are reported-source comparisons, not worldwide-priority or optimality
claims. The literature witnesses were not independently reconstructed here.

## Replay

    python -m pip install -r requirements.txt
    python verify.py
    python audit_independent.py .

No solver or network is required after installing dependencies. The scripts
work from another working directory when supplied their absolute paths; the
secondary audit must also receive this package's directory as its argument.
Never use Python optimization flags; the verifiers reject python -O. Tested
versions are NumPy 2.3.5 and python-flint 0.9.0.

The full producer replay passes exact PSD, global radial/Fourier signs, origin
Schur saturation, 213 pair boxes, 4,359 far boxes and 87 compact cyclic boxes.
The Fourier function needs no Gaussian repair. Independent review is recorded
in INDEPENDENT_AUDIT.md when finalized.

The full outward density enclosure is centered at

    0.4091353980508900600594412074029310240928937752092733306302830174062684616768

with radius below 7.5e-77. The stated strict decimal is rounded upward.

The proof data are exact rational entries in witness.json,
witness_psd_exact.json and prerequisites.json. discovery_snapshot.json is
provenance only. It is a rational 99/100 mixture of a new stronger numerical
checkpoint with 1/100 of the previous source. All conditions are recomputed
directly on the final data. The unmixed new checkpoint is not a bound.

This is a computer-assisted certificate, not external peer review or a formal
proof-assistant verification. Published spherical SDP count estimates remain
external mathematical inputs.

The repository wrapper also verifies the registry entry and frozen source hashes, asserts both Gaussian repair weights are zero, and runs the separate full independent audit in `independent/`. The primary standalone entry point is `verify_primary.py`. Earlier proofs remain in the parent package.
