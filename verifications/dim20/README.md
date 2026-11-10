# Sphere packing in dimension20: verified upper bound

This package proves the unrestricted sphere-packing density bound

Delta_20 < 0.006979502.

It applies to all packings, not only lattices. Against the reference bounds U=0.0071270536033763 and L≈0.0040241650630681, it reduces the upper bound by about2.07% and closes about4.7553% of the additive gap. It does not achieve the requested20% gap reduction. No worldwide-priority claim is made.

## Reproduce

Requirements: Python3.10 or newer and python-flint0.9.0.

    python -m pip install -r requirements.txt
    python verify.py
    python independent_verify.py

Both commands must finish successfully. The verifiers use assertions, so do not run Python with -O or otherwise disable assertions. Typical verification time is seconds to a few minutes. No internet connection is needed after installing the dependency; no optimization solver is needed.

## Contents

- certificate.json: exact rational polynomials and17 rational Gram matrices.
- PROOF.md: definitions, theorem, full-domain verification argument and normalization.
- verify.py, verify_global.py, verify_annular.py: primary reproducible verifier.
- independent_verify.py: independently implemented algebraic and interval audit.
- verification.log, independent_verification.log: successful reproduction logs.
- independent_verification.json: independently verified exact objective and rational density upper bound.
- requirements.txt: pinned verification dependency.
- SHA256SUMS: file-integrity hashes.

The certificate is based on Cohn–de Laat–Salmon Theorem1.1: https://arxiv.org/pdf/2206.15373 . The comparison table is https://cohn.mit.edu/sphere-packing/ . Numerical discovery is not used as proof: both verification paths check all global hypotheses from the exact rational witness.
