# 22D broader-profile refinement, 18 November 2026

Strict unrestricted sphere-packing density upper bound: **<0.0035915037512**.
Previous repository certificate: <0.0036110963151, [retained proof](../refinement_profile20261114/).
This lowers the upper bound by 0.542566%; no exhaustive worldwide-priority claim is made.

The profile is h(v)=exp[−20(v−1)]. Exact coefficients, angular PSD/cyclic inequalities, global Fourier and physical signs, complete radial coupling and strict density rounding are checked by both implementations. No external neighbor counts or lattice assumptions are needed. Independent implementations share FLINT/Arb root isolation.

Run with Python and python-flint installed, with assertions enabled:

```
python verifications/dim22/refinement_broad20261118/verify.py
```

See [the shared proof](../../profile_enrichment_20261114/PROOF.md) and [release audit](../../../research/broad_profiles_20261118/README.md). The discovery angular basis was padded to degree 25; canonical package metadata trims exactly zero trailing Gegenbauer coefficients to actual degree 13, without changing the polynomial. The retained angular_file field in certificate.json is discovery provenance, not a runtime dependency; the verifier reads local angular.json. Discovery flags are never used as evidence of validity.
