# General sphere-packing upper bounds in dimensions 45 and 48

| Dimension | Strict density upper bound | Internal predecessor | Additional decrease |
|---:|---:|---:|---:|
| 45 | <1.2613300683e-6 | 1.2656297919e-6 | 0.3397299611% |
| 48 | <4.0691468600e-7 | 4.0843570355e-7 | 0.3724007321% |

These apply to **unrestricted packings**, not only lattices. The main repository table retains its attributable external comparison; the internal changes above are separate and make no exhaustive priority assertion.

## Method and attribution

Reoptimized nonnegative Gegenbauer coefficients in a cyclic three-point PSD kernel, with signed logarithmic radial profile h(v)=v^192 exp[-216(v-1)], retaining the previous Fourier witness. The exact identity C_new/B_new=C_old/B_old preserves the full pair inequality while reducing the origin cost. The applicable theorem is the cyclic strengthening in [Cohn, de Laat and Salmon, Section 2.2, equation (2.4) and Theorem 2.4](https://arxiv.org/abs/2206.15373), not the stronger pointwise premise of their Theorem 1.1.

## Reproduction

Install the repository's python-flint requirement, then run from any working directory:

```sh
python -B verifications/general_angular_20261117/run_all.py
```

The manifest binds the intact audited inputs and scripts. The runner checks exact angular reconstruction, cap, sign structure and infinite tail; the entire compact domain with interval Hessians; both radial/Fourier/rounding implementations; and exact substitution against the preserved predecessor. Assertions must be enabled. No search or network is needed. The repository aggregate runner includes this package.

[Audit details](AUDIT.md) specify domain coverage, scale extension, regularity, 768-bit compact replay and 1024-bit radial/structure checks. This is independent mathematical/code review and higher-precision replay of the reviewed compact implementation, **not a separately authored compact solver or different arithmetic backend**. Both paths use FLINT/Arb. Legacy angular subdivision runs were not completed; the complete Taylor and structure/tail gates replace them. Frozen discovery flags and density_interval metadata are historical, not verified bounds; strict_density_upper_bound and exact replay govern.
