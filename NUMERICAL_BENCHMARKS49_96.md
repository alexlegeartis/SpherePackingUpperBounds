# Numerical spinless-bootstrap comparison benchmarks: dimensions 49–96

Comparison method: **Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini (2020), numerical spinless modular bootstrap**. [Primary paper, JHEP 12 (2020) 066](https://doi.org/10.1007/JHEP12%282020%29066) · [Original MIT dataset](https://hdl.handle.net/1721.1/125646).

These values are **numerical benchmarks, not certified sphere-packing upper bounds and not achieved packing densities**. They are converted from the authors' highest available truncation order N in `dataminusone.txt`. The source supplies optimized spectra, not sign-certified Cohn–Elkies functions; the general same-radius lift is conjectural (paper Section 2.3). They are useful as expected Cohn–Elkies comparison targets, with that distinction intact.

For Euclidean dimension n = c_total and first scaling gap Delta1, the displayed density is `(pi*Delta1/2)^(n/2)/Gamma(n/2+1)`. Digits are approximate, not outward-rounded rigorous inequalities. Changes between truncation orders are not error certificates. The numerical rows below remain separate from the certified results listed after the table. The lattice-only table remains restricted to dimensions at most 48.

| Dimension | Numerical density benchmark ≈ | Truncation N |
|---:|---:|---:|
| 49 | 3.157963967863e-7 | 463 |
| 50 | 2.2003929580188e-7 | 460 |
| 51 | 1.5321362440043e-7 | 457 |
| 52 | 1.0661233909836e-7 | 466 |
| 53 | 7.4137674631496e-8 | 469 |
| 54 | 5.1522853741015e-8 | 484 |
| 55 | 3.57847516992e-8 | 478 |
| 56 | 2.483937099234e-8 | 1025 |
| 57 | 1.7231960044411e-8 | 493 |
| 58 | 1.1947773131075e-8 | 496 |
| 59 | 8.2794955140058e-9 | 502 |
| 60 | 5.734447036655e-9 | 529 |
| 61 | 3.9696833663941e-9 | 511 |
| 62 | 2.7466447157762e-9 | 523 |
| 63 | 1.8994891623302e-9 | 517 |
| 64 | 1.3129980220119e-9 | 1025 |
| 65 | 9.071711718055e-10 | 526 |
| 66 | 6.2649431263453e-10 | 532 |
| 67 | 4.3246640051436e-10 | 541 |
| 68 | 2.9840040445823e-10 | 562 |
| 69 | 2.0580812987163e-10 | 547 |
| 70 | 1.4188806244807e-10 | 586 |
| 71 | 9.7780754670973e-11 | 553 |
| 72 | 6.7357963275328e-11 | 592 |
| 73 | 4.6382724003958e-11 | 568 |
| 74 | 3.1927056089977e-11 | 592 |
| 75 | 2.1968493577489e-11 | 592 |
| 76 | 1.5110672690196e-11 | 604 |
| 77 | 1.0389931119345e-11 | 586 |
| 78 | 7.1415093921072e-12 | 589 |
| 79 | 4.907031568306e-12 | 595 |
| 80 | 3.3705604159378e-12 | 616 |
| 81 | 2.3144224796372e-12 | 601 |
| 82 | 1.5887046500543e-12 | 628 |
| 83 | 1.0902004609935e-12 | 622 |
| 84 | 7.4788496688616e-13 | 622 |
| 85 | 5.1289798963799e-13 | 625 |
| 86 | 3.5163926988721e-13 | 658 |
| 87 | 2.4101067529564e-13 | 637 |
| 88 | 1.6513913882455e-13 | 646 |
| 89 | 1.1312036569892e-13 | 643 |
| 90 | 7.7465939867258e-14 | 661 |
| 91 | 5.3034941781127e-14 | 658 |
| 92 | 3.6299174396606e-14 | 655 |
| 93 | 2.4838006235112e-14 | 659 |
| 94 | 1.6991197192395e-14 | 667 |
| 95 | 1.1620381091102e-14 | 667 |
| 96 | 7.9452530317283e-15 | 678 |

[Machine-readable values, source checksums and reproducible conversion checks](research/cohn_bootstrap_baseline_20261117/). The full original dataset was checksum-checked and independently re-extracted; arithmetic verification of this conversion is not certification of the underlying numerical optimization.

## Certified results compared with these numerical benchmarks

- **49D:** certified density **<2.85041752e-7**, versus numerical benchmark ≈3.1579639678630e-7 (N=463), a9.7387573447% decrease from that numerical value.
- **50D:** certified density **<2.04993038e-7**, versus numerical benchmark ≈2.2003929580188e-7 (N=460), a6.8379867092% decrease from that numerical value.

[Exact values, method, complete verification and audit scope](verifications/general49_50_20261117/). These comparisons do not certify the external numerical optimizer or establish global priority. No new51D result or improvement is claimed here.
