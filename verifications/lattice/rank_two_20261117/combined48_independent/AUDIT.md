# Independent combined Fourier-SOS 48D audit, recovered after rollback

The pre-rollback independent audit passed: lattice packing density <= 3.527782891e-7. Exact rational density is in certificate.json, approximately3.527782890320842e-7. This is only about0.00208128% below the previously audited3.527856315e-7. No unrestricted-packing or worldwide-priority claim is asserted.

Recovery verified exact byte hashes of all frozen inputs: verifier5f2b6cab94ad213e104bd51f721c4ba4438a705aa3ab7743ddf041364c66e70b; candidate53cbfead2c02fc1f30b9a74a6883df45fff448e29142730d364254257698e9d2; basise17512ea2482bbab1bdec84abb5b33684999d8b5c27cf52af9f9be57579f09c6; producer certificate616909d02a2032d9b1c926ce58cb69efc77ad4511812df96d14f59e677b174a7. The independent script is reconstructed from that frozen verifier with1024-bit Arb precision, explicit matrix/basis dimension assertions, and exact safe-rounding assertion. A fresh recovery replay is recorded in replay.log; only its final PASS denotes recovery execution completion.

The original independent replay passed all120 radial bins/143 cells,97 profile branches, global Fourier and physical-tail signs, origin PSD, profile dual constraints and rational rounding, with its entire certificate JSON equal to the producer768-bit output. Recovery must reproduce this equality through check_sos.py.

The new additive Fourier polynomial is m(s)=b4(s-1)^T FM b4(s-1)+s*b3(s-1)^T FN b3(s-1), s>=0. Rational symmetric5x5 and4x4 matrices have positive leading principal minors, hence m>=0. check_sos.py separately reconstructs every coefficient by direct symbolic matrix multiplication. The2x2 quadratic multiplier is positive definite after VM00 is reset to1-m(0). Thus oldQ times multiplier plus qzero*m is nonnegative, and its value at zero stays exactly qzero.

The independent verifier asserts5/4/2 matrix sizes and12 matching basis/coefficient entries to exclude silent zip truncation. Nine additive physical basis polynomials equal qzero times repeated applications of the audited inverse-Fourier operator to1. check_fourier_transform.py separately inverts the full degree83 physical polynomial in the normalized Laguerre basis and checks the intended product-plus-additive Fourier expression and unchanged zero normalization.

All remaining full-line rank-two, radial, profile-grid and density machinery is unchanged from the independently audited quadratic48D certificate; the source diff was reviewed before rollback. Q48 and profile prerequisites are inherited from completed independent audits, not newly asserted. The first recovery attempt failed solely because the restored dependency archive was extracted to a different directory; the original48 certificate and Q48 file were then copied into the expected parent directory.

## Fresh post-rollback recovery result

PASS. Fresh1024-bit replay completed all120 bins/143 cells and97 profiles, with exact rational rounding. Separate exact degree83 Laguerre transform and direct symbolic SOS/PD checks passed. The full recovered independent certificate again equals the frozen producer JSON. Recovery is complete, not merely a restored prose conclusion.
