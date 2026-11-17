# Independent compact harmonic-matrix transfer audit:46D and47D

PASS lattice-only upper bounds46D<7.790880151e-7 and47D<5.186745468e-7.

The exact angular matrix polynomial is unchanged from the independently audited48D compact kernel (source SHAaca9b48134ba9456995cd6a0cfc0aafcc419f8acdf6fe10f42af5477ca59cc86). A separate exact normalized Gegenbauer recurrence verifies every entry of its target-dimension expansion and positive semidefiniteness of every coefficient matrix. Positive amplitude scaling preserves the already verified rank-two-plane sign inequalities. No new dimensional restriction on those geometric inequalities is introduced.

The constant harmonic matrix changes under restriction. The target origin Schur complement uses the scaled target constant matrix, not the old48D matrix. This and the common positive normalization q are independently checked. The complete compact-support collinear correction is reconstructed by the scalar verifier with amplitude included exactly once; no higher multiples or omitted root vectors lie in its support.

Each dimension uses its own Fourier-positive degree71 witness, target-dimensional mod2 annular caps and scalar budget. Fresh1024-bit replays passed633 annuli with4903 nodes in46D and8361 in47D, and all71 spatial roots prove the negative tail beyond180. The exact cutoff418/13 is an annular edge. Separate exact Fourier-factor/inversion, cap, contiguous-coverage, source-binding and normalization checks passed. Even46D volume is rational; odd47D uses a separately implemented alternating Machin lower bound and a squared rational comparison. These checks establish the displayed safe bounds strictly.

Audit limitations: scalar intervals share the producer implementation and FLINT/Arb backend. The algebraic checker is separately authored. Universal geometry is inherited from the complete independently replayed48D main+cap package, not merely from numerical samples or a claimed result. No unrestricted-packing or exhaustive global-priority claim is made.
