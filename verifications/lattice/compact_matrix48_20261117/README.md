# Certified compact harmonic-matrix lattice bound in48D

Lattice packing density is **<3.499004606e-7**. This does not bound unrestricted packings. Method: compactly supported two-feature positive-semidefinite harmonic matrices, complete rank-two lattice-plane inequalities, exact collinear correction, and mod-2 annular profiles with a positive-factor Fourier witness.

[Independent theorem review](THEOREM_REVIEW.md) · [Audit receipts](AUDIT_RESULT.json). The fresh audit verified all61,305 main and144,732 complementary cap leaves at384bits, and633 radial annuli at1024bits. Shared Jet/Arb implementation is disclosed; separate trace coverage, exact algebra and domain review are not described as another interval backend. Earlier rejected matrix candidates are not used.

Run `python run_all.py` for complete geometry/radial/tail replay. `--quick` checks exact algebra using archived interval receipts and does not freshly replay interval signs. Install `python-flint` and `sympy`, keeping assertions enabled. The wrapper hash-checks the package and runs in a temporary copy. Only the radial verifier's output-directory handling was changed for portability; frozen mathematical inputs and recorded audit receipts are unchanged.

External comparison: Zubrilina(2026), [Table1](https://arxiv.org/html/2610.35830v1#S4.T1), lattice-only upper4.19948220e-7. Our preceding3.527782891e-7 is internal development history, not the external comparison. No global optimality claim is made.
