# Exact origin saturation of cyclic kernels

Each source package proves a full cyclic density bound with its positive-definite degree-zero block A=[[a,b^T],[b,B]]. Exact Fraction elimination proves B positive definite, solves Bv=b and checks alpha=a-b^T v>0. The verifier independently checks every entry of

    A-alpha e0 e0^T = W B W^T,  W=[v^T;I].

The derived block is positive semidefinite, generally singular. The isolated origin feature is nonzero only at displacement zero. Subtraction therefore changes only K(0,0). Every annular pair expression at a nonzero displacement, every compact cyclic sum, every far term, the physical tail and Fourier transform are unchanged. Closed boundaries are unchanged. The same periodic finite-sum cyclic argument remains valid, decreasing the normalized objective by alpha/H(0).

The certificate is the unchanged source witness plus this exact derived transformation. Keep the positive-definite source unchanged for full-domain replay. A strictly-positive-pivot source verifier is not an appropriate test of the derived singular matrix.

The independent Fraction implementation reconstructs the normalized objective from the source polynomials and matrices, checks it equals exact_objective, and encloses the final density at 384-bit Arb precision. Each requested outward decimal is asserted strictly. Exact source hashes prevent substitution. Source theorem, all domain conditions and the original computational assumptions remain inherited, and verify.py explicitly replays them before the algebra.

These ordinary Schur-complement refinements do not add a new density theorem or a worldwide-priority claim. The single-Gaussian implementation here applies only to dimensions 18-20. Repaired 17D uses a changed Fourier normalization and is separately proved in ../dim17/refinement_repaired/repair_checks.py.

A second independent standard-library Fraction implementation, verify_rational.py, reconstructs every saturated residual entry and proves the outward thresholds using rational Machin bounds for pi (squaring odd-dimensional densities). It is invoked after the source hashes and registry checks; it has no Arb dependency. Its source continuum dependencies are still mandatory.
