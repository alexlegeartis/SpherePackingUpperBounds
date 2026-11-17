# Verified 82D all-deep binary-plane A2 gluing

Ordinary density ≥ 9.17033181889403e-18, minimum 8, Gram determinant 3, center-density square 2^82/3. This improves the preceding published mixed-depth 82D result by 7.0436048222%. Against the derived external MW80 × A2 example the improvement is 260.77894495%, not a strongest-prior claim.

Three exact parent automorphisms transfer already audited depth≥6 half-cosets to three distinct binary classes whose sum is zero. They form a binary plane. Glue its two generators to the A2 basis with vertical Gram [[2,-1],[-1,2]]. Every nonzero parity class has vertical minimum≥2 and horizontal minimum≥6. Nonzero even coefficient vectors have vertical minimum≥8. Parent vectors attain8. Schur coarea gives determinant3.

The finite arithmetic verifier reconstructs each full80-dimensional composite action and checks integrality, determinant±1, Gram preservation, exact source lift and binary class. Thus the group enumeration used in discovery is not needed for proof.

Run python run_all.py for all exact checks and archived search integrity. --full additionally replays inherited primary and independent exhaustive searches and controls. No new exhaustive absence is needed for this improvement. Requirements: Python3, numpy, python-flint and g++, assertions enabled. All replays use temporary copies.

The original81/82 source proof packages are retained intact, including historical pending labels, resolved by their complete independent audits. Parent minimum8 inherits Bachoc–Nebe Theorem4.7 with exact constructive hypotheses checked.
