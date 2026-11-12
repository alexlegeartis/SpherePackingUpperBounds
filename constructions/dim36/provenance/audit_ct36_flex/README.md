# Independent CT36 metric-flex audit

No new density bound is certified by this audit. The input contact list is incomplete.

`audit.py` independently verifies the integral symmetric direction D and G=G4/4, exact positive LDL decomposition, determinant 3^18(5/4)^12, and all 35,864 supplied vectors having G-norm 8 and D-norm zero. For A=G^-1 D, exact computation gives rank24, trace0, A^3=(16/15)A, and characteristic polynomial lambda^12(lambda^2-16/15)^12. Hence G+tD is positive definite precisely for |t|<sqrt(15)/4 and its determinant is det(G)(1-16t²/15)^12.

**Conditional only:** if its minimum remains8, its density relative to KP36 is (12288/15625)/(1-16t²/15)^6. This exceeds1 when |t|>0.19182618074183197. At t=1/5 the conditional ratio is exactly2187000000000000/2140753641621841, approximately1.02160284.

`neighbors.py` and `pairs.py` inspect supplied contacts perturbed by one/two unit coordinates. They do not certify absence of shorter vectors globally. The first contact event among the tested two-coordinate neighbors is |t|=1/4, with a norm17/2 vector having D-norm-2, saved in pairs.json. All this says is that this restricted test did not reject t=1/5. `old_witnesses.py` checks the old min6, newshort, kernel and representative witnesses after exact unimodular coordinate conversion; none rejects either sign of1/5.

Required remaining step: full minimum certification or a single rigorously checked obstruction. The saved shell was capped at100million nodes and is not a global minimum proof, even at t=0.

## Full minimum verifier audit

The subsequent `research/certify_ct36_minimum/interval_enum.cpp` run reports COMPLETE1,108910425 nodes and one leaf (zero). `audit_interval_setup.py` independently verified Q=20G(1/5), Q36.txt, unimodular transformation and all exact-LDL interval enclosures. The recursive state is a lower bound on the already assigned LDL norm contribution. Outward centers, lower d and outward radius prevent pruning any norm≤159 integer vector. The exact integral Gram has even diagonal, so excluding norms≤159 proves minimum≥160. Existing contacts prove equality.

The audit found no mathematical flaw, assuming standard IEEE floating behavior with no unsafe fast-math transformations. Crude coordinate envelopes remain below2^38, and exact leaf sums are far below signed128-bit overflow. Input/source hashes were saved to interval_audit_sha256.txt. This minimum certification does not itself audit the contemporary literature baseline or assert novelty.
