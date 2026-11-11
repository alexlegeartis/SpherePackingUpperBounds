# 15D curved-profile cyclic certificate

Strict target: density < 0.033671603. Run `python verify.py` with Python3.11+ and python-flint0.9.0. Numerical discovery is not used for acceptance.

This kernel uses quartic radial damping, not a single Gaussian rate. The full angular proof includes an exact one-sign-transition condition and an all-dilation lemma. The radial polynomial has degree63, with complete coupled pair verification through X=12 and exact physical negativity thereafter. The normalization is exact and the origin-cost guard is1000001/1000000. See THEOREM.md for all geometric and finite-sum assumptions.

The discovery directory preserves the initial sampled angular result, final radial coefficients, root-guided rejected checkpoints and fixed inputs. Optimizer replay is optional, may vary numerically, and is not implied by proof replay. `screen.py` first runs a=12,b=16 and two angular separation rounds; `search_radial.py` and the saved cut12 history then produced the final polynomial. Earlier cut40 radial trials failed complete pair signs and are retained as rejected checkpoints.

The proof package is frozen for independent review. Do not interpret primary verification as an independent audit.
