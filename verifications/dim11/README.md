# 11D curved-profile cyclic certificate

Strict target: density < 0.106638997. Run `python verify.py` with Python 3.11+ and python-flint 0.9.0. Full exact angular and radial gates are replayed; numerical optimizer feasibility is not used for acceptance.

Parameters: minimum squared separation 11/8, profile h(u)=exp[-8(u/s-1)-16(u/s-1)^2], angular degree at most19, radial degree63, physical negative cutoff X=12, and exact M guard1000001/1000000. The all-dilation theorem requires and proves an exact one-sign-transition condition on Q. No external near-neighbor count is used. See THEOREM.md.

The discovery folder includes all final coefficients, exact-export code, numerical settings, source scaling inputs and root-guided rejected checkpoints. A fresh optimizer run is not implied by proof replay and may not reproduce identical floating coefficients. The primary package is frozen for independent review; that status is separate from primary verification.
