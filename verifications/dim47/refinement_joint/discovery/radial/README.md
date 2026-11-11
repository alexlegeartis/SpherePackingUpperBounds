# Radial discovery and exact export

The final numerical model uses a nonnegative multiple of the saved normalized degree-111 source plus a degree-31 dimension-47 Laguerre correction, with Fourier origin normalized to 1. The shared new angular cost is frozen in cyclic_angular40_exact.json.

The first sampled candidates failed exact global sign checks. Exact root-isolation counterexamples were added repeatedly. refinement_history.json records each stage; the initial rejected data and final checkpoint are retained. M receives the exact multiplicative guard 1000001/1000000 before the pair proof.

radial_search.py --dimension 47 --degree 31 reproduces the final sampled model and its saved root cuts. It requires NumPy, SciPy and mpmath. audit_radial.py radial_perturb47_degree31.json reproduces rational export and radial verification with python-flint. The copied source file makes these two scripts portable. A new optimizer run need not reproduce identical floating-point coefficients. The exact export has been replayed; no new optimizer replay is implied.

The complete certificate is verified from the parent folder without any optimizer or discovery dependencies. Do not use the source-candidate values as bounds.
