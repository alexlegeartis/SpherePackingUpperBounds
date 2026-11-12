# Signed-profile 9D discovery checkpoint

Numerical discovery, not a density proof. The corresponding fully verified
strict bound is 0.190941549; see ../../verifications/dim9/.

This directory preserves the frozen step-11 angular/profile/radial inputs,
radial constraint rows, sensitivity data, and joint optimizer. The frozen
fit9.json and refit9.json agree exactly with the verification package inputs.
The monotone basis uses nonnegative initial value, jumps and interval slopes
of g(v)=v*phi'(v). It spans the admissible signed-hinge profiles.

Python dependencies: NumPy, SciPy and python-flint. In a disposable copy:

- `python basis.py` tests basis identities.
- `python joint.py` proposes a new joint angular/profile step using frozen
  radial LP sensitivities; output is not certified.
- `python refit_count.py 9` refits and root-cuts a radial candidate for the
  current fit9.json and inputs/angular_dim9.json. Floating output varies with
  solver/library versions and is not a rigorous result.
- `knots.py` is an exploratory knot-location optimizer, not used in this
  registered certificate.

No claim is made that rerunning optimization reproduces identical coefficients.
The deterministic exact export and complete primary/independent certificate
replay are in ../../verifications/dim9/. Any new candidate must pass those
checks, including all-scale monotonicity and the full cyclic domain.
