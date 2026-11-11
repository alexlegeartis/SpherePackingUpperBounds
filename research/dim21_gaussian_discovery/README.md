# Fresh 21D Gaussian cyclic search

A fresh dimension 21 angular optimization and joint radial Fourier search produced the exact witness in ../../verifications/dim21/refinement_gaussian/. Complete producer verification proves density <0.004985977000, versus previous repository 0.005149549312 and pinned published 0.0051596603948176. Independent verification has passed; see the ../../verifications/dim21/refinement_gaussian/independent/ audit package.

This directory is self-contained for discovery: inputs/ contains the exact 44D radial source and 48D angular scaling input; templates/ contains the exact exporter and angular/domain checkers. Search uses Python 3.12 with numpy 2.3.5, scipy 1.17.0, mpmath 1.3.0 and python-flint 0.9.0. Set OPENBLAS_NUM_THREADS=1.

    python search_angular.py --dimension 21 --degree 19 --alpha 40 --rounds 8
    python export_angular_n.py 21
    python dim21/independent_cyclic40_angular.py
    python dim21/independent_cyclic40_domain.py
    python search_robust.py --dimension 21 --degree 31 --r2 2.625
    python dim21/audit_radial.py radial_scaled21_degree31_r2p625.json

Run from this directory, or use complete paths. Repeat search_robust.py followed by audit_radial.py when forbidden roots are reported; the rejection file supplies the next separation rows. The original pass required seven radial solves, with the first six rejected. The archived last rejection file is stale discovery information, not a failed final proof.

Final angular degree 13 within degree 19 search; radial degree 111. The frozen proof package needs no optimizer or these source files. Run `python replay_frozen_exports.py` to re-export the frozen numerical checkpoints in a temporary directory, replay every angular/radial proof and compare both exact exported files byte for byte. This does not rerun the optimizer. Numerical search logs/candidates are not independently sufficient for any density claim.
