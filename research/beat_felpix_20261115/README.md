# Discovery behind the 4D improvement

The [certified result](../../verifications/dim4/refinement_energy20261115/)
is **density <0.629520189662**, below the pinned Felpix-Studios upper bound
0.629521681653. The base certificate and its authors are explicitly
[attributed](../../verifications/dim4/refinement_energy20261115/NOTICE.md).

Two ideas produced the improvement: replace the quartic far-pair energy
by a general degree-16 polynomial, then allow five negative Gegenbauer
coefficients by extracting exact credits from angular Gram energies that
the count census discarded. Neither changes the underlying kernel.

`checkpoint.json` preserves the sampled LP solution. `check_export.py`
reconstructs the frozen rational coefficients using decimal quantization,
scale 999/1000, and zeroing the negligible negative coefficients above
degree 9. `extract_credits.py` independently reconstructs the exact range
solutions from the upstream Gram rays. Both scripts assert agreement with
the published inputs and leave them untouched:

```sh
python research/beat_felpix_20261115/check_export.py
python research/beat_felpix_20261115/extract_credits.py
```

The optional numerical search needs SciPy in addition to the proof
dependencies. For example:

```sh
python -m pip install -r research/beat_felpix_20261115/requirements.txt
python research/beat_felpix_20261115/search_energy.py 1 --full --credits
```

The search writes caches and candidates in its ignored `scratch/` folder.
Cached radial quantities accelerate discovery only; the proof programs
run the original interval engine without these modifications. A solver
or sampled result is never a certificate. Different solver versions may
choose different nearby solutions; the frozen export is deterministic.

`EXPERIMENTS.json` includes rejected attempts. In particular, the
larger-radius polynomial passed the sampled constraints but failed the
full interval check near a one-sided radial boundary. The final theorem
retains the original larger-radius quartic coefficient. Non-finite
endpoint samples are omitted from numerical discovery; every closed
endpoint is still covered by the exact and interval proof replays.

A further **uncertified** direction is preserved in
`conditional_checkpoint_UNCERTIFIED.json`: groups (N,1) have no points in
(R1,R5], so their Gram forms can be restricted to the 21 active radial
coordinates. The sampled fit improves the energy objective, but no
conditional box certificate is included or used. The unchanged (26,5)
profile would become the next bottleneck. This is a research lead, not a
stronger registered density claim.
