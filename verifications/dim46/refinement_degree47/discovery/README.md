# Portable discovery/export provenance for dimension 46

All needed source data are local. `inputs/adaptive_cyclic40_radial111.json` is the exact source48 seed. `cyclic_angular46_exact.json` is the unchanged final angular polynomial. The baseline degree31 candidate, exact radial certificate, old root cuts, all new rounds, exact root-isolation failures and rational sign counterexamples are retained.

## Safe reproduction

Work in a copy of this directory to preserve the frozen history. The first attempt used the old root cuts and reference-degree31 row normalization. Reset `radial_perturb46_degree47_extra.json` from `initial_cuts46.json`, and archive the current degree47 rejection file in the working copy before reproducing the initial search. Run `python radial_search.py --dimension 46 --degree 47`, followed by `python audit_radial.py radial_perturb46_degree47.json`. Subsequent searches incorporate the latest exact-root rejection points.

The initial stage allowed one solve and four root-guided repairs. The continuation allowed six more, but stopped as soon as the full radial audit passed. Continuation searches used `--control-safe-margin`. The saved final coefficient vector, not rerunning a floating solver, defines the certificate.

## Control-padding caveat

The mathematical degree31 family embeds exactly by padding the Laguerre coefficient vector with sixteen zeros. Reference-degree31 row scaling preserves the old normalization instead of allowing new columns to change the artificial margin. The objective and normalization identities were checked numerically.

The sampled LP still used an artificial negative RHS of −1e−7. The saved47D control missed that margin by about6.5e−9 initially despite all sampled signs remaining negative; later root cuts increased the mismatch to about3.0e−8. The46D mismatch was about1.1e−9, within the declared numerical tolerance. These are numerical margin effects, not failures of the certified old polynomials or of algebraic span inclusion.

During continuation, each RHS is max(−1e−7, the padded control row value+1e−12), and is checked strictly negative. This explicitly admits the benchmark coefficient vector. Final exact root, interval pair, and full-package proofs do not rely on this numerical margin, optimizer status, or sampled feasibility.

The audit multiplies M by exactly1000001/1000000 for rigorous pair slack. Reported certified objectives use this inflated M and exact Fourier normalization.
