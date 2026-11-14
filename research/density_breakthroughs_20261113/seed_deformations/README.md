# Seed deformations and the 44-dimensional periodic contact cone

**Outcome: no new sphere-packing density record.** This directory records concrete searches, exact rejected proposals, and a numerical obstruction to local periodic improvement of a published-record-density 44D realization. It does not claim that any searched family exhausts all packings.

## Main reusable result: full periodic cone in published P48n

`periodic44_p48n/periodic44_eutaxy_sample.npz` and `audit_periodic44.py` are the compact, self-contained witness and checker. The ambient Gram equals the archived catalogue `CQ48a.html` **entry for entry**, so its minimum 6 is inherited from the published P48n lattice; no unproved isometry to another Gram is used.

The exact rank 4 normal Gram is

```
 6 -1 -1 -2
-1  6 -2 -2
-1 -2  6 -1
-2 -2 -1  6
```

It has determinant 605. The five normal-coordinate labels are 0, e1, e2, e3, e4. The five longest differences form a cycle; each has squared normal-projection length 16/55. Thus norm 6 mother-lattice vectors on these differences give projected contacts of squared length 314/55. The usual antipode bound gives center density

`5 * (157/110)^22 / sqrt(605) = 509.618982642864...`,

matching the Chen–Hu–Li–Wang–Wu 44D record parameters. The checker computes the exact Smith normal form of the four normal generators and verifies invariant factors 1,1,1,1. Consequently the normal lattice is primitive, the five coset labels are realized, and the period determinant is605.

The screen allows **all 990 symmetric metric variables and 176 relative coset-displacement variables**, rather than only the six ambient flex parameters. The selected 10000 exact contacts give a nonnegative numerical eutaxy representation of the metric trace functional, with full 1166 numerical rank on positive-weight support. The detailed residuals and thresholds are in `periodic44_eutaxy_audit.json`. This is strong numerical evidence against first-order local improvement. The integer contact identities are exact; the eutaxy weights and rank are floating-point results, not a formal local-optimality theorem.

The full pool of 55964 active contacts was selected from the independently prepared 300457-vector P48n pool in `../gluing/P48n_glue_pool.npz`. The compact checker does **not** need that full pool. Its dependency on the published minimum is explicit.

## Other experiments

| Experiment | Result | Scope |
|---|---|---|
| Q32a contact screen | Rank 528/528 on a sample of its norm 6 shell | Numerical rigidity of the sampled lattice contact constraints |
| 48D3-modular seed | Six integer Gram directions preserve 12485 exact norm 8 contact rays | Exact annihilation identities; nullity computed numerically |
| 48D finite necessary-cut polytope |546 numerical vertices; largest conditional center density 16834.1121957 | Numerical relaxation only; reproduces the known48D density value |
| CTturyn1/3/4/5 | Contact ranks 505,638,494,553 out of 666 | Floating shell enumeration/counts and numerical ranks |
| CTturyn1 three-scale proposals | Exact vectors of squared norms 6,13/2,5 reject three proposed minimum 8 metrics | Exact counterexamples |
| 44D five-point normal-Gram scan | Recovered determinant 605, beta16/55 | Bounded numerical finite search |
| 44D larger dual clusters | A tempting conditional density 561.9449 was rejected by an exact norm 2 normal-lattice vector | No valid improvement; rejection retained explicitly |

For the48D family, write `Q(x)=H+sum(x_i D_i)/206`. The endpoint

`x=(-2,4,6,4,2,28)/3`

has `Q=(4/3)P`, where P is exactly even integral with determinant 1 and has exhibited norm 6 vectors. Its minimum 6 and an isometry to catalogue P48n were **not** independently established. An attempted norm≤4 exclusion hit a billion-node floating search cap; it was not a certificate. All published-P48n 44D conclusions above use a separate, exact catalogue Gram and avoid this unresolved premise.

The older `dim48mod3/periodic44_eutaxy_sample.npz` is a second numerical contact-cone calculation on this recovered endpoint. Prefer `periodic44_p48n/` for the logically complete ambient construction.

## Reproduce the compact checks

From the repository root, with NumPy, SciPy and SymPy installed:

```bash
python research/density_breakthroughs_20261113/seed_deformations/audit_seed_results.py
python research/density_breakthroughs_20261113/seed_deformations/audit_periodic44.py research/density_breakthroughs_20261113/seed_deformations/periodic44_p48n
```

For expensive discovery, run `../build_tools.sh` and install fpylll/cysignals in the chosen environment. The scripts locate executables in the sibling `../bin/` directory. To reconstruct the published-P48n screen from the larger shared contact pool:

```bash
python research/density_breakthroughs_20261113/seed_deformations/periodic44_p48n_setup.py
python research/density_breakthroughs_20261113/seed_deformations/periodic44_adaptive.py research/density_breakthroughs_20261113/seed_deformations/periodic44_p48n
```

To reproduce seed reductions, `prepare_sources.py` and `bkz_sources.py` take catalogue names as positional arguments. All congruences use exact integer transforms after numerical basis selection. C++ shell routines use floating pruning and exact integer leaf checks: their completed or capped runs are **discovery**, not minimum certificates. `metric_descent.py` now reconstructs rational parameters and checks stored exact cuts before enumeration; its older rounded-parameter caches are excluded.

## Retention and sources

Retain the archived catalogue pages, Python/C++ source, compact NPZ inputs/witnesses, and JSON reports. Every retained artifact is under 1MB; the published-P48n 44D eutaxy sample is about312KB. Large raw shell lists, dense cone matrices, iterative state, and rejected-vector caches are ignored. They account for most of the discovery workspace. The nested upstream clone is ignored; `external_source.json` pins its URL and commit.

Primary sources:

- Nebe–Sloane lattice catalogue: <https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/>; the exact archived pages are in `sources/`.
- Chen, Hu, Li, Wang and Wu, *New Sphere Packings from the Antipode Construction*: <https://arxiv.org/abs/2505.02394>.
- Dorofeev–Sun–Wang antipode constructions and machine-readable data: <https://github.com/wcgbg/sphere-packing>, pinned in `external_source.json`.

No existing construction, record registry or certificate was modified by this work.
