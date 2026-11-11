# Portable discovery and exact export for dimensions 22 and 23

This companion preserves the search needed to regenerate the Gaussian cyclic witnesses. It is separate from the frozen proof packages and does not modify them. No network or other repository checkout is needed after obtaining this directory.

## Tested environment and dependencies

Python 3.12.14; python-flint 0.9.0; NumPy 2.3.5; SciPy 1.17.0; mpmath 1.3.0. The search uses scipy.optimize.linprog with bundled HiGHS and L-BFGS-B. Set OPENBLAS_NUM_THREADS=1. Exact export/replay requires only python-flint; no optimization solver is used during frozen export replay.

    python -m pip install -r requirements.txt
    python replay_frozen_exports.py
    python verify_known_lower.py

For export-only use requirements-export.txt. Fresh end-to-end tests for both dimensions also passed after seven radial solves each; see fresh_search_replay.json and fresh_search22/23.log.

The first replay works in an isolated temporary directory, rebuilds both angular polynomials and both full radial certificates from the frozen numerical checkpoints, performs complete angular/domain and radial checks, and compares all four exported files byte-for-byte with the frozen proof inputs. EXPECTED_EXPORTS.json pins those hashes. INPUT_MANIFEST.json pins every bundled source input.

## Fresh bounded discovery

    python run_search.py --dimension 22 --output /tmp/new-packing22
    python run_search.py --dimension 23 --output /tmp/new-packing23

The output directory must not exist. The orchestrator runs at most eight angular separation rounds and ten exact-root-guided radial repair rounds by default; each subprocess has a 180-second timeout. The limit is adjustable with --repair-rounds, up to 20. A failed check or exhausted budget is not a result. Floating optimizer outputs can differ across systems; only the subsequent exact global checks establish feasibility. No claim of optimality is made.

The lower-level scripts can also be run directly in a working copy:

    python search_angular.py --dimension 22 --degree 19 --alpha 40 --rounds 8
    python export_angular_n.py 22
    python search_robust.py --dimension 22 --degree 31 --r2 2.75
    python dim22/audit_radial.py radial_scaled22_degree31_r2p75.json

For 23D use dimension 23 and r2=2.875. Alternate search and radial audit when forbidden roots are reported. Rejection roots and midpoints feed the next LP. The archived final *_rejection.json files are from the preceding failed iteration and are discovery diagnostics, not proof inputs.

## Mathematical inputs and provenance

- inputs/angular_scaling48.json supplies positive column scales for angular numerical conditioning only.
- inputs/source44/certificate.json supplies the exact degree 111 physical polynomial and normalization for the radial search basis.
- inputs/source44/continuum_three_point/cyclic_angular40_exact.json supplies the old angular fallback data used by the generic search, although fresh target angular data are preferred when present.
- Final numerical checkpoints and repair-point lists are at the companion root. Provenance logs are under provenance/.
- templates/ contains the exact radial exporter and complete angular/domain verification scripts copied into each generated dimension folder.

Only path references were changed from the original search/export scripts. The exact export math is unchanged: angular rationalization subtracts 1/1000 from the constant harmonic and renormalizes it to 1; radial rounding uses denominators at most 10^14 and increases M by 501/500. The final angular polynomials have degree 14, within the searched degree 19 space. The physical/Fourier polynomials have degree 111. Both source polynomials are treated as exact discovery data, not as unverified mathematical assumptions in the final certificates.

Large numerical .npz caches are omitted and regenerate automatically. The proof packages themselves require none of these search inputs.

## Known lower bounds and gap comparisons

KNOWN_LOWER_BOUNDS.json records Henry Cohn's indexed table formulas and independent 256-bit arithmetic. The table calls the following quantities center densities for unit-radius spheres:

- 22D: 2^(-23) * 3^(-21/2) * 11^11.
- 23D: 2^(-34) * 3^(-12) * 23^(23/2).

Multiply by pi^(n/2)/Gamma(n/2+1) to obtain the known packing densities 0.0024510340441211… and 0.0019071945388448… respectively. The source is https://cohn.mit.edu/sphere-packing/, retrieved through its search index on 2026-11-11; the index was four months old and the live page returned 502. These are sourced published constructions, not construction-coordinate replays or an exhaustive current-best claim.

The computed reductions of the pinned upper-to-known-packing gaps are approximately 5.6913% and 5.3713%. These differ from reductions in the upper bounds themselves, approximately 1.9474% and 1.5550%. Full enclosures are in the JSON. The companion does not alter the proof registries or their comparator policy.
