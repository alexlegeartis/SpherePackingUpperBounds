# Six-dimensional discovery and export

This is the minimal self-contained warm-start search and prerequisite-export
companion to the frozen six-dimensional verification package. It has no
repository, environment, network or absolute-path dependencies after installing
the listed Python packages. It does not replace the proof package.

The frozen numerical output produced the verified strict density upper bound
**0.410189987**. A new solver run is numerical research only: its output may
differ, and a small objective or a successful solver status is not a proof.
Full exact checking is mandatory. Published average-neighbor-count bounds
remain external mathematical inputs, as explained in `METHOD.md` and the
verification package's proof.

## Files and provenance

- `search_exact_guided_v2.py` is the unchanged successful numerical search
  source. It freezes a degree-31 radial polynomial and optimizes degree-four
  Chebyshev radial kernels, angular degrees 0 through 12, and four nonnegative
  count weights. It adds pair and triangle cuts, entire antipodal Bernstein
  constraints, two exact-failure-guided local far boxes, critical angular
  fibers, and a proximal objective.
- `seed.json` contains the complete preceding checkpoint, fixed radial
  polynomial and 736 separation cuts. This seed itself is an invalid
  intermediate candidate, not a certified bound. Its inclusion makes the final
  search stage independently runnable without replaying the exploratory chain.
- `candidate_frozen.json` is the successful search's preserved output, including
  756 cuts and two rounds of diagnostics. Its unsaturated numerical density objective is
  0.4103164518985832; Gaussian repair and exact origin Schur saturation yield
  the final certified result. The snapshot retains `certified: false` because
  certification is established separately by the exact witness and verifier.
- `export_cheb.py` and `gaussian_guards.py` deterministically export exact
  rational matrices, Fourier polynomials, Gaussian guards and Schur arithmetic.
  Export alone does not check complete pair/far/cyclic domains.
- `provenance.json` records hashes, parameters, versions and the original
  successful run's history. `SHA256SUMS` binds this research package.

The candidate SHA-256 is
`797128e7ada65976ce25926ee84f626a9607208013bfb420243215df2f871606`.
The original successful solve took about 199 seconds on its test machine and
stopped after two rounds. Timing and floating-point solver output are not
guaranteed reproducible across machines, BLAS versions or solver releases.

## Run discovery

Python 3.12 was used. Install pinned discovery/export dependencies:

    python -m pip install -r requirements.txt
    python run_search.py --output candidate_new.json --rounds 5

The wrapper explicitly supplies the successful six-dimensional configuration
and resolves its default seed relative to this folder, so it also works from a
different working directory. It refuses to overwrite existing output.

This is a bounded restart from the supplied checkpoint, not a claim that the
entire original exploratory search is reproduced from an empty initialization.
The original numerical solver reported `optimal_inaccurate`; that warning does
not weaken the separate exact verification of the frozen output.

## Export the frozen checkpoint

    python export_and_check.py --output exported_frozen

This creates `witness.json`, `witness_psd_exact.json`, `prerequisites.json` and
an export report in a new directory. It must not be described as full proof
verification. To export a new numerical result, add
`--candidate candidate_new.json`.

## Export and fully check

Supply the path to the separately published frozen verification package:

    python export_and_check.py --output replay_frozen --verification-package ../../verifications/dim6

The command copies only the verifier programs into the new output directory,
leaves the frozen package unchanged, recomputes every exact condition and
enforces its strict upper decimal `0.410189987`. A new candidate above that
decimal will fail even if all its geometric inequalities pass. The proof
package can also be run on its own and needs no search dependency or solver.

Never run the exporter or verifiers with `python -O`. Export and full replay
from the frozen checkpoint were tested from another working directory. The
exported canonical witness and prerequisite data matched the frozen proof
package exactly as parsed JSON. A fresh two-round warm-start search also
reproduced every coefficient, matrix entry, count weight and separation cut
exactly in the tested environment. This does not guarantee bitwise numerical
reproduction in other environments. See `TEST_RESULTS.json` for the checks.

Repository integration preserves the portable primary verifier as `verify_primary.py`; the export wrapper selects it instead of the registry-checking wrapper when present. The standalone export replay checks the primary full proof, while the registered package wrapper additionally runs the independent audits.
