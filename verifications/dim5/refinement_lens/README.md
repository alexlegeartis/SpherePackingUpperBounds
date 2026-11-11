# Five-dimensional tail-funded lens certificate

This frozen, independently audited certificate proves unrestricted sphere-packing
density **strictly below 0.509996147** in dimension 5. It uses the elementary
centroid lens caps only. There are no Gaussian corrections and no newer
transverse, Rankin, or spherical-code lens budgets.

## Replay

    python -m pip install -r requirements.txt
    python verify_primary.py
    python independent/verify_all.py

The package is self-contained and may be copied anywhere. The entrypoint resolves
all inputs relative to itself, not the working directory. No solver or network is
needed after installing dependencies. Do not use Python optimization flags.
`verify_primary.py` rejects `python -O` and `python -OO` before reading proof inputs.

Every run recomputes exact source and parameter binding, all 66 positive rational
PSD pivots, direct Fourier identity, global Sturm signs, exact origin Schur
saturation, nonnegative count costs, outward density rounding, step-budget
inequalities, explicit half-open count endpoints, the complete pair domain,
compact cyclic domain and the new lens far domain. Cached success flags do not
replace fresh checks. Compact cyclic initialization is rebuilt by its checker.

Fresh portable producer replay: pair 121 boxes, compact cyclic 115 boxes,
new far 2,879 boxes (11 relaxed, 526 geometry exclusions). Strict decimal:
0.509996147. Density enclosure center:

    0.5099961464057082278011891612604254673934488294265848385588116530194451139022

The frozen exact source SHA-256 is
`235920216afbd10b38be5efd3672e9d94b5ec5af8fc54aca3a63d77d9997721f`.
`source_exact.json` is byte-for-byte unchanged. The sole change to the frozen
`witness.json` is metadata: its absolute `lens_budget` path is replaced by the
relative `step_allowances_005.json`. No mathematical payload changes.

`SHA256SUMS` freezes the shipped inputs, code and documentation; outputs that
contain runtime paths, timings or generated arrays are intentionally excluded.
Run `sha256sum -c SHA256SUMS` before verification. Verification writes reports and
its cyclic initialization in this package directory, which must be writable.

See `PROOF.md` for the actual six-cap theorem, parameter conventions and trusted
components, and `INDEPENDENT_AUDIT.md` for the completed source-reviewed audit.
`HISTORICAL_LENS_THEOREM.md` preserves the research theorem discussion: its
illustrative three-cap table does not define this package's six-cap input.
The standalone independent continuum audit remains a separately preserved
companion, not a claim that every producer checker was independently rewritten.
The numerical warm-start and exact export workflow is in the sibling
`lens_discovery_package`; optimizer output is never itself a certificate.

The count bounds 45, 46, 47, 48 use published spherical SDP inputs, not rerun here:
https://arxiv.org/html/1206.2608 (Table 2). No external code bound enters the lens
caps. Compared with the paper-reported 0.5122645 at
https://arxiv.org/html/2512.04829v2#S3.T1, this outward bound reduces the upper by
0.44280894% and closes 4.82557591% of the gap to the D5 density pi^2/(15 sqrt(2)).
The comparator certificate is not replayed. No optimality, exhaustive literature
priority, or formal proof-assistant certification is asserted.

Packaging adds fail-fast parameter and canonical payload checks to the frozen
far checker, without changing its polynomial construction, allowance selection,
geometry, subdivision, or sign acceptance rules. It also adds a fail-fast
optimization guard to the query helper. The pair and cyclic source files remain
byte-for-byte identical to the established producer checkers.

The repository wrapper checks the registered threshold and frozen source/budget identity, then runs the portable primary proof (`verify_primary.py`) and separate full independent audit. Earlier certified versions remain in sibling packages.

Within the repository, `python verify.py` additionally checks RESULTS.json and the registered package path before running both complete chains. That registry wrapper requires the repository layout; the two commands above are portable standalone replays.
