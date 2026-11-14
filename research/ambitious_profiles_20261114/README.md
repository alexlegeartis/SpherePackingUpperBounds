# Ambitious profile campaign, 14 November 2026

Sixteen frozen certificates improve dimensions **22, 25, 32 and 35–47**.
The [dated change record](../../CHANGES_20261114.md) gives exact bounds and a
before/after table. [The shared proof](../../verifications/profile_enrichment_20261114/PROOF.md)
derives every mathematical gate. The registered data are in
`verifications/dimN/refinement_profile20261114/`; this directory retains
discovery and audit provenance.

The selected profile is `h(v)=v^(-a) exp[-b(v-1)]`, `v=|x|^2/s`.
Fifteen dimensions use `a=-96,b=120`; 45D uses `a=-192,b=216`.
Each angular polynomial is reoptimized and each radial search uses the
degree-111 source polynomial together with degree-95 Laguerre corrections.
The source is transformed exactly into the target dimension. The admissibility
condition `a+b>0,b>0` admits these negative values of `a`, with a monotone
logarithmic slope and revised polynomial factors in the infinite-tail proof.
No external neighbor counts enter the sixteen new certificates.

## Frozen proof replay

Only `python-flint==0.9.0` and the Python standard library are needed for
these new proof gates. Do not disable assertions.

```bash
python research/ambitious_profiles_20261114/replay_registered.py --jobs 4
python comparison.py
```

For one result, for example:

```bash
python verifications/dim38/refinement_profile20261114/verify.py
```

The registered wrapper first checks its source hashes and registry, then
runs primary angular/cyclic, primary radial/Fourier, independent angular/cyclic
and independent radial verification. It writes full logs and JSON results.
The implementations use different reconstruction/enclosure formulas and share
FLINT/Arb, including root isolation. No separate backend or external audit is
claimed. The new proof programs were all replayed directly before registration;
the 38D wrapper was additionally run end to end. Changed historical wrappers
were checked through their registry/hash prefixes. This campaign did not rerun
the complete historical aggregate; see `CAMPAIGN_AUDIT.json` for its exact scope.

## Discovery and exact exports

Install [requirements-discovery.txt](requirements-discovery.txt) for the
floating-point search tools. The matrix dependencies are optional for the
selected scalar-profile search. `RUNTIME.json` records the actual environment.

```bash
python research/ambitious_profiles_20261114/search_angular.py 38 --a=-96 --b=120 --rounds 6
python research/ambitious_profiles_20261114/search_radial.py 38 --s 19/4 --a=-96 --b=120 --degree 95 --cut 60 --angular research/ambitious_profiles_20261114/n38_a-96_b120_ang19/angular.json
python research/ambitious_profiles_20261114/continue_radial.py research/ambitious_profiles_20261114/n38_d95_s19_4_a-96_b120_q1_c60/candidate.json --rounds 20
```

These commands modify discovery files. Use a fresh checkout or a new `--tag`
for new searches, or restore the checkpoint archive before resuming its history.
Floating-point optimizers need not reproduce identical coordinates across
versions; the frozen witness verification is exact and independent of solver
success. `certified:false` in a raw candidate means its numerical objective
is not a density theorem, even when its final exported witness later passes.

The selected `candidate.json`, `cuts.json`, angular source and successful
`radial_certificate.json` remain directly readable. To check that the frozen
polynomials, guarded origin and objective are exactly those exported from the
retained binary64 coordinates:

```bash
python research/ambitious_profiles_20261114/audit_exports.py 22 25 32 35 36 37 38 39 40 41 42 43 44 45 46 47
```

`EXACT_EXPORT_AUDIT.json` records that check. `CAMPAIGN_RESULTS.json` distinguishes
improvements over preceding main `8a080d6` from comparisons with the external
displayed table. `sources/COMPARISON_SNAPSHOT.json` extracts the unrestricted
packing rows from the earlier same-day archived response. A follow-up request
returned HTTP 502. This source check does not establish exhaustive priority.

## Retained failures and broader experiments

`DISCOVERY_INDEX.json` lists the numerical candidate families. Full continuation
rounds retain the root sign islands, actual positive pair counterexamples,
added cuts and final successful checks. To keep repeated large rational
polynomials out of thousands of loose files, these checkpoints and unselected
experiments are in `discovery_checkpoints.tar.xz`. Every member was read back
and SHA256-checked before compression replaced the loose copy; the complete
index is `ARCHIVE_MANIFEST.json`. Extract in this directory with:

```bash
tar -xJf discovery_checkpoints.tar.xz
```

Exploration included pure degree-127 radial spaces, sharper scalar profiles,
two/three/four-feature harmonic PSD matrix kernels, and Gram-ray LP column
generation. Matrix candidates produced tempting numerical objectives, but
some had genuine positive triangle counterexamples and others had solver
failures or incomplete continuum proofs. None contributes to a registered
density result. `interrupted_discovery.json` identifies deliberately stopped
exploratory runs. Partial PSD, cap, tail or sampled checks are not complete
triangle certificates.

The first 44D angular witness failed exact Sturm verification: it had a
positive island about `1e-5` wide near `t=-0.97073059`, with local maximum
about `0.00515768`, between sampled constraints. Additional angular cutting
planes repaired it before registration, followed by fresh radial export and
all four frozen proof paths. `dim44_rejected_angular/` retains the rejected
data and log; `audit_campaign.py` checks that the real failure is still
rejected. The repaired angular/radial discovery folders have `_repaired`
suffixes. The supplementary repair archive retains their continuation rounds.
Its member hashes are in `ARCHIVE_MANIFEST_repair_dim44.json`; restore it
with `tar -xJf dim44_repair_checkpoints.tar.xz` after the main archive.
