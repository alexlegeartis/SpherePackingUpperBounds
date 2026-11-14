"""Write the dated, source-attributed change record from frozen proof results."""
import json,platform,sys,hashlib
from pathlib import Path
from decimal import Decimal as D,localcontext
from importlib.metadata import version
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
dimensions=[22,25,32,*range(35,48)]
rows=json.loads((ROOT/'RESULTS.json').read_text());facts=[]
snapshot=json.loads((HERE/'sources/COMPARISON_SNAPSHOT.json').read_text())
external={int(r[0]):r[4] for r in snapshot['rows']}
with localcontext() as ctx:
    ctx.prec=80
    for n in dimensions:
        P=ROOT/f'verifications/dim{n}/refinement_profile20261114'
        cfg=json.loads((P/'config.json').read_text());result=json.loads((P/'full_replay_result.json').read_text())
        assert result['pass'] and result['complete_primary_and_independent_replay']
        row=next(r for r in rows if r['dimension']==n);assert row['package']==cfg['package']
        old,new=cfg['previous_repository_strict_upper'],cfg['strict_density_upper_bound']
        assert D(new)<D(old) and D(new)<D(external[n])
        gain=100*(D(old)-D(new))/D(old)
        facts.append({'dimension':n,'previous_main_strict_upper':old,'new_strict_upper':new,
                      'decrease_from_previous_main_percent':format(gain,'.8f'),
                      'live_table_displayed_upper':external[n],
                      'decrease_from_live_displayed_upper_percent':format(100*(D(external[n])-D(new))/D(external[n]),'.8f'),
                      'package':cfg['package'],'a_exact':cfg['profile_a_exact'],'b_exact':cfg['profile_b_exact']})
table='| Dimension | Preceding main upper | New strict upper | Further decrease | Certificate |\n| --- | --- | --- | --- | --- |\n'
table+=''.join(f"| {f['dimension']} | {f['previous_main_strict_upper']} | **{f['new_strict_upper']}** | {D(f['decrease_from_previous_main_percent']):.4f}% | [proof]({f['package']}/) |\n" for f in facts)
(HERE/'CAMPAIGN_RESULTS.json').write_text(json.dumps({'base_commit':'8a080d6a00909391a91ecbac097e4004c404c40d','results':facts},indent=2)+'\n')
changes='''# Sixteen stronger sphere-packing upper bounds — 14 November 2026

This campaign improves the unrestricted packing-density upper bound in
dimensions **22, 25, 32 and every dimension 35–47**. The table compares with
the preceding main commit `8a080d6a00909391a91ecbac097e4004c404c40d`, which
already includes the earlier 14 November improvements in 45D–47D. All bounds
are strict rational decimals; smaller values are better.

'''+table+'''
The largest further reductions are about 4.03% in 39D, 3.97% in 38D and
3.54% in 37D. The earlier same-day 45D, 46D and 47D values are superseded;
their exact witnesses and replay evidence remain available.

## What changed

The search now optimizes signed logarithmic radial profiles
`h(v)=v^(-a)*exp[-b(v-1)]`, where `v=|x|^2/s`, `a+b>0` and `b>0`.
Fifteen selected witnesses use `a=-96,b=120`; 45D uses `a=-192,b=216`.
Allowing negative `a` changes the contact slope and radial curvature while
preserving monotonicity and the all-dilation cyclic argument. Each profile
has a newly optimized dimension-specific angular polynomial. The physical
and Fourier polynomials have degree 111, obtained from a retained seed with
degree-95 Laguerre corrections. This expands both the radial profile and
the correction space beyond the preceding Gaussian continuations.

The [shared proof](verifications/profile_enrichment_20261114/PROOF.md)
includes a smooth global PSD feature extension, the exact origin Schur
identity, and the revised infinite-tail bounds needed when `a` is negative.
It uses the existing cyclic Cohn–de Laat–Salmon continuum framework. These
sixteen packages require no external neighbor counts.

`RESULTS.json`, both result tables, method labels and `run_all.py` include
the new packages. Historical wrappers accept a stronger registered successor
while continuing to check their original frozen witnesses. Their mathematical
inputs and proof conditions are unchanged. Discovery checkpoints, including
rejected sign islands and positive pair counterexamples, are retained apart
from the certificates. Matrix-kernel experiments that did not pass the full
triangle condition remain unregistered research candidates.

The first sampled 44D angular candidate was rejected by the exact Sturm
gate for a tiny positive interval near `t=-0.97073059`. Additional cutting
planes repaired it, followed by fresh radial export and complete replay.
The original data and failure log remain a real rejection control.

## Validation and reproduction

Every new package passed all four proof programs: primary angular/cyclic,
primary radial/Fourier, independent angular/cyclic, and independent radial.
Checks cover exact Fourier reconstruction, all roots on both unbounded
sign rays, the complete pair interval, the full triangle domain, all
dilations, infinite tails, origin PSD and rational strict density rounding.
The independent paths use different reconstructions and enclosures but share
FLINT/Arb and its root isolation; no separate arithmetic backend or external
audit is claimed.

Exact export replay reconstructs both rational polynomials, guarded origin
coefficient and objective from the retained binary64 discovery coordinates.
Registry/hash gates and rejection controls check tampered inputs, incorrect
Fourier data, negative harmonic coefficients, an inadmissible profile slope,
unrelated registry paths and disabled assertions. Generated comparisons agree
at 70 and 100 decimal digits. Full proof logs and result files are in each
new package; aggregate evidence is in
[`research/ambitious_profiles_20261114/`](research/ambitious_profiles_20261114/).

The new proof programs were replayed directly before registration. Every
changed historical wrapper's registry/hash prefix was checked; the full
historical `run_all.py` aggregate was not rerun in this campaign. The 38D
registered wrapper was also exercised end to end. The earlier historical
replay scope remains in the preceding campaign's `REPLAY_STATUS.json`.

With Python and `python-flint==0.9.0` installed, assertions enabled:

```bash
python research/ambitious_profiles_20261114/replay_registered.py --jobs 4
python research/ambitious_profiles_20261114/audit_exports.py 22 25 32 35 36 37 38 39 40 41 42 43 44 45 46 47
python comparison.py
```

For one dimension, run its `verify.py` wrapper. The first command above
replays all sixteen complete registered packages. `python run_all.py --jobs 4`
also replays the retained historical chains, so it takes substantially longer.

## Comparison scope

All sixteen new values are below the corresponding entries in the first,
unrestricted-packing table of [Henry Cohn's page](https://cohn.mit.edu/sphere-packing/)
in the archived HTTP-200 response retrieved on 14 November. A later fetch
returned HTTP 502, so the comparison uses that earlier same-day response.
The [extracted rows and source hash](research/ambitious_profiles_20261114/sources/COMPARISON_SNAPSHOT.json)
are retained. External witnesses were not replayed and this is not an
exhaustive literature or worldwide-priority search. The README's main table
continues to use its attributed historical published comparators; its
percentages therefore differ from the further reductions above.
'''
(ROOT/'CHANGES_20261114.md').write_text(changes)
readme=ROOT/'README.md';text=readme.read_text();start=text.index('> [!IMPORTANT]');end=text.index('> [!IMPORTANT]',start+1)
callout='''> [!IMPORTANT]
> **Sixteen stronger upper bounds — 14 November 2026.** New signed logarithmic profiles and degree-95 radial enrichment improve **22D, 25D, 32D and every dimension 35D–47D**. The largest further reductions from the preceding main are **4.03% in 39D**, **3.97% in 38D** and **3.54% in 37D**. All sixteen complete primary and independent proof chains pass. [Exact bounds, before/after table and what changed](CHANGES_20261114.md) · [Proof and verification scope](verifications/profile_enrichment_20261114/PROOF.md). The earlier same-day 45D–47D results are superseded and retained as historical witnesses.

'''
readme.write_text(text[:start]+callout+text[end:])
methods=ROOT/'METHODS_AND_EXPERIMENTS.md';text=methods.read_text()
text+='''
## Signed logarithmic profiles and degree-95 enrichment, 14 November 2026

The sixteen packages in dimensions 22, 25, 32 and 35–47 replace the selected
Gaussian witnesses with `h(v)=v^(-a) exp[-b(v-1)]`. The admissibility condition
is `a+b>0,b>0`, allowing negative `a`. The contact slope and curvature can
therefore change independently. Selected parameters are `(-96,120)`, except
45D's `(-192,216)`. Dimension-specific angular kernels and degree-95 Laguerre
corrections are reoptimized; the full radial polynomials have degree 111.

The powers of the dilation cancel in each cyclic profile ratio, leaving two
ordered exponential rates. The exact single-transition angular condition
therefore extends the compact proof to every dilation. Infinite-tail bounds
include the additional polynomial factors `4^max(0,-a)` and
`U^(2*max(0,-a))`. Complete compact subdivision covers `t in [-1,1/2]`,
including negative cosines, and the independent checker uses monotone ratio
endpoint bounds. Complete ray and pair proofs, PSD origin normalization and
rational Machin rounding are replayed by two implementations sharing FLINT/Arb.

The [dated change record](CHANGES_20261114.md) gives every exact bound and
comparison, the [proof](verifications/profile_enrichment_20261114/PROOF.md)
derives all gates, and the [research archive](research/ambitious_profiles_20261114/)
retains successful and rejected discovery checkpoints. Multiple-feature PSD
matrix searches produced numerical candidates but no complete new bound;
they remain outside the result registry.
'''
methods.write_text(text)
comparisons=ROOT/'COMPARISONS.md';text=comparisons.read_text();heading='# Comparison sources and conventions\n'
addition='''
**Sixteen-profile follow-up, 14 November 2026:** the selected bounds in
22D, 25D, 32D and 35D–47D are now stronger than the preceding main at `8a080d6`.
The [dated change record](CHANGES_20261114.md) gives exact before/after values;
[machine-readable comparisons](research/ambitious_profiles_20261114/CAMPAIGN_RESULTS.json)
also compare each new strict decimal with the live table's displayed upper.
All sixteen improve the first unrestricted-packing table in the earlier
same-day HTTP-200 snapshot. The later retrieval attempt returned HTTP 502;
the archived response, rather than an inferred fresh value, is used.
The earlier 45D–47D continuation described below is superseded by these
new witnesses. Historical published comparators in the generated result
table remain unchanged.
'''
comparisons.write_text(text.replace(heading,heading+addition,1))
runtime={'python':sys.version,'platform':platform.platform(),'packages':{name:version(name) for name in ['python-flint','numpy','scipy','mpmath','cvxpy','clarabel']}}
(HERE/'RUNTIME.json').write_text(json.dumps(runtime,indent=2)+'\n')
print('WROTE README, change record, method/comparison notes, campaign results and runtime')
