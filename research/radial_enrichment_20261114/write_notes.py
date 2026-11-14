"""Write dimension-specific proof notes and the bounded-search record."""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal,localcontext
import json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
summary=[]
for n,degree in [(45,47),(46,63),(47,47)]:
    D=ROOT/f'verifications/dim{n}/refinement_radial20261114'
    cfg=json.loads((D/'config.json').read_text())
    primary=json.loads((D/'primary_radial_result.json').read_text())
    independent=json.loads((D/'independent/radial_result.json').read_text())
    cy=json.loads((D/'independent/cyclic_result.json').read_text())
    full=json.loads((D/'full_replay_result.json').read_text())
    assert primary['pass'] and independent['pass'] and full['pass']
    new=cfg['strict_density_upper_bound'];old=cfg['previous_repository_strict_upper']
    with localcontext() as ctx:
        ctx.prec=40
        improvement=100*(Decimal(old)-Decimal(new))/Decimal(old)
    summary.append({'dimension':n,'strict_density_upper_bound':new,
                    'previous_repository_strict_upper':old,'relative_decrease_percent':str(improvement),
                    'correction_degree':degree,'package':cfg['package'],
                    'full_primary_and_independent_pass':True,'density_interval':independent['density_interval']})
    (D/'README.md').write_text(f'''# Dimension {n}: radial continuation of 14 November 2026

The exact witness proves **density <{new}**, improving the previous
repository bound {old} by **{improvement:.6f}%**. This is an unrestricted
packing-density upper bound. The first table of [Cohn's live sphere-packing
page](https://cohn.mit.edu/sphere-packing/), retrieved on 14 November 2026,
still displays the preceding upper value. The preserved source snapshot and
the distinction from the separate lattice-only table are recorded in
[`research/radial_enrichment_20261114/sources`](../../../research/radial_enrichment_20261114/sources/).
This comparison does not assert an exhaustive worldwide priority search.

Run from a checkout with Python 3.11+ and `python-flint==0.9.0`:

```sh
python {cfg['package']}/verify.py
```

Never use Python `-O`. No solver or network is needed for proof replay.
The entrypoint checks frozen input hashes and the current registry, then runs
all five primary and independent angular, cyclic and radial gates. Actual
completed replay status is in `full_replay_result.json`.

The angular polynomial is unchanged from `{cfg['previous_repository_package']}`.
The new radial function is a nonnegative multiple of the original degree-111
48D seed plus a degree-{degree} correction in normalized dimension-{n}
Laguerre Fourier eigenfunctions. Pair coupling is proved through X=40, and
physical negativity only from X=40 onward. The exact M guard is
1000001/1000000. Read [PROOF.md](PROOF.md) for the global argument and
[independent/AUDIT.md](independent/AUDIT.md) for verification scope.

All preceding witnesses remain available. The discovery history and
deterministic exact-export checks are separate from these proof inputs.
''')
    (D/'PROOF.md').write_text(f'''# Proof of the {n}D upper bound

The conclusion is

    Delta_{n} < {new}.

The underlying theorem is the existing cyclic Cohn–de Laat–Salmon continuum
framework ([arXiv:2206.15373](https://arxiv.org/abs/2206.15373), Section 2.2,
equation (2.4), and Theorem 2.4), together with its periodic PSD summation
argument. The cyclic sum condition is essential; the individual-sign premise
of the printed Theorem 1.1 is not substituted for it. This certificate has no
external neighbor-count input.

## Exact functions and normalization

Let s=6, kappa=17/8, X=pi*|x|^2, and h(x)=exp[-kappa*(X-6*pi)].
`certificate.json` gives rational degree-111 polynomials P and H and a positive
rational M. Set f(x)=P(X)*exp(-X). Both the direct coefficient-sum formula
and an independent Laguerre three-term recurrence prove exactly that H is
the dimension-{n} Fourier polynomial of f. In particular, fhat(0)=H(0)>0.

The angular polynomial Q has degree 19, normalized dimension-{n} Gegenbauer
coefficients b_l>=0, b_0=1 and B=Q(1)=sum(b_l)>0. Primary closed-form and
independent recurrence/moment calculations reconstruct these identities.
For allowed nonzero vectors put

    K(x,y) = (M/B)*h(x)*h(y)*Q(xhat dot yhat).

The addition theorem gives positive semidefinite nonorigin angular features.
In degree zero choose the origin feature -sqrt(M*B) and the nonorigin feature
sqrt(M/B)*h(x). Positive-degree origin features vanish. Thus

    K(0,0)=M*B,  K(0,x)=-M*h(x),  K(x,x)=M*h(x)^2.

The exact identity (M*B)*(M/B)=M^2 proves the rank-one Gram condition. A smooth
cutoff chi, zero near the origin and one for |x|>=sqrt(6), extends the features
globally: use sqrt(M/B)*chi*h-sqrt(M*B)*(1-chi) in degree zero and chi*h times
the positive-degree angular features. This preserves every required value
and makes the origin extension smooth. Features decay rapidly at infinity.

## Complete cyclic condition

Both angular implementations prove Q<-1/100000 on [-1,1/2], Q<0 on
[-1,5001/10000], and precisely one subsequent root before 1, using exact
rational Sturm chains. Nonnegative harmonic coefficients imply |Q|<=B.

For sorted normalized triangle coordinates 0<=t<=1/2, 1<=u<=5, 2*t*u<=1,
define w=sqrt(1+u^2-2*t*u), z1=(1-u*t)/w, z2=(u-t)/w. Complete interval
subdivision proves

    Q(t)+exp[-40*(1-2*t*u)]*Q(z1)
        +exp[-40*(u^2-2*t*u)]*Q(z2) < 0.

The primary check uses the inherited natural interval radicand. The independent
check uses exact monotone rational radicand bounds and nonnegative sorted-gap
bounds. It certifies {cy['certified_rectangles']} rectangles and discards only
rectangles with 2*t*u>1. Subdivision is a full-domain proof.

The single-sign-transition scale lemma extends the result to every larger
shortest side and every alpha>=40. The actual parameter is
alpha=6*pi*kappa=51*pi/4>40. To see the lemma, write the secondary terms as
A*exp(-p*T)+C*exp(-q*T), 0<=p<=q and T>=1. The ordered angles give either two
nonnegative terms, two nonpositive terms, or A<=0<=C. In the mixed case,
factoring exp[-p*(T-1)] shows the sum is at most the maximum of its value at
T=1 and zero. Since Q(t)<0, each case preserves strict negativity.

For t<=0 the secondary exponential gaps are at least 1. For u>=5 a positive
first secondary contribution forces its gap above 4.001; the second gap is
at least 24. Exact arithmetic checks 40^30/30!>2*B*100000, bounding these
positive tails below the strict negative cap margin. These analytic arguments
cover all scales and the infinite triangle tail; the compact rectangle proof
alone is not the global theorem.

## Global radial and pair gates

Complete certified complex-root isolation accounts for all 111 roots of each
polynomial and excludes real roots of H on [0,infinity) and of P on
[40,infinity). Exact H(0)>0 and P(40)<0 establish the ray signs. Both
implementations use FLINT/Arb root isolation; they do not constitute separate
root-isolation backends.

The required pair inequality is

    P(X)*exp(-X) <= M*(2*h(X)-h(X)^2),  X>=6*pi.

It is equivalent, after division by positive h, to

    P(X)*exp[(9/8)*X-(51/4)*pi]
        <= M*(2-exp[-(17/8)*(X-6*pi)]).

The primary proof uses centered polynomial Taylor upper bounds and monotone
exponentials on an enclosing rational interval [6*pi_lower,40]. It accepts
{primary['pair_intervals']} intervals, maximum depth {primary['pair_depth']}.
The independent proof encloses the complete signed expression and explicitly
checks that its {independent['pair_intervals']} accepted leaves tile that
interval, maximum depth {independent['pair_depth']}. Machin bounds obtained
from alternating arctangent sums enclose pi. Beyond 40, P<0 and 0<h<=1 make
the pair inequality automatic.

`tail_escape.json` gives an exactly positive P value at a rational point
beyond 7*pi. Thus the earlier stronger physical-negativity-from-7*pi template
cannot certify this witness; the extended pair proof is necessary. The
positive value is compatible with the actual coupled inequality.

## Density conclusion and strict rounding

The cyclic periodic summation and PSD hypotheses bound center density by
(P(0)+M*B)/H(0). Absolute convergence follows from rapid decay. Periodic
approximation extends the bound to unrestricted packings. Multiplying by
the radius-sqrt(6)/2 ball volume gives

    Delta_{n} <= (3*pi/2)^({n}/2)/Gamma({n}/2+1)
                 * (P(0)+M*B)/H(0).

For even n the upper Machin bound for pi directly gives rational strict
rounding. For odd n=2m+1 the squared volume is
4*pi^(n-1)*6^n*((m+1)!)^2/((n+1)!)^2, so rational comparison of squared
positive quantities proves the same conclusion without relying on displayed
floating digits. The independent density enclosure is

    {independent['density_interval']}

Both strict-rounding checks prove the stated decimal upper bound. Numerical
LP objectives, dense samples and the search checkpoint's solver status are
not proof inputs. This is a computer-assisted certificate, without a claim
of optimality, formal proof-assistant verification or external peer review.
''')
    (D/'independent/AUDIT.md').write_text(f'''# Independent implementation audit

Fresh `audit_angular_cyclic.py` and `audit_radial.py` replays both pass for
dimension {n}. `angular_result.json`, `cyclic_result.json` and
`radial_result.json` record their actual results. The frozen witness data are
copied under `inputs/` and hash-linked by the parent entrypoint.

The angular audit reconstructs the target-dimensional Gegenbauer basis by a
three-term recurrence, checks its constant coefficient by spherical moments,
uses rational Sturm chains, and proves the complete compact cyclic domain
with changed radicand enclosures. The scale and infinite-tail arguments are
given in [the proof](../PROOF.md).

The radial audit reconstructs the Fourier transform by recurrence, independently
of the primary coefficient-sum formula. It checks the full signed pair gap,
the exact union of all accepted intervals, the origin Gram identity, positive
normalization, and rational strict rounding. It shares FLINT/Arb root
isolation with the primary implementation. Independence here concerns code
and enclosures; it does not mean an independent arithmetic library or an
external human audit. No optimization solver is used.
''')
(HERE/'NEW_BOUNDS.json').write_text(json.dumps(summary,indent=2)+'\n')
table='| Dimension | Previous repository strict upper | New strict upper | Decrease |\n| --- | --- | --- | --- |\n'
for row in summary:
    table+=f"| {row['dimension']} | {row['previous_repository_strict_upper']} | {row['strict_density_upper_bound']} | {Decimal(row['relative_decrease_percent']):.6f}% |\n"
(HERE/'README.md').write_text('''# Radial enrichment: 14 November 2026

Three newly computed exact certificates improve the current unrestricted
sphere-packing upper bounds. The complete primary and independent proof
chains pass for every new witness.

'''+table+'''
The old angular kernels remain exactly unchanged. The successful changes are
degree-47 corrections with extended pair coupling in 45D and 47D, and a
degree-63 correction with the same extended coupling in 46D. Minimum squared
separation is 6, kappa=17/8, coupling covers X in [6*pi,40], and physical
negativity is imposed only from 40 onward. Every frozen M receives the exact
guard 1000001/1000000. Exact positive points beyond 7*pi show why the previous
physical-tail sign template does not cover these witnesses.

## Reproduce proofs

Install the root verification requirements and run each new package:

```sh
python verifications/dim45/refinement_radial20261114/verify.py
python verifications/dim46/refinement_radial20261114/verify.py
python verifications/dim47/refinement_radial20261114/verify.py
python research/radial_enrichment_20261114/verify_exports.py
python research/radial_enrichment_20261114/test_registration.py
```

`python run_all.py --jobs 4` includes every new proof and all preserved upper
certificates. Assertions must remain enabled. The new proof packages depend
only on Python and python-flint 0.9.0.

## Bounded discovery and failures

`search.py` is a finite sampled LP, not a verifier. Its span is a nonnegative
multiple of the original normalized degree-111 48D source plus normalized
target-dimensional Laguerre Fourier eigenfunctions. `model.py` constructs
the exact source and basis. Degrees 47, 63 and 47 were selected for 45D, 46D
and 47D respectively. No wider dimension sweep or hidden optimizer run is
claimed. `continue_search.py` adds exact-root-guided points and numerical pair
separators, retaining checkpoints and rejection logs under each `history/`.

Initial sampled witnesses failed on narrow Fourier/physical sign islands and
pair peaks between grid nodes. Root cuts repaired those failures. One 45D
HiGHS attempt returned an unknown/infeasible status; retry with the alternate
LP method recovered the feasible continuation. The first rational exports
used a 100ppm M guard; the final frozen exports use 1ppm and pass complete
interval verification. Earlier checkpoints remain research, even when their
numerical objectives are lower than the certified final objective.

For optional fresh discovery install NumPy, SciPy and mpmath alongside
python-flint; tested versions are recorded in `environment.json`. Example:

```sh
python research/radial_enrichment_20261114/search.py 46 --degree 63
```

`export_check.py` performs rationalization and radial proof gates. Its output
alone is not a density certificate: complete angular/cyclic replay is also
required. `verify_exports.py` deterministically reconstructs the frozen final
certificates and proves the new radial span degrees exactly. The proof
packages read frozen rational inputs rather than trusting discovery floats.

## Comparisons and replay scope

The source commit is 72028bb0c5eaaaeede50cd418ca09339ee404285. The first
(unrestricted-packing) table at https://cohn.mit.edu/sphere-packing/ was
successfully fetched on 14 November 2026. Its 45D, 46D and 47D upper entries
match the preceding repository results to the displayed precision. Raw HTML
and row extracts are preserved under `sources/`. The separate lattice-only
table is not a general-packing comparator. External comparison witnesses are
not independently replayed, and this is not an exhaustive worldwide priority
search.

Historical proof inputs remain unchanged. Registration wrappers accept the
newer package while continuing to verify their own frozen bounds. The full
replay also exposed a pre-existing 16D historical registration assertion that
still required the superseded package as current; its compatibility check was
updated without changing its witness, hashes, mathematical gates or bound.
Completed full-suite status and the initial failed attempt are recorded in
`REPLAY_STATUS.json` and the retained logs, rather than conflated with the
successful new-certificate replays.
''')
print('WROTE proof notes and comparison summary')
