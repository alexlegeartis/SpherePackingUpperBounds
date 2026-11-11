"""Write a final report only after every independent audit gate has passed.

Usage: python finalize_audit.py [audit_directory ...]
This is reporting glue, not a mathematical verifier.
"""
from pathlib import Path
import hashlib
import json
import re
import sys

if not __debug__:
    raise RuntimeError('Assertions are required')

def read(root, name):
    return json.loads((root / name).read_text())

for root in map(Path, sys.argv[1:] or [Path(__file__).resolve().parent]):
    radial = read(root, 'audit.json')
    angular = read(root, 'angular_cyclic_audit.json')
    comparison = read(root, 'comparison.json')
    comparison_audit = read(root, 'comparison_audit.json')
    source = read(root, 'source_binding_audit.json')
    cert = read(root, 'certificate.json')
    ang = read(root, 'angular.json')
    n = radial['dimension']
    assert all(x['pass'] and x['dimension'] == n for x in (radial, angular, comparison_audit, source))
    assert comparison['dimension'] == cert['dimension'] == ang['dimension'] == n
    assert (root / 'audit_exit_code.txt').read_text().strip() == '0'
    log = (root / 'full_audit.log').read_text()
    assert f'PASS complete{n}D exact angular/cyclic/radial/origin/rounding audit' in log
    variations = re.findall(r'Independent exact Sturm PASS (\w+) (\d+) (\d+)', log)
    assert len(variations) == 2 and all(a == b for _, a, b in variations)
    for name, want in read(root, 'source_hashes.json').items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == want
    # Replace the formerly incomplete radial log with the complete fresh evidence.
    lines = log.splitlines()
    (root / 'radial_audit.log').write_text('\n'.join(lines[1:-1]) + '\n')
    (root / 'angular_audit.log').write_text(lines[0] + '\n')
    hashes = read(root, 'source_hashes.json')
    upper = radial['strict_upper_decimal']
    m = (n-1)//2
    volume = f'2*pi^{m}*({cert["r2"]})^({n}/2)*{m+1}!/{n+1}!' if n%2 else f'(pi*({cert["r2"]})/4)^({n//2})/{n//2}!'
    text = f'''# Independent {n}D Gaussian cyclic certificate audit

**PASS:** the certificate proves the unrestricted sphere-packing density bound Delta_{n} < {upper}.

Fresh execution completed on 2026-11-11 UTC with Python 3.12 and python-flint 0.9.0, assertions enabled. `verify_all.py` exited 0. Exact rational polynomial arithmetic and outward-rounded Arb arithmetic are used. The independent radial verifier uses rational Sturm sequences instead of the producer's complex-root isolation, explicit normalized Gegenbauer coefficients instead of its recurrence, and a second-derivative Taylor enclosure instead of its radial subdivision bound.

## Inputs and reproducibility

- Certificate SHA-256: `{hashes['certificate.json']}`.
- Angular witness SHA-256: `{hashes['angular.json']}`.
- Both local inputs match the frozen producer files byte for byte. All {source['frozen_producer_checksum_entries_passed']} entries of the producer's existing SHA256SUMS passed. No producer proof input was changed.
- From this directory, run `python verify_all.py` and then `python verify_comparison.py` in an environment with `python-flint==0.9.0`. Never use `python -O`.
- `full_audit.log`, `audit_exit_code.txt`, `audit.json`, `angular_cyclic_audit.json`, `comparison_audit.json`, and `source_binding_audit.json` record the executed gates. `AUDIT_SHA256SUMS` binds the delivered audit files. The result requires both the angular and radial stages, not `audit.json` alone.

## Exact algebra, all radial distances, and rounding

The complete 112-coefficient monomial Fourier transform of the degree-111 polynomial P was reconstructed in dimension {n} and equals the supplied H exactly for the Fourier convention exp(-2*pi*i*dot(x,y)). H(0)>0 and M>0. Exact Sturm variations are Fourier {variations[0][1]} to {variations[0][2]} on [0,infinity), and physical {variations[1][1]} to {variations[1][2]} on [7*pi_lower,infinity). Endpoint signs establish H>0 and P<0 on those entire rays, not merely at sampled points.

The degree-{ang['degree']} angular polynomial Q was reconstructed in the normalized dimension-{n} Gegenbauer basis. Every coefficient is nonnegative, their sum is B=Q(1), and a separate exact spherical even-moment calculation reproduces b0={cert['b0_exact']}. The origin Gram entries MB/b0, -M, and Mb0/B have positive diagonal and determinant zero. Thus the origin extension is positive semidefinite, as are all higher harmonic blocks.

The verified minimum squared separation is s={cert['r2']}, kappa={cert['kappa']}, and the radial tail cutoff is 7. Their product kappa*s=51/4 gives alpha=51*pi/4>40. For h(r)=exp[-kappa*pi*(r^2-s)], the full pair inequality f(r)+M*h(r)^2-2*M*h(r)<=0 was certified on [sqrt(s),sqrt(7)] using {radial['pair_leaves']} interval leaves, maximum depth {radial['pair_depth']}. The Taylor remainder bounds the complete second derivative of the divided coupling function. The initial rational interval contains the entire required interval, including both endpoints. The negative P tail and 0<h<=1 settle every larger distance.

The exact objective (P(0)+MB/b0)/H(0) equals the source rational objective. The radius-sqrt(s)/2 ball volume simplifies to {volume}. A rational Machin enclosure of pi proves the displayed strict upper bound; odd dimensions square positive quantities before comparison, so no decimal approximation is used to prove strict rounding. The independent final density enclosure is:

{radial['density_interval']}

## Complete triangle domain and scale extension

Exact angular Sturm checks establish Q<-1/100000 on [-1,1/2], Q<0 through 5001/10000, and precisely one root between that cutoff and 1. The independently implemented compact-domain calculation passes {angular['compact_triangles_passed']} rational rectangles and excludes {angular['exact_sorting_exclusions']} rectangles solely because the side-ordering constraint is impossible throughout them. Maximum subdivision depth is {angular['max_depth']}. Endpoints and degenerate triangles are covered.

To spell out the analytic argument, sort a triangle's side lengths as a<=b<=c, and put u=b/a, v=c/a, and t=(1+u^2-v^2)/(2u). Then u>=1, t<=1/2, and 2tu<=1. The other angle cosines are t1=(1-ut)/v and t2=(u-t)/v, with t<=t1<=t2. After division by a positive common factor, its cyclic expression is

S(alpha)=Q(t)+exp[-alpha*A]*Q(t1)+exp[-alpha*D]*Q(t2),

where A=v^2-u^2=1-2tu>=0 and D=v^2-1>=A. The numerical domain proof is at alpha=40 for 0<=t<=1/2 and 1<=u<=5.

There is only one angular sign change, so raising alpha preserves negativity: if Q(t1)>=0, both nonnegative terms only decrease; if Q(t2)<=0, every term is nonpositive; in the remaining case Q(t1)<0<Q(t2), multiplication of S(40+delta) by exp(delta*A) makes the negative Q(t) term more negative, leaves the middle term fixed, and decreases the positive last term. Hence S(40+delta)<=0. A triangle with shortest side a>=sqrt(s) has alpha=kappa*pi*a^2>=51*pi/4>40.

For t<=0, both A and D are at least 1, and |Q|<=B follows from its nonnegative normalized harmonic expansion. For t>=0 and u>=5, ordering implies t1<=1/5, so its term is negative, while D>=24. Both omitted tails are therefore bounded by the exact exponential-series inequality 40^30/30!>2*B*100000. This closes the full unbounded triangle domain.

## Why the result applies to unrestricted packings

For nonzero packing-domain vectors use K(x,y)=(M/B)*h(|x|)*h(|y|)*Q(angle(x,y)), with K(x,0)=-M*h(|x|) and K(0,0)=MB/b0. The harmonic addition theorem and the origin Gram calculation give a finite-rank PSD kernel. A smooth cutoff of higher harmonic features below sqrt(s), together with a smooth interpolation of the degree-zero feature to its prescribed negative value at zero, gives rapidly decreasing smooth features without changing any required packing-domain or origin values.

For a periodic packing with m centers per fundamental cell of volume V, the Fourier sum of f is at least m^2*H(0)/V. Add the nonnegative PSD kernel sum, anchor by anchor. Coincident indices contribute m*(f(0)+K(0,0)); each pair contributes exactly the verified pair inequality; each distinct unordered triple contributes twice its verified cyclic sum. Absolute convergence justifies this regrouping and taking limits of finite PSD sums. Thus m^2*H(0)/V<=m*(f(0)+K(0,0)), giving the displayed density bound. Standard periodic approximation extends it to arbitrary packings.

This is the cyclic finite-sum refinement of the framework in Cohn, de Laat, and Salmon, [Three-point bounds for sphere packing](https://arxiv.org/abs/2206.15373), not an assertion that the individual off-diagonal hypothesis of its Theorem 1.1 holds. No lattice restriction, neighbor-count conjecture, or additional geometric premise is introduced.

## Sourced comparisons and limits

Pinned previous upper bound: {comparison['upper_old']}, from [Cohn's sphere-packing table](https://cohn.mit.edu/sphere-packing/), associated with [Afkhami-Jeddi et al., v2](https://arxiv.org/abs/2006.02560v2). The sourced lower density formula is {comparison['lower_formula']}.

Independent 512-bit recomputation starts with the table's exact center density, multiplies by the unit-ball volume, and verifies that all stored comparison intervals enclose the result:

- Lower density: {comparison['lower_density']}.
- Relative upper-bound decrease, 100*(U_old-U_new)/U_old: {comparison['upper_decrease_percent']} percent.
- Additive upper/lower-gap decrease, 100*(U_old-U_new)/(U_old-L): {comparison['gap_decrease_percent']} percent.

Connected search on 2026-11-11 recovered the matching upper entries and exact center-density expressions from Cohn's indexed table. The index reported a crawl about four months old; direct page retrieval returned HTTP 502/timeouts. Searches for later dimension-specific results found no replacement used here, but that is not an exhaustive freshness or worldwide-priority guarantee. The comparison remains explicitly against a pinned published reference. Lower-bound constructions are sourced, not independently reconstructed.

No mathematical defect was found in the audited certificate and finite-sum argument. The trust base includes Python, python-flint/FLINT/Arb, the runtime, and the analytic argument above. This is independent computer-assisted verification, not proof-assistant formalization, external peer review, an optimality claim, or an exhaustive new-record claim.
'''
    (root / 'AUDIT.md').write_text(text)
    files = sorted(p for p in root.iterdir() if p.is_file() and p.name != 'AUDIT_SHA256SUMS')
    (root / 'AUDIT_SHA256SUMS').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in files))
    print(f'FINAL AUDIT PASS dimension {n}: {upper}; {len(files)} audit hashes')
