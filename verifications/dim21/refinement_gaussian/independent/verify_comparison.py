"""Exact-source and rigorous arithmetic comparison for the 21D certificate.

--write regenerates comparison.json, then verifies it at higher precision.
No claim of lower-construction certification or exhaustive worldwide priority.
"""
from pathlib import Path
from fractions import Fraction
import json, math, sys
from flint import arb, fmpq, ctx
if not __debug__:
    raise RuntimeError('Assertions are required')
root = Path(__file__).resolve().parent
source = 'https://cohn.mit.edu/sphere-packing/'
published, previous, new = '0.0051596603948176', '0.005149549312', '0.004985977000'
formula = '2^-19*pi^10*11!*43^(21/2)/(22!*3^(23/2))'
def exact(s):
    f = Fraction(s)
    return arb(fmpq(f.numerator, f.denominator))
def calculate():
    # Indexed table center density: 2^-41 * 3^(-23/2) * 43^(21/2).
    center = arb(43)**10 * arb(43).sqrt() / (arb(2)**41 * arb(3)**11 * arb(3).sqrt())
    lower = center * arb.pi()**(arb(21)/2) / (arb(23)/2).gamma()
    # Independently use the half-integer gamma closed form.
    closed = arb.pi()**10 * math.factorial(11) * arb(43)**10 * arb(43).sqrt() / (arb(2)**19 * math.factorial(22) * arb(3)**11 * arb(3).sqrt())
    assert lower.overlaps(closed)
    old, prev, fresh = map(exact, [published, previous, new])
    assert 0 < lower < fresh < prev < old
    return {'lower_density': lower,
            'upper_decrease_percent': 100*(old-fresh)/old,
            'gap_decrease_percent': 100*(old-fresh)/(old-lower),
            'previous_repository_upper_decrease_percent': 100*(prev-fresh)/prev,
            'previous_repository_gap_decrease_percent': 100*(prev-fresh)/(prev-lower)}
if '--write' in sys.argv:
    ctx.prec = 256
    out = {'dimension': 21, 'upper_old': published, 'upper_new': new,
           'previous_repository_upper': previous,
           'lower_formula': formula, 'lower_source': source,
           'center_density_formula': '2^-41*3^(-23/2)*43^(21/2)',
           'status': 'Pinned published and prior repository comparisons; no exhaustive worldwide priority claim'}
    out.update({k: str(v) for k,v in calculate().items()})
    (root/'comparison.json').write_text(json.dumps(out,indent=2)+'\n')
ctx.prec = 512
stored = json.loads((root/'comparison.json').read_text())
assert stored['dimension'] == 21
assert stored['upper_old'] == published and stored['upper_new'] == new
assert stored['previous_repository_upper'] == previous
assert stored['lower_formula'] == formula and stored['lower_source'] == source
for k,v in calculate().items():
    assert arb(stored[k]).contains(v), k
out = {'pass': True, 'dimension': 21, 'precision_bits': 512,
       'lower_center_density_and_closed_form_agree': True,
       'scope': 'Pinned table and prior repository values, not independently reconstructed lower construction or current priority'}
out.update({k: str(v) for k,v in calculate().items()})
(root/'comparison_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
