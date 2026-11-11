"""Recompute sourced density comparisons with outward-rounded Arb arithmetic.

These are checks against pinned table entries, not proof of worldwide priority.
The lower constructions themselves are not reconstructed by this checker.
"""
from pathlib import Path
from fractions import Fraction
import json
from flint import arb, fmpq, ctx

if not __debug__:
    raise RuntimeError('Assertions are required')
ctx.prec = 512
root = Path(__file__).resolve().parent
comparison = json.loads((root / 'comparison.json').read_text())
n = comparison['dimension']
# Independently transcribed center densities from Cohn's indexed table.
centers = {
    29: lambda: 1 / arb(2).sqrt(),
    30: lambda: arb(1),
    31: lambda: arb(3)**15 / (arb(2)**23 * arb(2).sqrt()),
    33: lambda: arb(3)**16 * arb(3).sqrt() / arb(2)**25,
    35: lambda: 2 * arb(2).sqrt(),
    37: lambda: 4 * arb(2).sqrt(),
}
expected = {
    29: ('0.0003598581852089', '0.0003420173821', '2^30*pi^14*15!/(30!*sqrt(2))'),
    30: ('0.0002559028743732', '0.0002454113168', 'pi^15/15!'),
    31: ('0.0001817083813917', '0.0001763894399', '2^8*sqrt(2)*pi^15*3^15*16!/32!'),
    33: ('0.0000912356039023', '0.00009009682341', '2^9*pi^16*3^16*sqrt(3)*17!/34!'),
    35: ('0.0000455743843107', '0.00004463661144', '2^37*sqrt(2)*pi^17*18!/36!'),
    37: ('0.0000226586900106', '0.00002171160843', '2^40*sqrt(2)*pi^18*19!/38!'),
}
old, new, formula = expected[n]
assert Fraction(comparison['upper_old']) == Fraction(old)
assert Fraction(comparison['upper_new']) == Fraction(new)
assert comparison['lower_formula'] == formula
assert comparison['lower_source'] == 'https://cohn.mit.edu/sphere-packing/'
certificate = json.loads((root / 'certificate.json').read_text())
assert certificate['dimension'] == n
half_n = arb(n) / 2
lower = centers[n]() * arb.pi()**half_n / (half_n + 1).gamma()
def exact_decimal(value):
    rational = Fraction(value)
    return arb(fmpq(rational.numerator, rational.denominator))
uold, unew = exact_decimal(old), exact_decimal(new)
assert 0 < lower < unew < uold
upper_decrease = 100 * (uold - unew) / uold
gap_decrease = 100 * (uold - unew) / (uold - lower)
for key, value in [('lower_density', lower), ('upper_decrease_percent', upper_decrease), ('gap_decrease_percent', gap_decrease)]:
    assert arb(comparison[key]).contains(value), key
out = {'pass': True, 'dimension': n, 'precision_bits': ctx.prec,
       'lower_from_center_density_interval': str(lower),
       'upper_decrease_percent_interval': str(upper_decrease),
       'gap_decrease_percent_interval': str(gap_decrease),
       'scope': 'Pinned table entries; lower construction and worldwide priority not certified'}
(root / 'comparison_audit.json').write_text(json.dumps(out, indent=2) + '\n')
print(out)
