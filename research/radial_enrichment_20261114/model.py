"""Exact radial construction. Numerical discovery is kept separate from proof."""
from pathlib import Path
import json
from flint import fmpq as Q, fmpq_poly as P, arb, ctx

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'verifications/dim48/continuum_three_point/adaptive_cyclic40_radial111.json'
PACKAGES = {45: 'dim45', 46: 'dim46/refinement_extended_pair',
            47: 'dim47/refinement_degree47'}

def angular_path(n):
    return ROOT / 'verifications' / PACKAGES[n] / 'continuum_three_point/cyclic_angular40_exact.json'

def basis(n, degree, scale=2):
    a = Q(n, 2)
    x = P([0, scale])
    ls = [P([1])]
    if degree:
        ls.append(a-x)
    for k in range(1, degree):
        ls.append(((2*k+a-x)*ls[-1]-(k+a-1)*ls[-2])/(k+1))
    return ls

def fourier(poly, n):
    ls = basis(n, poly.degree(), 1)
    result, factorial = P([]), 1
    for k in range(poly.degree()+1):
        if k:
            factorial *= k
        result += poly[k]*factorial*ls[k]
    return result

def inputs(n, degree):
    raw = json.loads(SOURCE.read_text())
    p = P(list(map(Q, raw['physical_polynomial'])))/Q(raw['fourier_polynomial'][0])
    h = fourier(p, n)
    ls = [l/l(0) for l in basis(n, degree)]
    angular = json.loads(angular_path(n).read_text())
    assert angular['dimension'] == n
    assert Q(angular['gegenbauer_coefficients'][0]) == 1
    return p, h, ls, Q(angular['B_exact'])

def pi_bounds():
    def atan(x):
        lo = sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(128)), Q(0))
        return lo, lo+x**257/257
    a,b = atan(Q(1,5))
    c,d = atan(Q(1,239))
    return 16*a-4*d, 16*b-4*c

def density(objective, n):
    return arb(objective)*(arb.pi()*arb(3)/2)**(arb(n)/2)/(arb(n)/2+1).gamma()
