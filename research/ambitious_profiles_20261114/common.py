"""Exact algebra and normalized radial profiles for the larger cyclic family."""
import json,sys
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q, fmpq_poly as P, arb,ctx
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from research.radial_enrichment_20261114.model import basis,fourier,pi_bounds
HERE=Path(__file__).resolve().parent
SOURCE=ROOT/'verifications/dim48/continuum_three_point/adaptive_cyclic40_radial111.json'

def source(n,d):
    r=json.loads(SOURCE.read_text())
    p=P(list(map(Q,r['physical_polynomial'])))/Q(r['fourier_polynomial'][0])
    return p,fourier(p,n),[l/l(0) for l in basis(n,d)]

def base_angular(n):
    if n in [45,46,47]:return ROOT/f'verifications/dim{n}/refinement_radial20261114/continuum_three_point/cyclic_angular40_exact.json'
    if n==48:return ROOT/'verifications/dim48/continuum_three_point/cyclic_angular40_exact.json'
    return ROOT/f'verifications/dim{n}/continuum_three_point/cyclic_angular40_exact.json'

def rational(s):return Q(str(F(str(s))))

def h_arb(v,a,b):
    return v**(-a)*(-b*(v-1)).exp()

def volume_arb(n,s):
    return (arb.pi()*arb(s)/4)**(arb(n)/2)/(arb(n)/2+1).gamma()

def normalized_gegenbauer(n,d):
    x=P([0,1]);g=[P([1])]
    if d:g.append(x)
    for ell in range(1,d):g.append(((2*ell+n-2)*x*g[-1]-ell*g[-2])/(ell+n-2))
    return g
