"""Primary verifier helpers. Exact rational data, with outward Arb enclosures."""
import json,math
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_poly as P,arb

def rational(z):return Q(str(Fraction(str(z))))
def normalized_gegenbauer(n,d):
    x=P([0,1]);g=[P([1])]
    if d:g.append(x)
    for k in range(1,d):g.append(((2*k+n-2)*x*g[-1]-k*g[-2])/(k+n-2))
    return g
def pi_bounds(terms=128):
    def atan(x):
        lo=sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(terms)),Q(0))
        return lo,lo+x**(2*terms+1)/(2*terms+1)
    al,au=atan(Q(1,5));bl,bu=atan(Q(1,239))
    return 16*al-4*bu,16*au-4*bl
def inputs(directory):
    D=Path(directory).resolve()
    r=json.loads((D/'certificate.json').read_text());cfg=json.loads((D/'config.json').read_text())
    angular=json.loads((D/'angular.json').read_text())
    n=r['dimension'];assert n==cfg['dimension']==angular['dimension'] and n>=3
    s,a,b,q=map(rational,[r['s_exact'],r['a_exact'],r['b_exact'],r['q_exact']])
    assert s>0 and a+b>0 and b>0 and q==1
    assert a==rational(angular['profile_a_exact']) and b==rational(angular['profile_b_exact'])
    assert rational(r['strict_density_upper_bound'])==rational(cfg['strict_density_upper_bound'])
    return D,r,cfg,angular,n,s,a,b,q
def rational_rounding(objective,n,s,pi_upper,strict):
    assert objective>0 and strict>0
    if n%2==0:assert objective*(pi_upper*s/4)**(n//2)/math.factorial(n//2)<strict
    else:
        k=(n-1)//2
        assert objective**2*4*pi_upper**(n-1)*s**n*Q(math.factorial(k+1),math.factorial(n+1))**2<strict**2
