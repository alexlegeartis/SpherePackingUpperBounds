"""Sharper derived section guarantees using a finite Rogers simplex estimate.

See SECTIONS_PROOF.md. These are derived classical controls, not numbers
reported as dimension-specific records by the parent-lattice authors.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,prod,isqrt
from decimal import Decimal,localcontext,ROUND_CEILING
from functools import lru_cache
import json,argparse,sys

if not __debug__:raise RuntimeError('Run without -O.')
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'range48_96_20261115'))
from verify_range import parent_inputs,pi_interval,section_floor,decsqrt


def sqrt_upper(x):
    scale=10**22;a=isqrt(x.numerator*scale*scale//x.denominator)
    if a*a*x.denominator<x.numerator*scale*scale:a+=1
    out=F(a,scale)
    assert out*out>=x and (out-F(1,scale))**2<x
    return out


def half_power_upper(x,n):
    return x**(n//2)*(sqrt_upper(x)if n%2 else 1)


@lru_cache(None)
def hermite(n):
    lo,hi=pi_interval()
    # E[S]=2/(n+1), Var(S)=4(n-1)/((n+1)^2(n+2)(n+3)).
    # Taylor's theorem with f'' decreasing gives the explicit upper mean.
    correction=F(n*(n-1),2*(n+1)**2*(n+3))
    second=F(24*n+4*n*(n-1),n*(n+1)*(n+2)*(n+3))
    variance=second-F(2,n+1)**2
    assert variance==F(4*(n-1),(n+1)**2*(n+2)*(n+3))
    assert correction==F(n*(n+2),8)*variance
    average=(half_power_upper(F(n+1,n+3),n)
             +correction*half_power_upper(F(n,n+1),n)*F(n,n+1)**2)
    upper=F(n+1)*half_power_upper(F(1,2),n)*average
    volume=(lo**(n//2)/factorial(n//2)if n%2==0 else
            F(2**(n//2+1),prod(range(1,n+1,2)))*lo**(n//2))
    bound=4**n*(upper/volume)**2
    with localcontext()as ctx:
        ctx.prec=80
        x=Decimal(bound.numerator)/Decimal(bound.denominator)
        a=int(((x.ln()/n).exp()*10**9).to_integral_value(rounding=ROUND_CEILING))
    c=F(a,10**9)
    assert c**n>=bound
    # Establish improvement over the UNROUNDED Blichfeldt expression,
    # rather than only over the preceding centesimal numerical enclosure.
    if n%2==0:blich_lower=(2/hi)**n*factorial(n//2+1)**2
    else:
        j=n//2+2;g=F(factorial(2*j),4**j*factorial(j))
        blich_lower=F(2)**n*g*g/hi**(n-1)
    assert c**n<blich_lower
    return c,upper


def audit():
    previous={a['dimension']:a for a in json.loads((P/'previous_bounds.json').read_text())['rows']}
    candidates=[]
    for p in parent_inputs():
        d2=F(p['density_squared'])
        for n in range(p['dimension'],48,-1):
            c,upper=hermite(n)
            d2=section_floor(d2,n,c)
            if not d2:break
            candidates.append({'dimension':n-1,'name':p['name']+f' / Rogers section to {n-1}D',
                               'density_squared':str(d2),'center_density':decsqrt(d2),
                               'parent_dimension':p['dimension'],'parent_source':p['evidence'],
                               'derived_not_reported_dimension_specific_record':True})
    best={}
    for a in candidates:
        n=a['dimension']
        if n not in best or F(a['density_squared'])>F(best[n]['density_squared']):best[n]=a
    improved=[a for n,a in sorted(best.items())if F(a['density_squared'])>F(previous[n]['density_squared'])]
    return {'method':'Rogers packing simplex bound with a finite Dirichlet second-moment upper estimate',
            'maximum_dimension':96,'hermite_upper':{str(n):str(hermite(n)[0])for n in range(49,97)},
            'all_gamma_bounds_stronger_than_unrounded_Blichfeldt':True,
            'improved_controls':improved,'all_density_and_root_inequalities_exact':True,
            'worldwide_record_claim':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args()
    result=audit();path=P/'sections_results.json'
    if args.write_result:path.write_text(json.dumps(result,indent=2)+'\n')
    else:assert json.loads(path.read_text())==result
    for a in result['improved_controls']:print(a['dimension'],a['center_density'],a['name'])
