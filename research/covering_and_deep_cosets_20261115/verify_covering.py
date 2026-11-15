"""Strengthen greedy torus offsets with the classical covering simplex bound.

Published input: Coxeter--Few--Rogers (1959). COVERING_PROOF.md derives
the finite bound C_n=(n+1)*(n/(n+3))^(n/2), not an asymptotic estimate.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import factorial
import json
import argparse
import sys

if not __debug__:
    raise RuntimeError('Run without -O.')
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'range48_96_20261115'))
from verify_range import parent_inputs, torus_scan, determinant, decsqrt, pi_interval


def gain_interval(ratio2):
    with localcontext() as ctx:
        ctx.prec=60
        r=(Decimal(ratio2.numerator)/Decimal(ratio2.denominator)).sqrt()
        f=int(((r-1)*100)*10**6)
    lo,hi=F(f,10**6),F(f+1,10**6)
    assert (1+lo/100)**2 <= ratio2 < (1+hi/100)**2
    return [str(lo),str(hi)]


def audit():
    parents=parent_inputs()
    old=json.loads((P/'previous_bounds.json').read_text())
    previous={a['dimension']:a for a in old['rows']}
    assert set(previous)==set(range(48,97))
    candidates=torus_scan(parents)
    _,pi_hi=pi_interval()
    out=[]
    for a in candidates:
        n=a['parent_dimension'];r=a['dimension']-n
        assert n%2==0 and 48 <= n < a['dimension'] <=96
        C=F(n+1)*F(n,n+3)**(n//2)
        assert C>1
        moment=F(n*2*factorial(n-1),factorial(n+1))
        assert moment==F(2,n+1)
        radial=F(2**(n//2-1)*factorial(n//2-1),factorial(n-1))
        normalized=(F((n+1)*factorial(n),factorial(n//2))
                    *F(n,2*(n+1))**(n//2)*radial/(1+moment)**(n//2))
        assert normalized==C
        parent=F(a['parent_density_squared']);N=a['colors']
        Vhi=pi_hi**(n//2)/factorial(n//2)
        def fits(d):
            return (N-1)**2*parent*Vhi**2*d**n < C*C
        lo,hi=0,4_000_000
        assert not fits(F(hi,10**6))
        while lo+1 < hi:
            mid=(lo+hi)//2
            if fits(F(mid,10**6)):lo=mid
            else:hi=mid
        d=F(lo,10**6)
        assert fits(d) and not fits(d+F(1,10**6))
        lam=F(a['same_color_ratio']);t=max(4-d,4/lam)
        assert d+t>=4 and lam*t>=4
        d2=parent*(2/t)**r/determinant(a['vertical_root_gram'])
        assert d2>=F(a['density_squared'])
        out.append({'dimension':a['dimension'],'name':a['name'].replace('torus','CFR torus'),
                    'density_squared':str(d2),'center_density':decsqrt(d2),
                    'parent_dimension':n,'parent_density_squared':str(parent),'colors':N,
                    'vertical_root_gram':a['vertical_root_gram'],
                    'same_color_ratio':str(lam),'torus_squared_separation':str(d),
                    'vertical_squared_minimum':str(t),'covering_density_lower_bound':str(C),
                    'strict_covering_inequality_exact':True,'explicit_offsets':False})
    best={}
    for a in out:
        n=a['dimension']
        if n not in best or F(a['density_squared'])>F(best[n]['density_squared']):best[n]=a
    improved=[]
    for n,a in sorted(best.items()):
        ratio2=F(a['density_squared'])/F(previous[n]['density_squared'])
        if ratio2>1:
            improved.append(dict(a,previous=previous[n],ratio_to_previous=decsqrt(ratio2),
                                 gain_percent_interval=gain_interval(ratio2)))
    c96=next(a for a in improved if a['dimension']==96)
    bcs=F(3**96,2**48);ratio2=F(c96['density_squared'])/bcs
    return {'published_theorem':'Coxeter--Few--Rogers, Covering space with equal spheres (1959), DOI 10.1112/S0025579300002059',
            'finite_bound_derived':'C_n=(n+1)*(n/(n+3))^(n/2), from simplex solid angles and Jensen',
            'maximum_dimension':96,'candidate_count':len(out),'candidates':out,
            'improved_dimensions':improved,
            'BCS96_comparison':{'reported_construction_center_density_squared':str(bcs),
                               'new_to_reported_ratio_squared':str(ratio2),
                               'new_to_reported_ratio':decsqrt(ratio2),
                               'gain_percent_interval':gain_interval(ratio2)},
            'all_arithmetic_inequalities_exact':True,'worldwide_priority_proved':False}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args()
    result=audit();path=P/'covering_results.json'
    if args.write_result:path.write_text(json.dumps(result,indent=2)+'\n')
    else:assert json.loads(path.read_text())==result,'Stored covering result is stale.'
    for a in result['improved_dimensions']:
        print(a['dimension'],a['center_density'],'gain over checkpoint',a['gain_percent_interval'])
    print('96D ratio to BCS1982:',result['BCS96_comparison']['new_to_reported_ratio'])
