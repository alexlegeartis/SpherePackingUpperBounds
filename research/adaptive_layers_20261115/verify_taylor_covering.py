"""Finite cubic-Taylor improvement to the covering-simplex lower estimate."""
from fractions import Fraction as F
from math import factorial,prod
from decimal import Decimal,localcontext
import json,argparse,sys
from common import P,determinant,finish,checkpoint
if not __debug__:raise RuntimeError('Run without -O.')
sys.path.insert(0,str(P.parent/'covering_and_deep_cosets_20261115'))
from verify_covering import audit as preceding_audit
from verify_to96 import pi_interval


def correction(n):
    mu=F(2,n+1)
    m2=F(24*n+4*n*(n-1),prod(range(n,n+4)))
    m3=F(720*n+144*n*(n-1)+8*n*(n-1)*(n-2),prod(range(n,n+6)))
    var=m2-mu**2;third=m3-3*mu*m2+2*mu**3
    assert var==F(4*(n-1),(n+1)**2*(n+2)*(n+3))and third>=0
    a=F(n,2)
    factor=1+a*(a+1)*var/(2*(1+mu)**2)-a*(a+1)*(a+2)*third/(6*(1+mu)**3)
    assert factor>1
    return factor,{'mean':str(mu),'second_moment':str(m2),'third_moment':str(m3),
                  'variance':str(var),'third_central_moment':str(third),'correction_factor':str(factor)}


def audit():
    old=preceding_audit();_,pi_hi=pi_interval();out=[];moments={};scale=10**9
    for a in old['candidates']:
        n=a['parent_dimension'];r=a['dimension']-n;N=a['colors'];parent=F(a['parent_density_squared'])
        assert n%2==0 and n<a['dimension']<=96
        factor,proof=correction(n);moments[str(n)]=proof
        C=F(n+1)*F(n,n+3)**(n//2)*factor
        Vhi=pi_hi**(n//2)/factorial(n//2)
        def fits(d):return(N-1)**2*parent*Vhi**2*d**n<C*C
        lo,hi=0,4*scale
        while lo+1<hi:
            mid=(lo+hi)//2
            if fits(F(mid,scale)):lo=mid
            else:hi=mid
        d=F(lo,scale);assert fits(d)and not fits(d+F(1,scale))
        lam=F(a['same_color_ratio']);t=max(4-d,4/lam)
        assert d+t>=4 and lam*t>=4
        d2=parent*(2/t)**r/determinant(a['vertical_root_gram'])
        assert d2>=F(a['density_squared'])
        out.append({'dimension':a['dimension'],'parent_dimension':n,'colors':N,'density_squared':str(d2),
                    'covering_lower_bound':str(C),'torus_squared_depth':str(d),'vertical_squared_minimum':str(t),
                    'vertical_root_gram':a['vertical_root_gram'],'same_color_ratio':str(lam),
                    'strict_covering_inequality_exact':True})
    best={}
    for a in out:
        n=a['dimension']
        if n not in best or F(a['density_squared'])>F(best[n]['density_squared']):best[n]=a
    before=checkpoint();rows=[]
    for n,a in sorted(best.items()):
        d2=F(a['density_squared'])
        if d2>F(before[n]['density_squared']):
            rows.append(finish(n,'Cubic covering-simplex bound and colored layers',d2,'new_periodic_existence_proof',a))
    return {'method':'Coxeter-Few-Rogers simplex bound; global cubic Taylor lower polynomial for (1+s)^(-n/2)',
            'moments_and_corrections':moments,'candidate_count':len(out),'rows':rows,
            'all_covering_and_density_inequalities_exact':True,'maximum_dimension':96,
            'explicit_offsets':False,'inherited_inputs':'Coxeter-Few-Rogers1959; named parent-lattice and color-period theorems from the preceding audited package'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args();r=audit()
    if args.write_result:(P/'covering_results.json').write_text(json.dumps(r,indent=2)+'\n')
    else:assert json.loads((P/'covering_results.json').read_text())==r
    for a in r['rows']:print(a['dimension'],a['center_density'],'gain',a['gain_over_checkpoint_percent'])
