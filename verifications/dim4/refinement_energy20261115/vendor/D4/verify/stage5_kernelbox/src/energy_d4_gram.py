"""Exact rational enclosure of all fixed V7 D4 Gram/count/spherical energies.

This only evaluates the accepted law on a specified lattice neighborhood. It
performs no coefficient optimization, numerical search or certificate search.
"""
from fractions import Fraction as Q
from itertools import product, combinations
from collections import Counter
from math import isqrt, comb
from pathlib import Path
import argparse
import json
import sys

PROJECT = Path(__file__).resolve().parents[3]  # the D4 package root
SOURCE = PROJECT / 'data/v7'
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(40000)

def add(a, b): return a[0]+b[0], a[1]+b[1]
def mul(a, b):
    v=[x*y for x in a for y in b]
    return min(v), max(v)
def scale(a, c): return mul(a, (c,c))
def power(a, n):
    out=(Q(1),Q(1))
    for _ in range(n): out=mul(out,a)
    return out
def square(a):
    return (Q(0) if a[0]<=0<=a[1] else min(a[0]**2,a[1]**2), max(a[0]**2,a[1]**2))
def isum(xs):
    out=(Q(0),Q(0))
    for x in xs: out=add(out,x)
    return out
def sqrt_interval(n):
    den=2**192; k=isqrt(n*den*den)
    ans=Q(k,den),Q(k if k*k==n*den*den else k+1,den)
    assert ans[0]**2<=n<=ans[1]**2
    return ans

def polynomial(cs,ns):
    pairs=[(i,j) for i in range(len(ns)) for j in range(i,len(ns))]
    assert len(cs)==1+len(ns)+len(pairs)
    return Q(cs[0])+2*sum(Q(cs[i+1])*n for i,n in enumerate(ns))+sum(
        Q(c)*ns[i]*ns[j]*(1 if i==j else 2) for c,(i,j) in zip(cs[1+len(ns):],pairs))

def zonal(t,l):
    prev,cur=Q(1),2*t
    if l==0: return prev
    for _ in range(1,l): prev,cur=cur,2*t*cur-prev
    return cur/(l+1)

def even_zonal(t2,l):
    assert l%2==0
    # U_l(t)=sum_j (-1)^j binom(l-j,j)(2t)^(l-2j).
    return sum((-1)**j*comb(l-j,j)*2**(l-2*j)*t2**(l//2-j) for j in range(l//2+1))/Q(l+1)

def evaluate(dilation=Q(1)):
    assert dilation in (Q(1),Q(501473,500000)), "Only two explicitly assigned configurations."
    kernel=json.loads((SOURCE/'kernel.json').read_text())
    origin=json.loads((SOURCE/'origin.json').read_text())
    assert kernel['dimension']==4 and kernel['local_degree']==6 and kernel['angular_degree']==18
    assert len(kernel['blocks'])==19 and len(origin['u'])==49
    assert kernel['radial_basis_intervals']==[['1','3/2']]*6+[['3/2','7/4']]
    breaks=list(map(Q,kernel['breaks']))
    assert breaks[0]==1 and breaks[-1]==Q(7,4)
    pts=[p for p in product(range(-2,3),repeat=4) if sum(p)%2==0 and sum(x*x for x in p) in (2,4,6)]
    shells={s:[p for p in pts if sum(x*x for x in p)==2*s] for s in (1,2,3)}
    assert [len(shells[s]) for s in shells]==[24,24,96]
    assert min(sum((x-y)**2 for x,y in zip(a,b)) for a,b in combinations([(0,)*4]+pts,2))==2
    assert Q(3)*dilation**2<Q(7,4)**2<4*dilation**2  # all supported D4 shells accounted for
    phi={}
    for s,region in [(1,0),(2,5),(3,6)]:
        left,right=map(Q,kernel['radial_basis_intervals'][region]); width=right-left
        r=scale(sqrt_interval(s),dilation); t=scale(add(r,(-left,-left)),1/width)
        assert 0<=t[0]<=t[1]<=1
        assert breaks[region]<=r[0] and r[1]<=breaks[region+1]
        omt=1-t[1],1-t[0]
        phi[s]=[(Q(0),Q(0))]*49
        for j in range(7): phi[s][7*region+j]=scale(mul(power(t,j),power(omt,6-j)),comb(6,j))
    z=[isum(scale(phi[s][i],len(shells[s])) for s in shells) for i in range(49)]
    den=int(kernel['common_denominator']); assert den>0
    u=list(map(Q,origin['u']))
    e0=(Q(0),Q(0))
    for ray in kernel['blocks'][0]['rays']:
        assert ray['vector'][0]==0 and ray['weight']>0 and len(ray['vector'])==50
        a=isum(scale(add(zz,(uu,uu)),v) for v,zz,uu in zip(ray['vector'][1:],z,u))
        e0=add(e0,scale(square(a),Q(ray['weight'],den)))
    hist={}
    for s in shells:
        for t in shells:
            hist[s,t]=Counter(sum(a*b for a,b in zip(p,q)) for p in shells[s] for q in shells[t])
            assert sum(hist[s,t].values())==len(shells[s])*len(shells[t])
            assert all(hist[s,t][d]==hist[s,t][-d] for d in hist[s,t])
    angular=[]
    angular_shell_matrices={}
    for l,block in enumerate(kernel['blocks'][1:],1):
        assert block['dimension']==49
        assert all(ray['weight']>0 and len(ray['vector'])==49 for ray in block['rays'])
        if l%2:
            energy=(Q(0),Q(0))
        else:
            H={(s,t):sum(count*even_zonal(Q(d*d,4*s*t),l) for d,count in hist[s,t].items()) for s in shells for t in shells}
            angular_shell_matrices[str(l)]=[[str(H[s,t]) for t in shells] for s in shells]
            energy=(Q(0),Q(0))
            for ray in block['rays']:
                a={s:isum(scale(x,v) for x,v in zip(phi[s],ray['vector'])) for s in shells}
                energy=add(energy,scale(isum(scale(mul(a[s],a[t]),H[s,t]) for s in shells for t in shells),Q(ray['weight'],den)))
            # The exact total is nonnegative by PSD addition theorem.
            energy=max(Q(0),energy[0]),energy[1]
            assert energy[1]>=energy[0]
        angular.append({'degree':l,'interval':list(map(str,energy)), 'decimal':list(map(float,energy))})
    count=polynomial(origin['nested_count_quadratic']['coefficients'],[24]*3)
    assert len(origin['conditional_count_polynomials'])==2
    for row in origin['conditional_count_polynomials']:
        count+=Q(row['weight'])*polynomial(row['coefficients'],[24]*5)
    spherical=Q(0)
    assert len(origin['appended_spherical_local_energy_rays'])==2
    for ray in origin['appended_spherical_local_energy_rays']:
        assert dilation<=Q(ray['R']) and Q(ray['R'])**2<2*dilation**2
        H=sum(n*sum(Q(c)*zonal(Q(d,2),l) for l,c in enumerate(ray['scalar_coefficients'])) for d,n in hist[1,1].items())
        spherical+=Q(ray['weight'])*(24*Q(ray['B'])-H)
    angular_total=isum(tuple(map(Q,row['interval'])) for row in angular)
    gram_count=add(add(e0,angular_total),(count,count))
    subtotal=add(gram_count,(spherical,spherical))
    raw=Q(origin['objective'])+Q(4769,200000000)
    target=8*(raw-Q(5,8))
    if dilation==1:
        assert Q(4383,100000)<subtotal[0]<=subtotal[1]<Q(43833,1000000)
        assert subtotal[0]>target
    return dict(status='EXACT_D4_GRAM_COUNT_SPHERICAL_PASS',source=str(SOURCE.relative_to(PROJECT)),
        dilation=str(dilation),d4_shell_counts=[24,24,96],bits=192,E0=list(map(str,e0)),E0_decimal=list(map(float,e0)),
        angular=angular,angular_shell_matrices=angular_shell_matrices,angular_total=list(map(str,angular_total)),angular_total_decimal=list(map(float,angular_total)),
        count=str(count),count_decimal=float(count),spherical=str(spherical),spherical_decimal=float(spherical),
        gram_count=list(map(str,gram_count)),gram_count_decimal=list(map(float,gram_count)),
        subtotal=list(map(str,subtotal)),subtotal_decimal=list(map(float,subtotal)),raw=str(raw),raw_decimal=float(raw),
        target_credit=str(target),target_credit_decimal=float(target),
        scope='Completed E0; every angular Gram block including already embedded broad19 rays once; three count polynomials; two signed spherical energies. Other local energies not included here.')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--dilation',choices=['1','501473/500000'],default='1')
    args=ap.parse_args(); result=evaluate(Q(args.dilation))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k=='status' or k.endswith('_decimal')},indent=2))
