"""Fixed-V7 pointed and graph energy on the actual normalized D4 star.

No search: exact shell incidence, rational interval evaluation of radicals.
Run from the project root with the project Python. Output is package-local.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product, combinations
from math import comb, isqrt
from pathlib import Path
import json
import argparse

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[2] / 'data/v7'  # the D4 package data folder
BITS = 160

def pt(x):
    return Q(x), Q(x)

def add(a,b):
    return a[0]+b[0],a[1]+b[1]

def neg(a):
    return -a[1],-a[0]

def mul(a,b):
    values=[x*y for x in a for y in b]
    return min(values),max(values)

def div(a,b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))

def scale(a,x):
    return mul(a,pt(x))

def root(x):
    x=Q(x); assert x>=0
    nr,dr=isqrt(x.numerator),isqrt(x.denominator)
    if nr*nr==x.numerator and dr*dr==x.denominator:
        return pt(Q(nr,dr))
    d=1<<BITS
    a=isqrt(x.numerator*d*d//x.denominator)
    low=Q(a,d)
    high=low if low*low==x else Q(a+1,d)
    assert low*low<=x<=high*high
    return low,high

def power(a,n):
    result=pt(1)
    for _ in range(n): result=mul(result,a)
    return result

def isum(items):
    result=pt(0)
    for item in items:result=add(result,item)
    return result

def clip(a,lo=Q(0),hi=Q(1)):
    return max(lo,min(hi,a[0])),max(lo,min(hi,a[1]))

def sign_radical(a,n,b):
    """Exact sign of a*sqrt(n)+b for rational a,b and n>=0."""
    a,b,n=Q(a),Q(b),Q(n)
    if not a:return (b>0)-(b<0)
    if a<0:return -sign_radical(-a,n,-b)
    if b>=0:return 1 if n or b else 0
    diff=a*a*n-b*b
    return (diff>0)-(diff<0)

def probability(m,n,h,cap):
    """P(x|y), |x|²=m, |y|²=n, x.y=h/2, in dimension four.

    On the transverse 2-sphere a coordinate is uniform on [-1,1].
    c=(m+d²-a²)/(2d sqrt(m)), t=h/(2sqrt(mn)).
    """
    d,a2,sigma=map(Q,(cap['d'],cap['a2'],cap['sigma']))
    t=div(pt(h),scale(root(m*n),2))
    c=div(pt(m+d*d-a2),scale(root(m),2*d))
    A=add(c,neg(scale(t,sigma)))
    D2=(1-sigma*sigma)*(1-Q(h*h,4*m*n))
    assert D2>=0
    if not D2:
        # The closed event includes equality A=0.
        s=sign_radical(m+d*d-a2,n,-sigma*d*h)
        return pt(int(s<=0))
    return clip(scale(add(pt(1),neg(div(A,root(D2)))),Q(1,2)))

def features(m,kernel):
    r=root(m); breaks=list(map(Q,kernel['breaks']))
    assert r[0]>=breaks[0] and r[1]<=breaks[-1]
    region=next(i for i in range(len(breaks)-1)
                if breaks[i]<=r[0]<=r[1]<=breaks[i+1])
    a,b=map(Q,kernel['radial_basis_intervals'][region]);degree=kernel['local_degree']
    t=scale(add(r,pt(-a)),1/(b-a));om=add(pt(1),neg(t))
    assert 0<=t[0]<=t[1]<=1
    z=[pt(0)]*((degree+1)*(len(breaks)-1))
    for j in range(degree+1):
        z[region*(degree+1)+j]=scale(mul(power(t,j),power(om,degree-j)),comb(degree,j))
    return z

def radial_weight(m,row,kernel):
    if 'pivot_lower_cutoff' in row:
        cutoff=Q(row['pivot_lower_cutoff'])
        if m<=cutoff*cutoff:return pt(0)
    w=list(map(Q,row['w']));z=features(m,kernel)
    assert len(w)==len(z) and min(w)>=0
    return isum(scale(a,b) for a,b in zip(z,w))

def encoded(interval):
    # Round once outwards for a compact portable certificate. Internal
    # arithmetic retains the full rational interval.
    scale=1<<128
    low=Q(interval[0].numerator*scale//interval[0].denominator,scale)
    high=Q(-((-interval[1].numerator*scale)//interval[1].denominator),scale)
    assert low<=interval[0]<=interval[1]<=high
    return {'interval':list(map(str,(low,high))),'decimal_diagnostic':list(map(float,interval))}

def evaluate(dilation=Q(1)):
    dilation=Q(dilation)
    assert dilation in (Q(1),Q(501473,500000))
    origin=json.loads((SOURCE/'origin.json').read_text())
    kernel=json.loads((SOURCE/'kernel.json').read_text())
    graphs=json.loads((SOURCE/'graph_energies.json').read_text())
    assert kernel['dimension']==4 and kernel['local_degree']==6
    lattice=[q for q in product(range(-2,3),repeat=4)
             if sum(q)%2==0 and 0<sum(x*x for x in q)<=6]
    shells={m*dilation**2:[q for q in lattice if sum(x*x for x in q)==2*m] for m in (1,2,3)}
    assert list(map(len,shells.values()))==[24,24,96]
    assert Q(3)*dilation**2<Q(7,4)**2<Q(4)*dilation**2
    assert min(sum((a-b)**2 for a,b in zip(x,y)) for x,y in combinations([(0,0,0,0)]+lattice,2))==2
    histogram=Counter((m,n,dilation**2*sum(a*b for a,b in zip(x,y)))
                      for m,xs in shells.items() for n,ys in shells.items() for x in xs for y in ys)
    assert sum(histogram.values())==144**2
    caps=[]
    for index,cap in enumerate(origin['pointed_caps']):
        by_shell=[]
        for n,ys in shells.items():
            occupancy=isum(scale(probability(m,n,h,cap),count)
                           for (m,nn,h),count in histogram.items() if nn==n)
            slack=add(pt(cap['M']*len(ys)),neg(occupancy))
            # Positivity is a previously certified geometry theorem; this
            # interval calculation does not clamp away a negative endpoint.
            by_shell.append(scale(mul(radial_weight(n,cap,kernel),slack),Q(cap['weight'])))
        caps.append({'index':index,**encoded(isum(by_shell))})
    tents=[]
    for index,tent in enumerate(origin['pointed_tents']):
        d=Q(tent['d']); assert Q(tent['sigma'])==-1
        by_shell=[]
        for n,ys in shells.items():
            def phi(m,h):
                # |x+d*yhat|²=m+d²+d*h/sqrt(n).
                raw=add(pt(1-2*(m+d*d)),neg(div(pt(2*d*h),root(n))))
                return max(Q(0),raw[0]),max(Q(0),raw[1])
            mass=isum(scale(phi(m,h),count) for (m,nn,h),count in histogram.items() if nn==n)
            slack=add(pt(Q(tent['capacity'])*len(ys)),neg(mass))
            by_shell.append(scale(mul(radial_weight(n,tent,kernel),slack),Q(tent['weight'])))
        tents.append({'index':index,**encoded(isum(by_shell))})
    # Full-pair count gauge = sum(w)*(1-N_near/24), hence exactly zero.
    gauge=origin['full_pair_count_gauge']
    assert gauge['gauge']=='indicator_over_count_bound' and gauge['count_bound']==24
    assert 1<Q(gauge['radius'])**2<2
    assert sum(len(v) for m,v in shells.items() if m<=Q(gauge['radius'])**2)==24
    conditional=[]
    for row in origin['conditional_radial_rays']:
        assert row['kind']=='one_outer' and row['maximum']==24
        counts=[sum(len(v) for m,v in shells.items() if m<=Q(row[key])**2)
                for key in ('inner_radius','outer_radius')]
        assert counts==[24,24]
        value=Q(row['weight'])*(23-counts[0])*(counts[1]-counts[0])
        assert value==0
        conditional.append(str(value))
    assert len(conditional)==2
    # Both unshifted allocation profiles have capacity(0)=1 and vanish at
    # distance²>=1. Every pivot y has -y in this centrally symmetric star.
    allocations=origin['independent_antipodal_allocations']
    assert len(allocations)==2
    assert {x['capacity_profile'] for x in allocations}=={'classical_improvement','two_level'}
    for row in allocations:
        assert Q(row.get('shift','0'))==Q(row.get('delta','0'))==0
        if row['capacity_profile']=='two_level':
            assert Q(row['inner_radius'])>0 and Q(row['outer_radius'])<1
    assert all(tuple(-a for a in x) in lattice for x in lattice)
    # Graph energies weight*(ordered capacity - ordered selected edges).
    # Active slots use negative angular intervals. Membership below is exact
    # via squared inequalities with the sign checked separately.
    def cosine_between(h,m,n,lo,hi):
        assert -1<=lo<=hi<0
        return h<0 and 4*m*n*lo*lo>=h*h>=4*m*n*hi*hi
    def count_edges(radial,angular):
        a,b=radial;lo,hi=angular
        return sum(count for (m,n,h),count in histogram.items()
                   if a*a<=m<=b*b and a*a<=n<=b*b and cosine_between(h,m,n,lo,hi))
    assert list(map(Q,graphs['weights'][:2]))==[0,0]
    specs=[('broad_antipodal',Q(graphs['weights'][2]),1358,(Q(33,28),Q(7,4)),(Q(-1),Q(-7,8))),
           ('near_antipodal',Q(graphs['outer_graph_weight']),120,(Q(1),Q(103,100)),(Q(-1),Q(-3,4))),
           ('narrow_antipodal',Q(graphs['narrow_graph_weight']),604,(Q(33,28),Q(7,4)),(Q(-1),Q(-183,200))),
           ('moderate',Q(graphs['moderate_graph_weight']),1880,(Q(7,5),Q(3,2)),(Q(-4,5),Q(-11,20)))]
    graph_rows=[]
    for name,weight,capacity,radial,angular in specs:
        edges=count_edges(radial,angular)
        assert 0<=edges<=capacity
        graph_rows.append({'name':name,'weight':str(weight),'capacity':capacity,'ordered_edges':edges,'energy':str(weight*(capacity-edges))})
    graph_total=sum((Q(x['energy']) for x in graph_rows),Q(0))
    cap_total=isum(tuple(map(Q,x['interval'])) for x in caps)
    tent_total=isum(tuple(map(Q,x['interval'])) for x in tents)
    total=add(add(cap_total,tent_total),pt(graph_total))
    assert 0<total[0]<=total[1]
    if dilation==1:assert total[1]<Q(651224,100000000)
    assert total[1]-total[0]<Q(1,10**30)
    report={'status':'EXACT_D4_POINTED_GRAPH_ENERGY_ENCLOSURE_PASS','source':str(SOURCE),'dilation':str(dilation),'sqrt_bits':BITS,
            'caps':caps,'cap_total':encoded(cap_total),'tents':tents,'tent_total':encoded(tent_total),
            'full_pair_count_gauge':'0','independent_antipodal_allocations':['0','0'],
            'conditional_radial_rays':conditional,
            'graphs':graph_rows,'graph_total':str(graph_total),'total':encoded(total),
            'scope':'Actual normalized D4 star. Pointed, tent, graph, two independent antipodal and full-pair gauge components only; no density theorem is claimed.'}
    return report

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dilation',choices=['1','501473/500000'],default='1')
    args=ap.parse_args()
    report=evaluate(Q(args.dilation))
    name='pointed_graph_energy.json' if args.dilation=='1' else 'pointed_graph_energy_r1.json'
    (HERE/name).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','cap_total','tent_total','graph_total','total']},indent=2))

if __name__=='__main__':main()
