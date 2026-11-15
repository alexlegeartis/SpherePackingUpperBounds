#!/usr/bin/env python3
"""Exact origin credit from a nonnegative radial piece polynomial.

The result preserves every physical pair and cyclic expression.  It does not
convert an unverified or rejected source into a global packing bound.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,math

def trim(p):
    p=list(p)
    while p and p[-1]==0:p.pop()
    return p
def derivative(p):return trim([i*p[i] for i in range(1,len(p))])
def evaluate(p,x):
    y=Q(0)
    for c in reversed(p):y=y*x+c
    return y
def remainder(a,b):
    a=trim(a);b=trim(b);assert b
    while len(a)>=len(b):
        c=a[-1]/b[-1];off=len(a)-len(b)
        for i,v in enumerate(b):a[i+off]-=c*v
        a=trim(a)
    return a
def variations(values):
    s=[int(v>0)-int(v<0) for v in values if v]
    return sum(a!=b for a,b in zip(s,s[1:]))
def sturm_positive(p,lo,hi):
    p=trim(p)
    if not p or evaluate(p,lo)<=0 or evaluate(p,hi)<=0:return None
    if len(p)==1:return dict(endpoint_positive=True,root_count=0,sturm_length=1)
    seq=[p,derivative(p)]
    while seq[-1]:
        r=[-x for x in remainder(seq[-2],seq[-1])]
        if not r:break
        seq.append(r)
    left=variations([evaluate(s,lo) for s in seq]);right=variations([evaluate(s,hi) for s in seq])
    if left!=right:return None
    return dict(endpoint_positive=True,root_count=left-right,sturm_length=len(seq),
                left_variations=left,right_variations=right)
def bernstein_power(b):
    n=len(b)-1
    return [math.comb(n,k)*sum((-1)**(k-j)*math.comb(k,j)*b[j] for j in range(k+1)) for k in range(n+1)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('directory',type=Path)
    ap.add_argument('--baseline-factor',action='store_true');args=ap.parse_args()
    folder=args.directory.resolve();data=json.loads((folder/'exact_range.json').read_text());source=Path(data['source'])
    for n,h in data['source_sha256'].items():assert sha(source/n)==h
    W=[list(map(Q,row)) for row in data['W']];Y=[list(map(Q,row)) for row in data['Y']]
    k=data['rank'];degree=data['local_degree'];nc=degree+1;n=len(W);ns=n//nc
    factor=None
    if args.baseline_factor:
        assert k==13 and degree==4 and ns==5
        c=Q(253,500);aa=Q(1683024786574,10**20);bb=Q(125668197529249,10**20)
        p=[c*c,-2*c-c*c,1+2*c,Q(-1),Q(0)]
        power=[aa*p[i]+(bb*p[i-1] if i else 0) for i in range(5)]
        b=[sum(power[j]*Q(math.comb(i,j),math.comb(4,j)) for j in range(i+1)) for i in range(5)]
        z=[Q(0)]*5+b+[Q(0)]*3
        assert aa>=0 and bb>=0 and 0<=c<=1
        factor=dict(formula='(1-t)*(t-c)^2*(a+b*t)',t='4*r-6',c=str(c),a=str(aa),b=str(bb),zero_on_inner=True)
        eps=Q(0)
    else:
        numerical=json.loads((folder/'optimization.json').read_text())['passes'][-1]['coefficients']
        assert len(numerical)==k
        zbase=[Q(round(x*10**20),10**20) for x in numerical]
        accepted=None
        for eps in [Q(1,10**13),Q(3,10**13),Q(1,10**12),Q(3,10**12),Q(1,10**11),Q(3,10**11),Q(1,10**10)]:
            z=list(zbase)
            for j in range(10):z[j]+=eps
            coeff=[sum(x*y for x,y in zip(row,z)) for row in W]
            tests=[]
            for region in range(ns):
                lo,hi=map(Q,data['basis_intervals'][region]);left,right=map(Q,data['breaks'][region:region+2])
                power=bernstein_power(coeff[nc*region:nc*(region+1)])
                test=sturm_positive(power,(left-lo)/(hi-lo),(right-lo)/(hi-lo))
                if test is None:break
                tests.append(dict(region=region,power_coefficients=list(map(str,power)),left=str((left-lo)/(hi-lo)),right=str((right-lo)/(hi-lo)),**test))
            if len(tests)==ns:accepted=tests;break
        assert accepted is not None,'No exact positive radial weight found at the declared small shifts'
    a=[sum(x*y for x,y in zip(row,z)) for row in W]
    b=[sum(x*y for x,y in zip(row,z)) for row in Y]
    origin=json.loads((source/'origin.json').read_text());ker=json.loads((source/'kernel.json').read_text())
    den=int(ker['common_denominator']);M=[[0]*n for _ in range(n)]
    for ray in ker['blocks'][0]['rays']:
        assert ray['weight']>0
        vec=ray['vector'][1:]
        for i,x in enumerate(vec):
            for j,y in enumerate(vec):M[i][j]+=ray['weight']*x*y
    assert all(sum(Q(M[i][j],den)*b[j] for j in range(n))==a[i] for i in range(n))
    u=list(map(Q,origin['u']));v=list(map(Q,origin['v']));kap=Q(origin['kappa'])
    assert all(sum(Q(M[i][j],den)*u[j] for j in range(n))==v[i] for i in range(n))
    assert sum(x*y for x,y in zip(v,u))==kap
    gain=2*sum(x*y for x,y in zip(u,a))-sum(x*y for x,y in zip(a,b))
    assert gain>0
    unew=[x-y for x,y in zip(u,b)];vnew=[x-y for x,y in zip(v,a)];knew=kap-gain
    assert knew==sum(x*y for x,y in zip(vnew,unew))>=0
    objective=Q(origin['objective'])-gain/8
    threshold=Q(math.floor(objective*10**12)+1,10**12)
    assert objective<threshold
    if factor:
        assert a[:20]==[0]*20
        assert bernstein_power(a[20:])==power
        tests=[dict(region=4,proof='EXACT_FACTORIZATION_ON_CLOSED_[0,1]',**factor)]
    record=dict(status='EXACT_NONNEGATIVE_RADIAL_ORIGIN_CREDIT_VERIFIED',
        source=str(source),source_sha256=data['source_sha256'],
        source_global_feasibility='UNCHANGED; no sign condition is proved by this origin-only calculation',
        source_global_bound_proved=False,all_nonorigin_kernel_entries_unchanged=True,
        total_physical_origin_coupling_unchanged=True,all_pair_and_cyclic_expressions_unchanged=True,
        source_objective=origin['objective'],origin_credit=str(gain),density_credit=str(gain/8),
        new_formal_objective=str(objective),new_formal_threshold=str(threshold),
        source_objective_decimal=float(Q(origin['objective'])),density_credit_decimal=float(gain/8),
        new_formal_objective_decimal=float(objective),new_formal_threshold_decimal=float(threshold),
        epsilon_added_to_radial_weight=str(eps),functional_coefficients=list(map(str,z)),
        weight_Bernstein_coefficients=list(map(str,a)),range_solution=list(map(str,b)),
        new_psd_origin_u=list(map(str,unew)),new_psd_origin_v=list(map(str,vnew)),new_psd_kappa=str(knew),
        new_linear_energy='2*sum_x w(|x|)>=0, with w=a^T radial_features; zero added origin and nonorigin kernel',
        local_energy_identity='Enew(C)=Eold(C)-g=(u-b+z_C)^T M (u-b+z_C)+2*a^T*z_C+unchanged_other_local_energies',
        radial_nonnegativity=tests,
        exact_range_equations=n,origin_objective_identity=True)
    for file,h in data['source_sha256'].items():assert sha(source/file)==h
    (folder/'ORIGIN_CREDIT_VERIFIED.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({key:record[key] for key in ['status','density_credit_decimal','new_formal_objective_decimal','new_formal_threshold_decimal','epsilon_added_to_radial_weight','exact_range_equations']}))

if __name__=='__main__':main()
