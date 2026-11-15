"""Build literal V7Expression ray-based smooth-core syntax for root0.

The output remains a candidate until reflected against the owner's Lean AST.
No precontracted Gram matrices or externally asserted polynomial identities.
"""
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

DATA=Path(__file__).resolve().parent.parent/'data/v7'

def build(dag,lam,a):
    kd,od,sd=[json.loads((DATA/n).read_text()) for n in ('kernel.json','origin.json','scalar.json')]
    c=dag.rat;add=dag.add;mul=dag.mul;neg=dag.neg;div=dag.div;power=dag.pow
    def sub(x,y):return add(x,neg(y))
    def total(xs):
        r=c(0)
        for x in reversed(list(xs)):r=add(x,r)
        return r
    def basis(i,r):
        lo,hi=map(Q,kd['radial_basis_intervals'][0]);z=div(sub(r,c(lo)),sub(c(hi),c(lo)))
        return mul(mul(c(comb(6,i)),power(z,i)),power(sub(c(1),z),6-i))
    radial_cache={}
    def radial(v,r):
        key=(tuple(map(Q,v)),r)
        if key not in radial_cache:radial_cache[key]=total(mul(c(v[i]),basis(i,r)) for i in range(7))
        return radial_cache[key]
    cheb_cache={}
    def cheb(n,t):
        key=(n,t)
        if key not in cheb_cache:cheb_cache[key]=c(1) if n==0 else mul(c(2),t) if n==1 else sub(mul(mul(c(2),t),cheb(n-1,t)),cheb(n-2,t))
        return cheb_cache[key]
    def zonal(n,t):return div(cheb(n,t),c(n+1))
    def gram(r,s,t):
        return total(div(mul(zonal(ell,t),total(mul(mul(c(ray['weight']),radial(ray['vector'][int(ell==0):],r)),radial(ray['vector'][int(ell==0):],s)) for ray in block['rays'])),c(kd['common_denominator'])) for ell,block in enumerate(kd['blocks']))
    counts=[dict(od['nested_count_quadratic'],weight=1)]+od['conditional_count_polynomials']
    def matrix(row):
        n=len(row['radii']);cs=row['coefficients'];mat=[[Q(0)]*n for _ in range(n)];k=1+n
        for i in range(n):
            for j in range(i,n):mat[i][j]=mat[j][i]=Q(cs[k]);k+=1
        return mat
    count_pair=total(mul(c(row['weight']),total(total(c(x) for x in r) for r in matrix(row))) for row in counts)
    def sph(t):return total(mul(neg(c(row['weight'])),total(mul(c(x),zonal(i,t)) for i,x in enumerate(row['scalar_coefficients']))) for row in od['appended_spherical_local_energy_rays'])
    def core(r,s,t):return add(add(gram(r,s,t),sph(t)),count_pair)
    def linear(r):
        out=radial(od['v'],r)
        for values in (
            [div(mul(mul(c(row['weight']),c(row['M'])),radial(row['w'],r)),c(2)) for row in od['pointed_caps']],
            [div(mul(mul(c(row['weight']),c(row['capacity'])),radial(row['w'],r)),c(2)) for row in od['pointed_tents']],
            [div(mul(c(row['weight']),radial(row['w'],r)),c(2)) for row in od['independent_antipodal_allocations']],
            [div(mul(c(row['weight']),c(row['B'])),c(2)) for row in od['appended_spherical_local_energy_rays']],
            [mul(c(row['weight']),total(c(x) for x in row['coefficients'][1:1+len(row['radii'])])) for row in counts]):
            out=add(out,total(values))
        return out
    def scalar(r):
        x=mul(c(2),power(r,2));poly=c(0)
        for coeff in reversed(sd['F']):poly=add(c(coeff),mul(x,poly))
        return mul(dag.node('exp',mul(neg(c(2)),power(r,2))),poly)
    def q0(r):return add(add(scalar(r),mul(c(2),linear(r))),core(r,r,c(1)))
    r,s,t=[dag.node('var',i) for i in range(3)]
    return sub(sub(neg(core(r,s,t)),div(add(q0(r),q0(s)),c(46))),mul(c(lam),sub(power(t,4),mul(power(c(a),2),power(t,2)))))
