"""Independent 5D Fourier, guards, count conventions, PSD and exact rounding."""
from pathlib import Path
from fractions import Fraction as F
import json,math,sys
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
ctx.prec=384;D=Path(__file__).resolve().parent;s=json.loads((D/'witness.json').read_text());pre=json.loads((D/'prerequisites.json').read_text());source=json.loads((D/'witness_psd_exact.json').read_text())
q=lambda x:Q(x.numerator,x.denominator)
a=lambda x:arb(q(x))
assert (s['dimension'],s['p'],s['K'],s['J'])==(5,4,12,1) and F(s['R'])==F(3,2)
c=[F(str(x))for x in s['x']];assert len(c)==32
def rising(a,k):return math.prod((a+j for j in range(k)),start=F(1))
def binom(a,k):return math.prod((a-j for j in range(k)),start=F(1))/math.factorial(k)
half=F(5,2)
Pcs=[sum((c[k]*(-2)**j*binom(k+half-1,k-j)/math.factorial(j)/binom(k+half-1,k) for k in range(j,32)),F(0))for j in range(32)]
Hcs=[sum(((-1)**k*c[k]*(-2)**j*binom(k+half-1,k-j)/math.factorial(j)/binom(k+half-1,k) for k in range(j,32)),F(0))for j in range(32)]
assert Hcs==[sum((Pcs[k]*math.factorial(k)*(-1)**j*rising(half+j,k-j)/math.factorial(k-j)/math.factorial(j)for k in range(j,32)),F(0))for j in range(32)]
As=[[[F(z)for z in row]for row in A]for A in source['matrices']];assert len(As)==13;pivots=0
for ell,A in enumerate(As):
 n=6 if ell==0 else 5;assert len(A)==n and all(len(row)==n for row in A) and all(A[i][j]==A[j][i]for i in range(n)for j in range(n));S=[row[:]for row in A]
 for k in range(n):
  assert S[k][k]>0;pivots+=1
  for i in range(k+1,n):
   for j in range(k+1,n):S[i][j]-=S[i][k]*S[k][j]/S[k][k]
# Explicit external-data binding and conservative threshold monotonicity.
weights=list(map(lambda x:F(str(x)),s['count_weights']));bounds=list(map(F,s['count_bounds']));radii=list(map(F,s['count_radii']))
assert len(weights)==4 and min(weights)>=0 and bounds==[45,46,47,48]
published=list(map(F,['1.003013','1.008097','1.013259','1.017446']))
assert all(r==pub-F(1,10**6)for r,pub in zip(radii,published));count_cost=sum((w*b for w,b in zip(weights,bounds)),F(0))
def atan_bounds(x):
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1)for j in range(80));return lo,lo+x**161/161
l5,u5=atan_bounds(F(1,5));l239,u239=atan_bounds(F(1,239));pl=16*l5-4*u239;pu=16*u5-4*l239
P=Poly(list(map(q,Pcs)));H=Poly(list(map(q,Hcs)))
def chain(p):
 ss=[p,p.derivative()]
 while ss[-1].degree()>0:
  z=-(ss[-2]%ss[-1]);assert z;v=z.numer();ss.append(Poly(v/abs(v.content())))
 return ss
def variation(ss,t):
 vals=[p[p.degree()]if t is None else p(q(t))for p in ss];sg=[int(x>0)-int(x<0)for x in vals if x];return sum(x!=y for x,y in zip(sg,sg[1:]))
def no_roots(ss,lo,hi):assert variation(ss,lo)==variation(ss,hi)
def bernstein_upper(cs,lo,hi):
 degree=len(cs)-1;power=[sum((cs[k]*math.comb(k,j)*lo**(k-j)*(hi-lo)**j for k in range(j,degree+1)),F(0))for j in range(degree+1)]
 return max(sum((power[j]*F(math.comb(i,j),math.comb(degree,j))for j in range(i+1)),F(0))for i in range(degree+1))
g=pre['guards']['reciprocal_gaussian_guards'];eta=F(g['eta']);zeta=F(g['zeta']);guardc=F(g['c'])
assert eta==0 and zeta==0 and guardc==F(157,50) and F(5,8)<guardc<pl
cutoff=(10-guardc)/16;cost=eta*32*(10-guardc);low=H-q(cost)
assert low[0]>0 and low(q(cutoff))>0;no_roots(chain(low),F(0),cutoff)
assert H(q(cutoff))>0 and not g['fourier_intervals'];no_roots(chain(H),cutoff,None)
left=F(9,4)*pl;assert left>guardc;pc=chain(P);guard_reports=[]
for lo,hi in [tuple(map(F,z))for z in g['physical_intervals']]:
 assert left<lo<hi and P(q(left))<0 and P(q(lo))<0 and P(q(hi))<0;no_roots(pc,left,lo)
 upper=bernstein_upper(Pcs,lo,hi);allowance=a(eta*(lo-guardc))*(a(F(3,4)*lo)).exp();assert a(upper)<allowance
 guard_reports.append({'interval':[str(lo),str(hi)],'upper':str(a(upper)),'allowance':str(allowance)});left=hi
no_roots(pc,left,None);assert P(q(left))<0
norm=Hcs[0]-cost;assert norm>0 and norm==F(g['corrected_normalization'])
A=As[0];B=[row[1:]for row in A[1:]];b=[row[0]for row in A[1:]];m=len(b);aug=[B[i][:]+[b[i]]for i in range(m)]
for k in range(m):
 assert aug[k][k]>0
 for i in range(k+1,m):
  t=aug[i][k]/aug[k][k]
  for j in range(k,m+1):aug[i][j]-=t*aug[k][j]
x=[F(0)]*m
for i in reversed(range(m)):x[i]=(aug[i][m]-sum(aug[i][j]*x[j]for j in range(i+1,m)))/aug[i][i]
alpha=A[0][0]-sum(b[i]*x[i]for i in range(m));assert alpha==F(pre['origin_slack']) and alpha>0
W=[x]+[[F(i==j)for j in range(m)]for i in range(m)]
assert all(A[i][j]-(alpha if i==j==0 else 0)==sum(W[i][k]*B[k][ell]*W[j][ell]for k in range(m)for ell in range(m))for i in range(m+1)for j in range(m+1))
objective=(Pcs[0]+A[0][0]-alpha+count_cost+eta*guardc)/norm
assert objective==F(pre['saturated_objective']) and objective>0
upper=objective*pu**2/60;assert upper<F('.509996147')
out={'pass':True,'positive_PSD_pivots':pivots,'direct_Fourier_identity':True,'guard_intervals':guard_reports,'exact_count_cost':str(count_cost),'guarded_normalization':str(norm),'saturated_objective':str(objective),'density_interval':str(a(objective)*arb.pi()**2/60),'strict_upper':'.509996147','scope':'Independent exact radial/PSD/count arithmetic; full pair and cyclic/far replays and published count bounds remain dependencies'}
(D/'algebra_guards_audit.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items()if k not in ['guard_intervals','saturated_objective']})
# Export exact radial polynomials for independent pair verification.
(D/'radial_exact.json').write_text(json.dumps({'P':[str(x)for x in Pcs],'H':[str(x)for x in Hcs]},indent=2))
