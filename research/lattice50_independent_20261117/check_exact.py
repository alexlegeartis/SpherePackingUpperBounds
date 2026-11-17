from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import sys,json,hashlib
from flint import fmpq as Q,fmpq_poly as Poly
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;P=D/'producer';R=P;pub=D/'dependency'
r=json.load(open(P/'repaired.json'));cand=json.load(open(P/'input_candidate.json'));ann=json.load(open(P/'annular_final_certificate.json'))
for name in ['repaired.json','input_candidate.json','scalar48_kernel_upward_extensions.json']:
 assert hashlib.sha256((P/name).read_bytes()).hexdigest()==json.load(open(P/'FROZEN_HASHES.json'))[name]
source=pub/'rank2_candidate_B783489561_20261117.json';assert hashlib.sha256(source.read_bytes()).hexdigest()=='119dc3c5ef12da19e6edb01332e87e4c472da94d6c60c23fa438f74b7da147ae';co=[Q(str(F(str(x)))) for x in json.load(open(source))['even_coefficients']];co[0]-=Q(1,5000)
def basis(n,N):
 x=Poly([0,1]);b=[Poly([1]),x]
 for k in range(2,N+1):b.append(((2*k+n-4)*x*b[-1]-(k-1)*b[-2])/(k+n-3))
 assert all(p(Q(1))==1 for p in b)
 return b
b48=basis(48,2*len(co)-2);b50=basis(50,2*len(co)-2)
raw=sum((a*b48[2*i] for i,a in enumerate(co)),Poly([0]));rem=raw;exp=[Q(0)]*len(co)
for i in range(len(co)-1,-1,-1):
 k=2*i;exp[i]=rem[k]/b50[k][k];rem-=exp[i]*b50[k]
assert rem.is_zero() and min(exp)>=0 and exp[0]>0
normalized=[a/exp[0] for a in exp]
while normalized[-1]==0:normalized.pop()
saved=next(z for z in json.load(open(P/'scalar48_kernel_upward_extensions.json')) if z['n']==50)
assert normalized==list(map(Q,saved['coefficients']))
q=sum(normalized);assert q==Q(ann['q'])
# Reconstruct all Fourier-positive factors, rejecting unknown types.
f=Poly([1]);counts={}
for z in r['factors']:
 typ=z['type'];counts[typ]=counts.get(typ,0)+1
 if typ=='negative_linear':
  v=Q(z['a']);assert v>0;factor=Poly([1,1/v])
 elif typ=='positive_double':
  v=Q(z['mid']);assert v>0;factor=Poly([1,-1/v])**2
 elif typ=='complex_pair':
  re,im=Q(z['re']),Q(z['im']);assert im>0;den=re*re+im*im;factor=(Poly([-re,1])**2+im*im)/den
 else:raise AssertionError(typ)
 f*=factor
assert f==Poly(list(map(Q,r['Q']))) and f(0)==1 and f.degree()==71
# Independent direct Gaussian monomial transform, without Laguerre recurrence.
p=Poly(list(map(Q,r['P'])))
ft=[]
for j in range(72):
 ft.append(sum((p[k]*(-1)**j*factorial(k)*comb(k+24,k-j)/factorial(j) for k in range(j,72)),Q(0)))
assert Poly(ft)==f
# Exact primal Schur normalization and restored primitive-line bound.
a=F(str(cand['Z'][1][1]));z=F(str(cand['Z'][0][1]))-F(1,10**9);s=z*z/a+F(1,10**12)
rho=F(1,10**70);L=(1+rho)/(1-rho)
assert a>0 and s*a-z*z>0 and a*L+z<0
assert (a,z,s,L)==tuple(F(ann[k]) for k in ['a','z','s','L'])
assert F(168**168,factorial(168))>1/rho
# Recompute each exact spherical polynomial cap from its moment coefficients.
budget=F(0);last=F(22);ks=[]
for row in ann['rows']:
 l,r=F(row['l']),F(row['r']);h=F(row['height']);assert l==last and r-l==F(1,4) and h>=0;last=r
 radmax=r/22
 if radmax<2:C=1
 else:
  ss=(1-2/radmax)**2;mu=F(1);best=None;bestk=None
  for k in range(501):
   ratio=F(2*k+1,50+2*k);c0=mu*(ratio-ss)
   if c0>0:
    value=(1-ss)/c0
    if best is None or value<best:best=value;bestk=k
   mu*=ratio
  assert best is not None;C=best.numerator//best.denominator;ks.append(bestk)
 assert C==row['cap'] and C>=1
 budget+=2**51*C*h
assert len(ann['rows'])==632 and last==180 and budget==F(ann['budget'])
value=F(22,4)**25/F(factorial(25))*(F(str(p[0]))+s+budget)
safe=F('1.75102693e-7');assert 0<value<safe
# Producer and independent interval replays may have slightly different outward heights;
# each final rational objective must satisfy the same announced strict bound.
pa=json.load(open(R/'original_annular_final_certificate.json'));pv=F(22,4)**25/F(factorial(25))*(F(str(p[0]))+s+F(pa['budget']))
assert pv==F(json.load(open(R/'density_result.json'))['density_upper_exact']) and pv<safe
out=dict(status='PASS',dimension=50,safe_upper=str(safe),density_decimal=float(value),factor_counts=counts,annuli=632,budget_decimal=float(budget),max_selected_cap_k=max(ks),source_kernel_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_repo_commit='af863e2014561e060af7ee2b0eecffc468f4eca3',published_scalar_dependency=True,standalone_new_rank_two_replay=False)
(D/'exact_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
