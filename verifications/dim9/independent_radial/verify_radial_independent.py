"""Separate 9D monomial Fourier/Sturm/first-derivative pair audit."""
if not __debug__: raise RuntimeError('Assertions required')
import json, math, sys, time, hashlib
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_poly,fmpq,fmpq_poly,ctx
ctx.prec=384;sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;start=time.monotonic()
def q(x):return fmpq(str(F(x)))
def a(x):return arb(q(x))
j=json.loads((D/'certificate.json').read_text());assert j['dimension']==9
p=list(map(F,j['physical']));h=list(map(F,j['fourier']));C=F(j['C_exact']);lam=F(j['lambda_count_exact']);assert C>=0 and lam>=0
# Fourier(X^j exp(-X)) = j! L_j^(4)(Y) exp(-Y).
def binom(t,k):
 r=F(1)
 for i in range(k):r*=F(t-i,i+1)
 return r
hcheck=[sum((p[jj]*math.factorial(jj)*binom(F(2*jj+7,2),jj-k)*F((-1)**k,math.factorial(k)) for jj in range(k,len(p))),F(0)) for k in range(len(p))]
while len(hcheck)>1 and hcheck[-1]==0:hcheck.pop()
assert h==hcheck
P=fmpq_poly(list(map(q,p)));H=fmpq_poly(list(map(q,h)));assert H(0)==1

def variation(signs):
 signs=[z for z in signs if z];return sum(x!=y for x,y in zip(signs,signs[1:]))
def sign(x):return 1 if x>0 else -1 if x<0 else 0

def ray(poly,lo,want):
 seq=[poly,poly.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);assert rem
  # Positive rescaling keeps signs while controlling rational growth.
  rem=rem/abs(rem.leading_coefficient());seq.append(rem)
 v0=variation([sign(z(q(lo))) for z in seq]);vi=variation([sign(z.leading_coefficient()) for z in seq])
 assert v0==vi and want*poly(q(lo))>0
 return dict(start=str(lo),variations_at_start=v0,variations_at_infinity=vi)
rays={'fourier':ray(H,F(0),1),'physical':ray(P,F(7),-1)}
print('CLOSED FORM FOURIER AND STURM RAYS PASS',flush=True)
fit=json.loads((D/'fit.json').read_text(),parse_float=str);par=list(map(F,fit['parameters']));nv=len(fit['v_knots']);nl=len(fit['log_knots']);assert len(par)==2+nv+nl
ang=json.loads((D/'angular.json').read_text());B=F(ang['B_exact']);assert ang['dimension']==9 and B>0
VK=list(map(F,fit['v_knots']));LK=list(map(F,fit['log_knots']));assert all(1<k<=2 for k in VK+LK)
assert par[0]+par[1+nv]>=0
slope=par[0];assert slope>=0
for k in sorted(set(VK+LK)):
 da=sum((v for x,v in zip(VK,par[1:1+nv]) if x==k),F(0));db=sum((v for x,v in zip(LK,par[2+nv:]) if x==k),F(0));assert k*da+db>=0;slope+=da;assert slope>=0
threshold=F('1.003012')**2;assert threshold==F(j['threshold_v_exact'])
AP=arb_poly(list(map(a,p)));DP=AP.derivative();scale=arb.pi()*a(F(9,8));ratio=a(C/B)
def polynomial(poly,x):
 mid=(x.lower()+x.upper())/2;rad=(x.upper()-x.lower())/2
 terms=poly(arb_poly([mid,rad]));return terms[0]+arb(0,sum((abs(z) for z in terms.coeffs()[1:]),arb(0)).upper())
def evaluate(v,linear,logarithm,constant,count):
 phi=a(linear)*v+a(logarithm)*v.log()+constant;sl=a(linear)+a(logarithm)/v
 X=scale*v;e1=(X-phi).exp();e2=(X-2*phi).exp()
 return polynomial(AP,X)-ratio*(2*e1-e2)-(a(lam)*X.exp() if count else 0),scale*polynomial(DP,X)-ratio*(2*e1*(scale-sl)-e2*(scale-2*sl))-(a(lam)*scale*X.exp() if count else 0)
knots=sorted(set([F(1),threshold,F(2)]+[x for x in VK+LK if 1<x<2]));leaves=[];checks=0
for lo,hi in zip(knots,knots[1:]):
 mid=(lo+hi)/2;linear=par[0];logarithm=par[1+nv];constant=-a(par[0])
 for b,k in zip(par[1:1+nv],VK):
  if k<mid:linear+=b;constant-=a(b*k)
 for b,k in zip(par[2+nv:],LK):
  if k<mid:logarithm+=b;constant-=a(b)*a(k).log()
 stack=[(lo,hi,0)]
 while stack:
  l,r,depth=stack.pop();checks+=1;assert checks<1000000 and time.monotonic()-start<900
  v=a(l).union(a(r));value,deriv=evaluate(v,linear,logarithm,constant,hi<=threshold)
  if not value<0:
   center=evaluate(a((l+r)/2),linear,logarithm,constant,hi<=threshold)[0]
   value=center+abs(deriv)*a((r-l)/2)
  if value<0:leaves.append((l,r));continue
  assert depth<50
  m=(l+r)/2;stack.extend([(l,m,depth+1),(m,r,depth+1)])
 print('PAIR',lo,hi,len(leaves),flush=True)
leaves.sort();assert leaves[0][0]==1 and leaves[-1][1]==2 and all(x[1]==y[0] for x,y in zip(leaves,leaves[1:]));assert scale*2>7
# Beyond v=2, P is negative; phi>=0 gives 2h-h^2>=0.
value=a(p[0]+C+373*lam)*arb.pi()**a(F(9,2))/a(F(11,2)).gamma()*a(F(9,32))**a(F(9,2))
assert value<a('0.190926892')
out=dict(status='PASS',dimension=9,strict_density_upper_bound='0.190926892',density_interval=str(value),independent_fourier_closed_form=True,sturm_rays=rays,pair_intervals=len(leaves),pair_checks=checks,average_neighbor_count_external=373,seconds=time.monotonic()-start,sha256={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['certificate.json','fit.json','angular.json']})
(D/'RADIAL_INDEPENDENT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
