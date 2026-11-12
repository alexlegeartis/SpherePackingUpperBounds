"""Exact reconstructed 9D radial pair verification. Angular proof is separate."""
if not __debug__:raise RuntimeError('Assertions required')
import sys,json,time,math,hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
sys.set_int_max_str_digits(0);ctx.prec=384
D=Path(__file__).parent;start=time.monotonic()
def q(x):return Q(str(F(str(x))))
def af(x):return arb(Q(x.numerator,x.denominator)) if isinstance(x,F) else arb(x)
def write(name,j):(D/name).write_text(json.dumps(j,indent=2)+'\n')
j=json.loads((D/'certificate.json').read_text());assert j['dimension']==9
P=Poly(list(map(Q,j['physical'])));H=Poly(list(map(Q,j['fourier'])));C=Q(j['C_exact']);lam=Q(j['lambda_count_exact']);assert lam>=0 and C>=0 and H(0)==1
# Reconstruct the exact normalized Laguerre Fourier relation independently.
x=Poly([0,1]);G=[Poly([1]),Q(9,2)-2*x]
for k in range(1,P.degree()):G.append(((2*k+Q(9,2)-2*x)*G[-1]-(k+Q(7,2))*G[-2])/(k+1))
G=[g/g(0) for g in G];rest=P;coeff=[]
for k in reversed(range(len(G))):
 c=rest[k]/G[k][k];coeff.append((k,c));rest-=c*G[k]
assert not rest
assert H==sum(((-1)**k*c*G[k] for k,c in coeff),Poly())
gates={}
for name,p,lo,sign in [('fourier',H,Q(0),1),('physical',P,Q(7),-1)]:
 rr=p.numer().complex_roots();assert sum(m for z,m in rr)==p.degree()
 for z,m in rr:
  if z.imag.contains(0):
   assert z.imag.is_exact() and z.imag==0
   assert z.real<arb(lo),(name,str(z.real))
 assert sign*p(lo)>0
 gates[name]=dict(start=str(lo),sign=sign,degree=p.degree(),root_count=sum(m for z,m in rr))
fit=json.loads((D/'fit.json').read_text());par=list(map(q,fit['parameters']));nv=len(fit['v_knots']);nl=len(fit['log_knots']);assert len(par)==2+nv+nl
angular=json.loads((D/'angular.json').read_text());B=Q(angular['B_exact']);assert B>0
VK=list(map(lambda v:F(str(v)),fit['v_knots']));LK=list(map(lambda v:F(str(v)),fit['log_knots']))
assert all(1<k<=2 for k in VK+LK)
assert par[0]+par[1+nv]>=0
for k in [F(1)]+VK:
 assert par[0]+sum((par[1+i] for i,x in enumerate(VK) if x<=k),Q(0))>=0
for k in sorted(set(VK+LK)):
 assert sum((q(k)*par[1+i] for i,x in enumerate(VK) if x==k),Q(0))+sum((par[2+nv+i] for i,x in enumerate(LK) if x==k),Q(0))>=0
endpoint=F(2);threshold=F('1.003012')**2;assert threshold==F(j['threshold_v_exact']);cuts=sorted(set([F(1),threshold,endpoint]+[a for a in VK+LK if 1<a<endpoint]))
scale=arb.pi()*arb(Q(9,8));ratio=arb(C/B)
AP=arb_poly([arb(z) for z in P]);AP1=AP.derivative();AP2=AP1.derivative()
def profile_parts(mid):
 a=par[0];b=par[1+nv];constant=-arb(par[0])
 for knot,c in zip(VK,par[1:1+nv]):
  if mid>knot:a+=c;constant-=arb(c)*af(knot)
 for knot,c in zip(LK,par[2+nv:]):
  if mid>knot:b+=c;constant-=arb(c)*af(knot).log()
 return arb(a),arb(b),constant
def exp_parts(v,abc):
 a,b,c=abc;phi=a*v+b*v.log()+c;xx=scale*v
 return xx,(xx-phi).exp(),(xx-2*phi).exp(),a+b/v,-b/v**2
def bound(lo,hi,abc,count):
 mid=(lo+hi)/2;rad=(hi-lo)/2;v=arb(af(mid),af(rad))
 xx,t1,t2,pp,p2=exp_parts(v,abc)
 shift=AP2(arb_poly([scale*af(mid),scale*af(rad)]));poly2=shift[0]
 for k in range(1,len(shift)):poly2+=arb(0,abs(shift[k]))
 f2=scale**2*poly2-ratio*(2*t1*((scale-pp)**2-p2)-t2*((scale-2*pp)**2-2*p2))
 if count:f2-=arb(lam)*scale**2*xx.exp()
 xx,t1,t2,pp,p2=exp_parts(af(mid),abc)
 f0=AP(xx)-ratio*(2*t1-t2);f1=scale*AP1(xx)-ratio*(2*t1*(scale-pp)-t2*(scale-2*pp))
 if count:f0-=arb(lam)*xx.exp();f1-=arb(lam)*scale*xx.exp()
 return f0+abs(f1)*af(rad)+abs(f2)*af(rad)**2/2
accepted=[];nodes=0
for lo,hi in zip(cuts,cuts[1:]):
 abc=profile_parts((lo+hi)/2);stack=[(lo,hi,0)]
 while stack:
  a,b,depth=stack.pop();nodes+=1
  assert nodes<200000 and time.monotonic()-start<600
  count=(lo+hi)/2<threshold;upper=bound(a,b,abc,count)
  if upper<0:accepted.append(dict(lo=str(a),hi=str(b),upper=str(upper),count_active=count));continue
  assert depth<40,('UNRESOLVED',str(a),str(b),str(upper))
  m=(a+b)/2;stack.extend([(m,b,depth+1),(a,m,depth+1)])
 print('PAIR',lo,hi,'leaves',len(accepted),flush=True)
leaves=sorted(accepted,key=lambda z:F(z['lo']));assert F(leaves[0]['lo'])==1 and F(leaves[-1]['hi'])==2
assert all(F(a['hi'])==F(b['lo']) for a,b in zip(leaves,leaves[1:]))
assert sum((F(z['hi'])-F(z['lo']) for z in leaves),F(0))==1
assert scale*2>7
# For v>=2 the physical function is negative and h in (0,1], so pair holds.
assert any(F(z['lo'])==threshold and not z['count_active'] for z in leaves)
assert all(not F(z['lo'])<threshold<F(z['hi']) and z['count_active']==(F(z['hi'])<=threshold) for z in leaves)
density=arb(P(0)+C+373*lam)*arb.pi()**arb(Q(9,2))*arb(Q(9,32))**arb(Q(9,2))/arb(Q(11,2)).gamma()
assert density<arb(q('0.190941549'))
write('partition.json',leaves)
result=dict(status='PASS',scope='Exact conditional radial Fourier, signs, pair and objective; angular proof and published average count373 separate',dimension=9,density_interval=str(density),strict_bound='0.190941549',gates=gates,nodes=nodes,leaves=len(leaves),elapsed=time.monotonic()-start,inputs={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['certificate.json','fit.json','angular.json']})
write('RESULT.json',result);print(json.dumps(result),flush=True)
