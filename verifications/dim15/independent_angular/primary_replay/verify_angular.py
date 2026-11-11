"""Exact angular PSD, one-sign-transition and full curved-profile cyclic gate.
Discovery input is rationalized once; every acceptance below uses exact
rational arithmetic or rigorous Arb interval enclosures, never floating signs.
"""
if not __debug__:
    raise RuntimeError('Assertions are required; do not run this verifier with -O or -OO.')

import json,math,hashlib,time
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq,fmpq_poly,arb,arb_poly,ctx
ctx.prec=256
D=Path(__file__).resolve().parent
ff=lambda z:fmpq(F(z).numerator,F(z).denominator)
X=fmpq_poly([0,1]);n=15;a=fmpq(12);b=fmpq(16)
source=D/'a12_b16.json';j=json.loads(source.read_text());c=[ff(float(v)) for v in j['coefficients']]
assert c[0]==1 and min(c)>=0
guard=fmpq(1,10000);c[0]-=guard;c=[z/(1-guard) for z in c]
G=[fmpq_poly([1]),X]
for l in range(1,len(c)-1):G.append(((2*l+n-2)*X*G[-1]-l*G[-2])/(l+n-2))
Q=sum((z*g for z,g in zip(c,G)),fmpq_poly());B=sum(c)
cert={'dimension':n,'a_exact':str(a),'b_exact':str(b),'guard_exact':str(guard),'gegenbauer_coefficients':list(map(str,c)),'polynomial_coefficients':list(map(str,Q.coeffs())),'B_exact':str(B),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
(D/'angular_exact.json').write_text(json.dumps(cert,indent=2)+'\n')
def sturm(P):
 seq=[P,P.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);assert rem
  num=rem.numer();seq.append(fmpq_poly(num/abs(num.content())))
 return seq
def var(seq,t):
 sg=[int(v>0)-int(v<0) for p in seq if (v:=p(t))!=0]
 return sum(x!=y for x,y in zip(sg,sg[1:]))
seq=sturm(Q);assert Q(-1)<0 and Q(fmpq(1,2))<0 and Q(1)>0
assert var(seq,-1)==var(seq,fmpq(1,2))
assert var(seq,fmpq(1,2))-var(seq,1)==1
margin=fmpq(1,20000);cap=Q+margin;seq2=sturm(cap)
assert cap(-1)<0 and cap(fmpq(1,2))<0 and var(seq2,-1)==var(seq2,fmpq(1,2))
print('EXACT ANGULAR PASS',float(B),'cap margin',str(margin),flush=True)
AQ=arb_poly([arb(z) for z in Q.coeffs()]);DQ=AQ.derivative()
def qeval(P,v):
 lo=max(arb(-1),v.lower());hi=min(arb(1),v.upper());assert lo<=hi
 mid=(lo+hi)/2;rad=(hi-lo)/2
 cc=P(arb_poly([mid,rad]));err=sum((abs(cc[k]) for k in range(1,len(cc))),arb(0))
 return cc[0]+arb(0,err.upper())
def fun(T,U,der=False):
 d=1-2*T*U;e=U*U-2*T*U;W=(1+e).sqrt()
 if not W.is_finite():return None
 t1=(1-U*T)/W;t2=(U-T)/W
 E1=-arb(a)*d-arb(b)*d*(2*U*U-2+d);E2=-arb(a)*e-arb(b)*e*e
 r1=E1.exp();r2=E2.exp();Q0=qeval(AQ,T);Q1=qeval(AQ,t1);Q2=qeval(AQ,t2)
 V=Q0+r1*Q1+r2*Q2
 if not der:return V
 dt=-2*U;du=-2*T;et=-2*U;eu=2*(U-T)
 wt=et/(2*W);wu=eu/(2*W)
 t1t=-U/W-t1*wt/W;t1u=-T/W-t1*wu/W
 t2t=-1/W-t2*wt/W;t2u=1/W-t2*wu/W
 E1t=-(arb(a)+2*arb(b)*e)*dt
 E1u=-(arb(a)+2*arb(b)*e)*du-4*arb(b)*d*U
 E2t=-(arb(a)+2*arb(b)*e)*et;E2u=-(arb(a)+2*arb(b)*e)*eu
 Vt=qeval(DQ,T)+r1*(E1t*Q1+qeval(DQ,t1)*t1t)+r2*(E2t*Q2+qeval(DQ,t2)*t2t)
 Vu=r1*(E1u*Q1+qeval(DQ,t1)*t1u)+r2*(E2u*Q2+qeval(DQ,t2)*t2u)
 return V,Vt,Vu
start=time.time();stack=[(fmpq(-1),fmpq(1,2),fmpq(1),fmpq(4),0,0)];leaves=skipped=0;maxdepth=0;tested=0
while stack:
 tl,th,ul,uh,ld,md=stack.pop();tested+=1
 if 2*tl*ul>1:skipped+=1;continue
 T=arb(tl).union(arb(th));U=arb(ul).union(arb(uh));got=fun(T,U,True);okay=False
 if got is not None:
  V,Vt,Vu=got;assert all(z.is_finite() for z in got)
  if V<0:okay=True
  else:
   tm=(tl+th)/2;um=(ul+uh)/2
   center=fun(arb(tm),arb(um));assert center.is_finite()
   upper=center+abs(Vt)*arb((th-tl)/2)+abs(Vu)*arb((uh-ul)/2)
   okay=upper<0
 if okay:
  leaves+=1;maxdepth=max(maxdepth,ld,md);continue
 # Check corners first to expose a real positive lobe before deep subdivision.
 if max(ld,md)>=8:
  for tt in [tl,th]:
   for uu in [ul,uh]:
    if 2*tt*uu<=1:
     z=fun(arb(tt),arb(uu))
     assert not z>0,('POSITIVE',str(tt),str(uu),str(z))
 assert max(ld,md)<28,('DEPTH',tl,th,ul,uh)
 if ld<=md:
  mid=(tl+th)/2;stack.extend([(tl,mid,ul,uh,ld+1,md),(mid,th,ul,uh,ld+1,md)])
 else:
  mid=(ul+uh)/2;stack.extend([(tl,th,ul,mid,ld,md+1),(tl,th,mid,uh,ld,md+1)])
 if tested%5000==0:print('PROGRESS',tested,leaves,skipped,len(stack),'seconds',time.time()-start,flush=True)
# For u>=4 and t>=-1/4, second cosine <=1/u-t<=1/2;
# only third term may be positive; z>=u², so r2<=exp(-12*15-16*225).
# For t<=-1/4, d>=3, r1<=exp(-36), and r2 has same small bound.
small=(-arb(12*15+16*225)).exp()
assert arb(B)*small<arb(margin)
assert arb(B)*((-arb(36)).exp()+small)<arb(margin)
out={'pass':True,'dimension':15,'precision_bits':ctx.prec,'cap_margin_exact':str(margin),'single_sign_transition':True,'certified_rectangles':leaves,'discarded_sorting_rectangles':skipped,'max_depth':maxdepth,'tested':tested,'all_u_ge4_tail_pass':True,'seconds':time.time()-start,'angular_sha256':hashlib.sha256((D/'angular_exact.json').read_bytes()).hexdigest()}
(D/'angular_verified.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
