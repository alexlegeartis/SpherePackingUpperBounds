import sys;sys.set_int_max_str_digits(0)
import json,math,time
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=1024
if not __debug__: raise RuntimeError('Do not run this verifier with -O')
from pathlib import Path
p=str(Path(__file__).resolve().parent)+'/'
z=json.load(open(p+'input_candidate.json'));raw=json.load(open(p+'repaired.json'));PP=list(map(F,raw['P']));co=list(map(F,z['co']));q=F(z['q']);assert sum(co)==q and co[0]==1 and min(co)>=0
Bs=[[[F(x)for x in row]for row in B]for B in z['B']];zz=list(map(F,z['z']));ss=F(z['s']);LL=F(1)
assert len(Bs)==len(co)
import sympy as sy
for B in Bs:
 M=sy.Matrix([[sy.Rational(v.numerator,v.denominator)for v in row]for row in B]);assert M==M.T and all(M[:k,:k].det()>0 for k in [1,2,3])
M=sy.Matrix(Bs[0]);Z=sy.Matrix(zz);assert sy.Rational(ss.numerator,ss.denominator)>=(Z.T*M.inv()*Z)[0]
W=[[sum(co[k]*Bs[k][i][j]for k in range(len(co)))/q for j in range(3)]for i in range(3)]
from flint import fmpq_poly,fmpq
hp=fmpq_poly([fmpq(19,6),-fmpq(13,102)]);powers=list(map(int,z['powers']));corr=sum((hp**powers[i]*fmpq(str(zz[i]))for i in range(3)),fmpq_poly([0]))+sum((hp**(powers[i]+powers[j])*fmpq(str(W[i][j]))for i in range(3)for j in range(3)),fmpq_poly([0]));corr*=fmpq(str(2/q));cor=list(map(lambda x:F(str(x)),corr));assert len(PP)==82

def ar(x):x=F(x);return arb(x.numerator)/x.denominator
def cap(R):
 if R<2:return 1
 ss=1-2/R;mu=F(1);best=None
 for k in range(501):
  rat=F(2*k+1,40+2*k);den=mu*(rat-ss*ss)
  if den>0:
   v=(1-ss*ss)/den
   if best is None or v<best:best=v
  mu*=rat
 return best.numerator//best.denominator

def shift(P,c):
 v=[ar(x)for x in P];c=ar(c)
 for i in range(len(v)-2,-1,-1):
  for j in range(i,len(v)-1):v[j]+=c*v[j+1]
 return v

def mul(A,B):
 C=[arb(0)for _ in range(len(A)+len(B)-1)]
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]+=a*b
 return C

def horner(P,x):
 v=arb(0)
 for a in reversed(P):v=v*x+a
 return v

rows=[];total=F(0);nodes=0
edges=sorted(set([F(17)+F(bi,4)for bi in range(653)]+[F(323,13)]))
for bi,(l,r)in enumerate(zip(edges[:-1],edges[1:])):
 center=(l+r)/2;rad=(r-l)/2;assert l>=17 and r<=180 and 0<r-l<=68;C=cap(r/17);assert C>=1;target=F(1,10**8)/((2**41)*C)
 terms=[(PP,F(1),F(0))]+([(cor,F(0),F(0))]if r<=F(323,13)else[])
 combined=[arb(0)]*213;err=arb(0)
 for P,lam,offset in terms:
  pol=shift(P,center);factor=(ar(offset-lam*center)).exp();pol=[v*factor for v in pol];ee=[ar((-lam)**j/F(math.factorial(j)))for j in range(21)];prod=mul(pol,ee)
  for j,v in enumerate(prod):combined[j]+=v
  ball=arb(0,ar(rad));pm=horner(pol,ball).abs_upper();rem=(ar(lam*rad)).exp()*ar((lam*rad)**21/F(math.factorial(21)));err+=pm*rem
 for j in range(33):
  vv=(horner(combined,ar(-rad+2*rad*F(j,32)))+arb(0,err.abs_upper())).upper()
  if vv>0:target=max(target,F(str(vv.fmpq()))*F(101,100))
 stack=[(-rad,rad)];mx=arb(0);count=0
 while stack:
  lo,hi=stack.pop();count+=1;nodes+=1;ball=arb(ar((lo+hi)/2),ar((hi-lo)/2));v=horner(combined,ball)+arb(0,err.abs_upper())
  if v.upper()<=ar(target) or hi-lo<F(1,10**5):
   if v.upper()>mx:mx=v.upper()
  else:
   if hi-lo<F(1,10**15):raise Exception(('failed',bi,float(lo),float(hi),str(v),str(err)))
   mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
 # Conservative fixed per-bin cap, independent of any decimal rounding.
 height=F(str(mx.fmpq()));total+=height*(2**41)*C;rows.append(dict(l=str(l),r=str(r),cap=C,height=str(height),nodes=count))
 if bi%10==0:print(bi,'nodes',nodes,'upper',str(mx),'err',str(err),flush=True)
json.dump(dict(status='ALL FINITE ANNULI PASS; tail and angular dependencies separately required',B=z['B'],co=z['co'],z=z['z'],s=str(ss),q=str(q),L=str(LL),P=raw['P'],budget=str(total),rows=rows,nodes=nodes),open(p+'annular_final_certificate.json','w'));print('PASS',nodes,'budget',total,flush=True)
