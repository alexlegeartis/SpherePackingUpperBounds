import sys;sys.set_int_max_str_digits(0)
import json,math,time
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=768
p='research/lattice50_root_repair_20261117/'
z=json.load(open(p+'input_candidate.json'));raw=json.load(open(p+'repaired.json'));PP=list(map(F,raw['P']));a=F(str(z['Z'][1][1]));zz=F(str(z['Z'][0][1]))-F(1,10**9);s=zz*zz/a+F(1,10**12);q=F(1150750294032540977074575,1437673731494371);rho=F(1,10**70);LL=(1+rho)/(1-rho)
def ar(x):x=F(x);return arb(x.numerator)/x.denominator
def cap(R):
 if R<2:return 1
 ss=1-2/R;mu=F(1);best=None
 for k in range(501):
  rat=F(2*k+1,50+2*k);den=mu*(rat-ss*ss)
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
for bi in range(632):
 l=F(22)+F(bi,4);r=l+F(1,4);center=(l+r)/2;rad=(r-l)/2;C=cap(r/22);target=F(1,10**7)/((2**51)*C)
 terms=[(PP,F(1),F(0)),([F(0)]*96+[2*zz/q/F(22)**96],F(60,11),F(120)),([F(0)]*192+[2*a*LL/q/F(22)**192],F(120,11),F(240))]
 combined=[arb(0)]*213;err=arb(0)
 for P,lam,offset in terms:
  pol=shift(P,center);factor=(ar(offset-lam*center)).exp();pol=[v*factor for v in pol];ee=[ar((-lam)**j/F(math.factorial(j)))for j in range(21)];prod=mul(pol,ee)
  for j,v in enumerate(prod):combined[j]+=v
  ball=arb(0,ar(rad));pm=horner(pol,ball).abs_upper();rem=(ar(lam*rad)).exp()*ar((lam*rad)**21/F(math.factorial(21)));err+=pm*rem
 stack=[(-rad,rad)];mx=arb(0);count=0
 while stack:
  lo,hi=stack.pop();count+=1;nodes+=1;ball=arb(ar((lo+hi)/2),ar((hi-lo)/2));v=horner(combined,ball)+arb(0,err.abs_upper())
  if v.upper()<=ar(target):
   if v.upper()>mx:mx=v.upper()
  else:
   if hi-lo<F(1,10**15):raise Exception(('failed',bi,float(lo),float(hi),str(v),str(err)))
   mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
 # Conservative fixed per-bin cap, independent of any decimal rounding.
 total+=target*(2**51)*C;rows.append(dict(l=str(l),r=str(r),cap=C,height=str(target),nodes=count))
 if bi%10==0:print(bi,'nodes',nodes,'upper',str(mx),'err',str(err),flush=True)
json.dump(dict(status='ALL FINITE ANNULI PASS; tail and angular dependencies separately required',a=str(a),z=str(zz),s=str(s),q=str(q),L=str(LL),P=raw['P'],budget=str(total),rows=rows,nodes=nodes),open(p+'annular_certificate.json','w'));print('PASS',nodes,'budget',total,flush=True)
