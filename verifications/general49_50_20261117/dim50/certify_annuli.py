import sys;sys.set_int_max_str_digits(0)
import json,math,time
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=1024
if not __debug__: raise RuntimeError('Do not run this verifier with -O')
from pathlib import Path
p=str(Path(__file__).resolve().parent)+'/'
z=json.load(open(p+'general_joint_candidate.json'));raw=json.load(open(p+'repaired.json'));PP=list(map(F,raw['P']))
q=F(json.load(open(p+'angular.json'))['B_exact']);oldq=F(str(sum(json.load(open(p+'angular_lp_candidate.json'))['coefficients'])))
a=F(str(z['Z'][1][1]))*q/oldq;zz=F(str(z['Z'][0][1]))*q/oldq-F(1,10**8);s=zz*zz/a+F(1,10**12)
assert len(PP)==72 and a>0 and s*a-zz*zz>0 and a+2*zz<0
co=list(map(F,json.load(open(p+'angular.json'))['gegenbauer_coefficients']));assert co[0]==1 and min(co)>=0 and sum(co)==q

def ar(x):x=F(x);return arb(x.numerator)/x.denominator
def cap(R):
 ss=1-F(1,2)/R;n=50;kmin=max(0,math.floor(float(ss*(n-1)/(1-ss)))+1)
 choices=[]
 for k in range(kmin,kmin+501):
  gap=k/(k+n-1)-float(ss)
  if gap>0:choices.append((math.log((1-float(ss))/gap)+math.lgamma(n-1+k)-math.lgamma(n-1)+math.lgamma((n-1)/2)-math.lgamma((n-1)/2+k),k))
 k=min(choices)[1];gap=F(k,k+n-1)-ss;assert k>=0 and gap>0
 ratio=F(1)
 for j in range(k):ratio*=F(n-1+j)/F(n-1+2*j,2)
 bound=(1-ss)*ratio/gap
 return bound.numerator//bound.denominator,k

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
 l=F(22)+F(bi,4);r=l+F(1,4);center=(l+r)/2;rad=(r-l)/2;assert l>=22 and r<=180 and 0<r-l<=22;C,kcap=cap(r/22);assert C>=1;target=F(1,10**6)/C
 terms=[(PP,F(1),F(0)),([F(0)]*192+[2*zz/q/F(22)**192],F(108,11),F(216)),([F(0)]*384+[a/q/F(22)**384],F(216,11),F(432))]
 combined=[arb(0)]*405;err=arb(0)
 for P,lam,offset in terms:
  pol=shift(P,center);factor=(ar(offset-lam*center)).exp();pol=[v*factor for v in pol];ee=[ar((-lam)**j/F(math.factorial(j)))for j in range(21)];prod=mul(pol,ee)
  for j,v in enumerate(prod):combined[j]+=v
  ball=arb(0,ar(rad));pm=horner(pol,ball).abs_upper();rem=(ar(lam*rad)).exp()*ar((lam*rad)**21/F(math.factorial(21)));err+=pm*rem
 stack=[(-rad,rad)];mx=arb(0);count=0;leaves=0
 while stack:
  lo,hi=stack.pop();count+=1;nodes+=1;ball=arb(ar((lo+hi)/2),ar((hi-lo)/2));v=horner(combined,ball)+arb(0,err.abs_upper())
  if v.upper()<=ar(target) or hi-lo<=F(1,2**16):
   leaves+=1
   if v.upper()>mx:mx=v.upper()
  else:
   if hi-lo<F(1,10**15):raise Exception(('failed',bi,float(lo),float(hi),str(v),str(err)))
   mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
 assert not stack and count==2*leaves-1
 # Exact full binary coverage.
 height=F(str(mx.fmpq()));total+=height*C;rows.append(dict(l=str(l),r=str(r),cap=C,cap_degree=kcap,height=str(height),nodes=count))
 if bi%10==0:print(bi,'nodes',nodes,'upper',str(mx),'err',str(err),flush=True)
json.dump(dict(status='ALL FINITE ANNULI PASS; tail and angular dependencies separately required',a=str(a),z=str(zz),s=str(s),q=str(q),P=raw['P'],budget=str(total),rows=rows,nodes=nodes),open(p+'annular_final_certificate.json','w'));print('PASS',nodes,'budget',total,flush=True)
