import sys,json,math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=1024;n=51;D=Path(__file__).resolve().parent;raw=json.loads((D/'repaired.json').read_text());input=json.loads((D/'general_joint_candidate.json').read_text());contact=json.loads((D/'contact_cap_certificate.json').read_text());d=F(raw['contact_x']);cut=180 if n==51 else 250;PP=list(map(F,raw['P']));eps=F(contact['epsilon']);assert contact['dimension']==n and F(contact['s'])==1-F(1,2)/(1+eps)
if n==51:
 base=D;q=F(json.loads((base/'angular.json').read_text())['B_exact']);oldq=F(str(sum(json.loads((base/'angular_lp_candidate.json').read_text())['coefficients'])));a=F(str(input['Z'][1][1]))*q/oldq;z=F(str(input['Z'][0][1]))*q/oldq-F(1,10**8);s=z*z/a+F(1,10**12);assert a>0 and a+2*z<0 and s*a-z*z>0
else:a=z=s=F(0);q=F(1)
def ar(x):x=F(x);return arb(x.numerator)/x.denominator
def cap(R):
 ss=1-F(1,2)/R;kmin=max(0,math.floor(float(ss*(n-1)/(1-ss)))+1);choices=[]
 for k in range(kmin,kmin+501):
  gap=k/(k+n-1)-float(ss)
  if gap>0:choices.append((math.log((1-float(ss))/gap)+math.lgamma(n-1+k)-math.lgamma(n-1)+math.lgamma((n-1)/2)-math.lgamma((n-1)/2+k),k))
 k=min(choices)[1];gap=F(k,k+n-1)-ss;assert gap>0;ratio=F(1)
 for j in range(k):ratio*=F(n-1+j)/F(n-1+2*j,2)
 b=(1-ss)*ratio/gap;return b.numerator//b.denominator,k
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
edges=sorted(set([d,d*(1+eps)]+[d+F(j,4)for j in range(1,int(4*(cut-d))+1)]));assert edges[0]==d and edges[-1]==cut
rows=[];total=F(0);nodes=0
for bi,(l,r)in enumerate(zip(edges,edges[1:])):
 center=(l+r)/2;rad=(r-l)/2;assert 0<r-l<=d
 if bi==0:C=contact['cap'];kcap=None
 else:C,kcap=cap(r/d)
 target=F(1,10**7)/C;terms=[(PP,F(1),F(0))]
 if n==51:terms+=[([F(0)]*192+[2*z/q/d**192],F(216)/d,F(216)),([F(0)]*384+[a/q/d**384],F(432)/d,F(432))]
 combined=[arb(0)]*(max(len(p)for p,_,_ in terms)+20);err=arb(0)
 for P,lam,offset in terms:
  pol=shift(P,center);factor=ar(offset-lam*center).exp();pol=[v*factor for v in pol];ee=[ar((-lam)**j/F(math.factorial(j)))for j in range(21)];prod=mul(pol,ee)
  for j,v in enumerate(prod):combined[j]+=v
  pm=horner(pol,arb(0,ar(rad))).abs_upper();err+=pm*ar(lam*rad).exp()*ar((lam*rad)**21/F(math.factorial(21)))
 stack=[(-rad,rad)];mx=arb(0);count=leaves=0
 while stack:
  lo,hi=stack.pop();count+=1;nodes+=1;ball=arb(ar((lo+hi)/2),ar((hi-lo)/2));v=horner(combined,ball)+arb(0,err.abs_upper())
  if v.upper()<=ar(target) or hi-lo<=F(1,2**18):
   leaves+=1
   if v.upper()>mx:mx=v.upper()
  else:mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
 assert not stack and count==2*leaves-1;height=F(str(mx.fmpq()));total+=C*height;rows.append(dict(l=str(l),r=str(r),cap=C,cap_degree=kcap,height=str(height),nodes=count))
 if bi%10==0:print(n,bi,'nodes',nodes,'budget',float(total),flush=True)
out=dict(status='FINITE UNRESTRICTED ANNULI PASS; global dependencies separate',dimension=n,d=str(d),cut=cut,a=str(a),z=str(z),s=str(s),q=str(q),P=raw['P'],budget=str(total),rows=rows,nodes=nodes);(D/'annular_certificate.json').write_text(json.dumps(out,indent=2));print('PASS',n,nodes,float(total),flush=True)
