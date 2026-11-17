import sys,json,math,hashlib
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx,fmpq as Q,fmpq_poly as Poly
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=768;D=Path(__file__).parent;n=50;d=F(22);cut=F(180);R=F(19,13);support=d*R;raw=json.load(open(D/'repaired.json'));ker=json.load(open(D/'kernel.json'));assert ker['dimension']==n and ker['R']=='19/13' and ker['powers']==[12,24];assert raw['dimension']==n and F(raw['contact_x'])==d;PP=list(map(F,raw['P']));co=list(map(F,ker['co']));B=[[[F(v)for v in row]for row in b]for b in ker['B']];z=list(map(F,ker['z']));q=F(ker['q']);ss=F(ker['s']);assert q==sum(co)and co[0]==1 and min(co)>0
for b in B:assert b[0][0]>=0 and b[1][1]>=0 and b[0][0]*b[1][1]>=b[0][1]**2
b=B[0];det=b[0][0]*b[1][1]-b[0][1]**2;assert det>0 and ss>(b[1][1]*z[0]**2-2*b[0][1]*z[0]*z[1]+b[0][0]*z[1]**2)/det
W=[[sum(co[k]*B[k][i][j]for k in range(len(co)))/q for j in range(2)]for i in range(2)];qq=lambda v:Q(str(v));h=Poly([Q(19,6),-Q(13,6)/qq(d)])**12;kp=qq(2*z[0]/q)*h+qq((W[0][0]+2*z[1])/q)*h**2+qq(2*W[0][1]/q)*h**3+qq(W[1][1]/q)*h**4;KP=list(map(lambda v:F(str(v)),kp));assert kp(qq(support))==0

def ar(v):v=F(v);return arb(v.numerator)/v.denominator
def rat(v):return F(str(v.fmpq()))
def cap(radius):
 s=1-F(1,2)/radius;kmin=max(0,math.floor(float(s*(n-1)/(1-s)))+1);choices=[]
 for k in range(kmin,kmin+501):
  gap=k/(k+n-1)-float(s)
  if gap>0:choices.append((math.log((1-float(s))/gap)+math.lgamma(n-1+k)-math.lgamma(n-1)+math.lgamma((n-1)/2)-math.lgamma((n-1)/2+k),k))
 k=min(choices)[1];gap=F(k,k+n-1)-s;assert gap>0;ratio=F(1)
 for j in range(k):ratio*=F(n-1+j)/F(n-1+2*j,2)
 v=(1-s)*ratio/gap;return v.numerator//v.denominator,k

def shift(P,c):
 v=[ar(x)for x in P];c=ar(c)
 for i in range(len(v)-2,-1,-1):
  for j in range(i,len(v)-1):v[j]+=c*v[j+1]
 return v

def horner(P,x):
 v=arb(0)
 for c in reversed(P):v=v*x+c
 return v

def mul(A,B):
 C=[arb(0)for _ in range(len(A)+len(B)-1)]
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]+=a*b
 return C
edges=sorted(set([d,cut,support]+[d+F(j,4)for j in range(1,int(4*(cut-d))+1)]));assert edges[0]==d and edges[-1]==cut;rows=[];total=F(0);nodes=allleaves=0
for bi,(l,r)in enumerate(zip(edges,edges[1:])):
 assert 0<r-l<=d and not l<support<r;center=(l+r)/2;rad=(r-l)/2;C,kcap=cap(r/d);target=F(1,10**8)/C;pol=shift(PP,center);pol=[v*ar(-center).exp()for v in pol];combined=mul(pol,[ar(F((-1)**j,math.factorial(j)))for j in range(21)]);pm=horner(pol,arb(0,ar(rad))).abs_upper();err=pm*ar(rad).exp()*ar(rad**21/F(math.factorial(21)));E=rat(err.abs_upper())
 if r<=support:
  ks=shift(KP,center)
  for j,v in enumerate(ks):combined[j]+=v
 dp=[j*combined[j]for j in range(1,len(combined))];ddp=[j*dp[j]for j in range(1,len(dp))];stack=[(-rad,rad)];mx=F(0);bestlower=F(0);count=leaves=0
 while stack:
  lo,hi=stack.pop();count+=1;nodes+=1;mid=(lo+hi)/2;rr=(hi-lo)/2;ball=arb(ar(mid),ar(rr));natural=rat(horner(combined,ball).upper())+E;center_value=horner(combined,ar(mid));bestlower=max(bestlower,rat(center_value.lower())-E);C0=rat(center_value.upper());D0=rat(horner(dp,ar(mid)).abs_upper());M=rat(horner(ddp,ball).upper());assert D0>=0
  quadratic=C0-D0*D0/(2*M)if M<0 and D0< -M*rr else C0+D0*rr+M*rr*rr/2;upper=min(natural,quadratic+E)
  if upper<=bestlower+target or hi-lo<=F(1,2**18):leaves+=1;mx=max(mx,upper)
  else:mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
 assert not stack and count==2*leaves-1;allleaves+=leaves;total+=C*mx;rows.append(dict(l=str(l),r=str(r),cap=C,cap_degree=kcap,height=str(mx),nodes=count,leaves=leaves,kernel_active=bool(r<=support)))
 if bi%20==0:print('BIN',bi,'NODES',nodes,'BUDGET',float(total),flush=True)
assert nodes==2*allleaves-len(rows)
out=dict(status='FULL FINITE RADIAL COVERAGE PASS; kernel geometry and tail separate',dimension=n,d=str(d),cut=str(cut),support_cut=str(support),q=str(q),s=str(ss),P=raw['P'],kernel_polynomial=list(map(str,KP)),kernel_sha256=hashlib.sha256((D/'kernel.json').read_bytes()).hexdigest(),budget=str(total),rows=rows,nodes=nodes,leaves=allleaves,precision=ctx.prec,method='Uniform exponential Taylor remainder plus natural and concave-quadratic polynomial enclosures');(D/'annular_certificate.json').write_text(json.dumps(out,indent=2));print('PASS',len(rows),nodes,float(total),flush=True)
