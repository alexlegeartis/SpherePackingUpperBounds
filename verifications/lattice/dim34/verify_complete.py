import json,math,time,hashlib,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
if not __debug__:raise RuntimeError('Run without -O; assertions are required')
from fractions import Fraction as F
from pathlib import Path
src=Path(__file__).with_name('certificate.json')
c=json.loads(src.read_text());n=c['n'];tau=F(c['tau']);gap=F(c['gap']);d=list(map(F,c['d']));E=list(map(F,c['E']));Gam=list(map(F,c['Gam']));J=E.index(F(2));M=len(E)-1
assert n==34 and tau==36 and gap==0 and c['exp_terms']==128
assert E[0]==1 and E[-1]==F(9,4) and all(a<b for a,b in zip(E,E[1:]))
assert len(Gam)==M==120 and min(Gam)>=0
def add(a,b):
 r=[F(0)]*max(len(a),len(b))
 for i,v in enumerate(a):r[i]+=v
 for i,v in enumerate(b):r[i]+=v
 return r

def mul(a,b):
 r=[F(0)]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):r[i+j]+=u*v
 return r

def scale(a,v):return [x*v for x in a]
# Independent recurrence for ordinary generalized Laguerre polynomials alpha7.
Ls=[[F(1)],[F(34,2),F(-1)]]
for k in range(2,len(d)):
 Ls.append(scale(add(mul([(2*k+F(30,2)),F(-1)],Ls[-1]),scale(Ls[-2],-(k+F(30,2)))),F(1,k)))
p=[F(0)];q=[F(0)]
for k,coef in enumerate(d):
 cc=coef/math.prod((F(32,2)+j)/j for j in range(1,k+1));p=add(p,scale(Ls[k],cc));q=add(q,scale(Ls[k],(-1)**k*cc))
assert p==list(map(F,c['p'])) and q==list(map(F,c['q']))
print('independent Laguerre recurrence matches p,q',flush=True)

def exp_low(x,N=128):return sum((x**j/F(math.factorial(j)) for j in range(N+1)),F(0))

def transform(poly,a,w):
 # p(a+wz) via independent direct binomial expansion.
 return [sum((poly[k]*F(math.comb(k,j))*a**(k-j)*w**j for k in range(j,len(poly))),F(0)) for j in range(len(poly))]

def bern01(poly):
 deg=len(poly)-1
 return [sum((poly[j]*F(math.comb(i,j),math.comb(deg,j)) for j in range(i+1)),F(0)) for i in range(deg+1)]

def prove_nonnegative(b,depth=0):
 if min(b)>=0:return depth
 assert b[0]>=0 and b[-1]>=0 and depth<20
 levels=[b]
 while len(levels[-1])>1:levels.append([(x+y)/2 for x,y in zip(levels[-1],levels[-1][1:])])
 return max(prove_nonnegative([r[0] for r in levels],depth+1),prove_nonnegative([r[-1] for r in levels[::-1]],depth+1))

mx=0
for b in range(M):
 a=tau*E[b];w=tau*(E[b+1]-E[b]);pt=transform(p,a,w)
 envelope=[Gam[b],Gam[b]*w/2,Gam[b]*w*w/8]
 mx=max(mx,prove_nonnegative(bern01(add(envelope,scale(pt,-1)))))
print('ALL120 INTERVALS PASS',mx,flush=True)
def B(S):
 z=1-2/S;assert 1-n*z*z>0
 return (n*(1-z*z))//(1-n*z*z)
def A(a,S):
 if a+S<=4:return 0
 z=1-2/S;v=(a+S-4)**2/(4*a*S);z=(z+v)/(1-v)
 return min(B(S)-1,((n-1)*(1-z*z))//(1-(n-1)*z*z)) if (n-1)*z*z<1 else B(S)-1
heights=[v/exp_low(tau*e/2) for e,v in zip(E,Gam)]
for b in range(M-2,J-1,-1):heights[b]=max(heights[b],heights[b+1])
from gegenbauer import gegen
# Independently reconstruct projective polynomials and prove cap inequalities.
nu=list(map(F,c['nu']));assert len(nu)==len(c['angular_cuts']) and min(nu)>=0
cutweights=[];constants=[];cap_options=set()
for cut in c['angular_cuts']:
 U=E[1];assert list(map(F,cut['source_interval']))==[1,U] and cut['n']==n
 co=list(map(F,cut['gegenbauer_even']));assert co[0]==1 and min(co)>=0
 f=[F(0)]
 for k,v in enumerate(co):f=add(f,scale(gegen(n,2*k)[::2],v))
 while len(f)>1 and f[-1]==0:f.pop()
 assert f==list(map(F,cut['f_y']))
 prove_nonnegative(bern01(transform(scale(f,-1),0,(U-1)**2)))
 f1=sum(f);N0=F(cut['constant']);assert f1==F(cut['f1']) and N0>0
 rw=[]
 for j in range(J):
  a,b=E[j:j+2];lo=max(F(0),(2-b)/(2*U));hi=min(1-1/(2*U),(2*U-a)/2)
  cap=cut['caps'][j];assert list(map(F,cap['cosine_squared_interval']))==[lo*lo,hi*hi]
  mj=F(cap['M']);assert mj>=0
  # Direct adaptive Bernstein of M-f, independently of producer's fixed partition.
  pp=scale(f,-1);pp[0]+=mj
  prove_nonnegative(bern01(transform(pp,lo*lo,hi*hi-lo*lo)))
  r=U-a/4;z=1-1/(2*r);eta=(U-1)**2/(a*(4-a));z=(z+eta)/(1-eta);dim=n-1
  opts=[]
  for power in (1,2,3):
   # Verify positive Gegenbauer expansion of x^(2p)-z²*x^(2p-2).
   poly=[F(0)]*(power+1);poly[power]=1;poly[power-1]=-z*z
   remainder=poly[:];coeff=[F(0)]*(power+1)
   for k in range(power,-1,-1):
    gk=gegen(dim,2*k)[::2];coeff[k]=remainder[k]/gk[k]
    for h,v in enumerate(gk):remainder[h]-=coeff[k]*v
   assert all(v==0 for v in remainder)
   if coeff[0]>0 and min(coeff[1:])>=0:
    opts.append((1-z*z)//coeff[0])
  assert opts;C=2*min(opts);assert C==cap['C']
  rw.append((2-f1/N0 if j==0 else 0)-C*mj/N0)
 assert rw==list(map(F,cut['weights']));cutweights.append(rw);constants.append(N0)
print('ALL ANGULAR POLYNOMIALS, CONTINUOUS DOMAINS, AND CODE CAPS PASS',flush=True)
basebudgets=[];budgets=[]
for j in range(-1,J):
 caps=[B(E[h+1]) if j==-1 else A(E[j+1],E[h+1]) for h in range(J,M)]
 for h in range(len(caps)-2,-1,-1):caps[h]=min(caps[h],caps[h+1])
 value=F(0) if j==-1 else heights[j]
 for h,cap in enumerate(caps):value+=(cap-(caps[h-1] if h else 0))*heights[J+h]
 basebudgets.append(value)
 if j>=0:value-=sum(v*w[j] for v,w in zip(nu,cutweights))
 budgets.append(value)
assert basebudgets==list(map(F,c['basebudgets'])) and budgets==list(map(F,c['budgets']))
kappa=max(budgets);assert kappa==F(c['kappa']);extra=2*sum(v*b for v,b in zip(nu,constants));assert extra==F(c['extra_objective'])
ratio=(p[0]+2*(2**n-1)*kappa+extra)/q[0];assert ratio>0
density=(tau/8)**(n//2)/math.factorial(n//2)*ratio
assert density==F(c['density'])
bound=F('0.000053148401');assert 0<density<=bound
print('ALL97 PROFILE BRANCHES AND EVEN NORMALIZATION PASS',float(density),flush=True)
halfproofs={}
for name,pol,left in [('fourier',q,F(0)),('tail',[-z for z in p],tau*E[-1])]:
 deg=len(pol)-1
 shifted=transform(pol,left,F(1))
 bb=[v/F(math.comb(deg,i)) for i,v in enumerate(shifted)]
 stack=[(bb,0)];leaves=0;depthmax=0
 while stack:
  bb,depth=stack.pop();depthmax=max(depthmax,depth)
  if min(bb)>0:leaves+=1;continue
  assert bb[0]>0 and bb[-1]>0 and depth<80
  ls=[bb[0]];rs=[bb[-1]]
  for _ in range(deg):bb=[(u+v)/2 for u,v in zip(bb,bb[1:])];ls.append(bb[0]);rs.append(bb[-1])
  stack.extend([(ls,depth+1),(rs[::-1],depth+1)])
 halfproofs[name]={'positive':True,'leaves':leaves,'depth':depthmax}
 print('GLOBAL',name,halfproofs[name],flush=True)

result={'passed':True,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'bound':str(bound),'intervals':M,'branches':J+1,'angular_cuts':len(nu),'global_signs':halfproofs}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2))
print('FULL INDEPENDENT34D SCHUR AUDIT PASS',flush=True)
