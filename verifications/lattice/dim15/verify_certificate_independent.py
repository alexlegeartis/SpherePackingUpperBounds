if not __debug__:raise RuntimeError('Assertions required; do not run with python -O')
import json,math,time,hashlib
from fractions import Fraction as F
from pathlib import Path
src=Path(__file__).with_name('certificate.json')
c=json.loads(src.read_text());n=c['n'];tau=F(c['tau']);gap=F(c['gap']);d=list(map(F,c['d']));E=list(map(F,c['E']));Gam=list(map(F,c['Gam']));J=E.index(F(2));M=len(E)-1
assert n==15 and tau==17 and gap==F(1,5) and c['exp_terms']==128
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
Ls=[[F(1)],[F(15,2),F(-1)]]
for k in range(2,len(d)):
 Ls.append(scale(add(mul([(2*k+F(11,2)),F(-1)],Ls[-1]),scale(Ls[-2],-(k+F(11,2)))),F(1,k)))
p=[F(0)];q=[F(0)]
for k,coef in enumerate(d):
 cc=coef/math.prod((F(13,2)+j)/j for j in range(1,k+1));p=add(p,scale(Ls[k],cc));q=add(q,scale(Ls[k],(-1)**k*cc))
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
nn=M-J;ad=[E[J+i+1]-2 for i in range(nn)]
provided={}
for r in c['profile_duals']:
 key=(r['j'],r['h']);assert key not in provided;provided[key]=r
expected=set();budgets=[];affected=[]
for j in range(-1,J):
 caps=[B(E[J+b+1]) if j<0 else A(E[j+1],E[J+b+1]) for b in range(nn)]
 for b in range(nn-2,-1,-1):caps[b]=min(caps[b],caps[b+1])
 t=n if j<0 else n-1
 shortdev=F(0) if j<0 else E[j+1]-2
 shortheight=F(0) if j<0 else heights[j]
 def increments(limit):
  cumul=[min(x,limit) for x in caps]
  return [cumul[0]]+[b-a for a,b in zip(cumul,cumul[1:])]
 mintrace=shortdev+sum(a*b for a,b in zip(ad,increments(t)))
 if caps[-1]<t or mintrace>=gap:
  budget=shortheight+sum(a*b for a,b in zip(heights[J:],increments(caps[-1])))
 else:
  affected.append(j)
  budget=shortheight+sum(a*b for a,b in zip(heights[J:],increments(t-1)))
  for h in range(nn):
   trace_rhs=t*ad[h]+shortdev-gap
   if caps[h]<t or trace_rhs<0:continue
   key=(j,h);expected.add(key);r=provided[key];y=list(map(F,r['y']))
   assert len(y)==nn+3 and min(y)>=0
   lhs=[F(0)]*nn;tail=F(0)
   for b in range(nn-1,-1,-1):tail+=y[b];lhs[b]=tail
   for b in range(nn):
    if b<h:lhs[b]+=y[nn]+(ad[h]-ad[b])*y[nn+2]
    if b<=h:lhs[b]-=y[nn+1]
    assert lhs[b]>=heights[J+b]
   upper=shortheight+sum(v*yy for v,yy in zip(caps,y))+t*y[nn]-t*y[nn+1]+trace_rhs*y[nn+2]
   assert upper==F(r['bound']);budget=max(budget,upper)
 budgets.append(budget)
assert expected==set(provided) and len(expected)==46
assert affected==c['affected'] and budgets==list(map(F,c['budgets']))
kap=max(budgets);assert kap==F(c['kappa'])
ratio=(p[0]+2*(2**15-1)*kap)/q[0];assert ratio>0
D=math.prod(range(1,16,2))
squared_pi=tau**15*F(2)**(1-30)/D**2*ratio**2
assert squared_pi==F(c['density_squared_times_pi'])
# Independently certify producer's coarse rational pi lower bound.
def at(q,N):return sum(((-1)**j*F(1,(2*j+1)*q**(2*j+1)) for j in range(N)),F(0))
pilower=F(c['pi_lower']);assert 0<pilower<16*at(5,20)-4*at(239,21)
bound=F(c['density']);assert bound>0 and bound**2*pilower>=squared_pi
assert bound<=F('0.025511135431')
# I<=4 low branch: Delta=4*pi^5/(sqrt2*11!!).
assert 16**2*F(22,7)**14/(2*D**2)<bound**2
print('ALL46 PROFILEDUALS AND OBJECTIVE PASS',float(bound),flush=True)
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
result={'passed':True,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'exact_bound':str(bound),'decimal_upper':'0.025511135431','duals_checked':len(expected),'interval_count':M,'interval_depth':mx,'global_signs':halfproofs,'low_branch':'rho<=1024sqrt2','trace_gap':'gap1/5; independently audited in audit15_theorem.py'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2))
print('FULL INDEPENDENT15D CERTIFICATE AUDIT PASS',flush=True)
