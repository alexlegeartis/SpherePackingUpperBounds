import sys,math,json,time,hashlib
from pathlib import Path
from fractions import Fraction as F
import sympy as s
from flint import arb,ctx
from jet_fast import Jet,pe2
if not __debug__:raise RuntimeError('Do not use -O')
sys.set_int_max_str_digits(0);ctx.prec=384;n=int(sys.argv[1]);D=Path(__file__).parent/f'dim{n}';mode='main';z=json.load(open(D/'frozen.json'));co=list(map(s.Rational,z['co']));Bs=[[[s.Rational(x)for x in row]for row in B]for B in z['B']];x=s.symbols('x');pol=[]
indices=[(i,j)for i in range(3)for j in range(i,3)]
for i,j in indices:
 P=s.Poly(sum(co[k]*Bs[k][i][j]*s.gegenbauer(k,s.Rational(n-2,2),x)/s.gegenbauer(k,s.Rational(n-2,2),1)for k in range(len(co))),x);pol.append([arb(int(P.nth(k).p))/int(P.nth(k).q)for k in range(len(co))])
assert len(Bs)==len(co) and co[0]==1 and min(co)>0
import itertools
for B in Bs:
 M=s.Matrix(B);assert M==M.T and M.shape==(3,3)
 for size in [1,2,3]:
  for ii in itertools.combinations(range(3),size):assert M.extract(ii,ii).det()>=0

def ar(v):v=F(v);return arb(v.numerator)/v.denominator
def bil(t,u,v):
 up=[u**i for i in range(3)];vp=[v**i for i in range(3)]
 return sum((pe2(pol[k],t)*(up[i]*vp[j] if i==j else up[i]*vp[j]+up[j]*vp[i])for k,(i,j)in enumerate(indices)),Jet(0))

assert z["R"]=="19/13" and z["powers"]==[6, 10, 14]
def parts(x,y,z):
 e=ar(F(6,13));a=1+e*x;b=a+e*(1-x)*y;c=b+e*(1-x)*(1-y)*z;ha=(1-x)**12;hb=ha*(1-y)**12;hc=hb*(1-z)**12
 ua=(1-x)**4;ub=ua*(1-y)**4;uc=ub*(1-z)**4
 return bil((a+b-c)/(2*(a*b).sqrt()),ua,ub)+(1-z)**6*bil((a+c-b)/(2*(a*c).sqrt()),ua,uc)+((1-y)*(1-z))**6*bil((b+c-a)/(2*(b*c).sqrt()),ub,uc)
def check(box):
 mids=[(box[2*i]+box[2*i+1])/2 for i in range(3)];rads=[(box[2*i+1]-box[2*i])/2 for i in range(3)];args=[Jet(arb(ar(m),ar(r)),[arb(int(i==j))for j in range(3)])for i,(m,r)in enumerate(zip(mids,rads))];f=parts(*args);center=parts(*[Jet(ar(m),[arb(int(i==j))for j in range(3)])for i,m in enumerate(mids)])
 if center.v>0:raise RuntimeError(('COUNTEREXAMPLE',list(map(str,mids)),str(center.v)))
 first=center.v+sum((abs(f.g[i])*ar(rads[i])for i in range(3)),arb(0));second=center.v+sum((abs(center.g[i])*ar(rads[i])for i in range(3)),arb(0))+sum((abs(f.h[i][j])*ar(rads[i]*rads[j]/2)for i in range(3)for j in range(3)),arb(0));upper=min(f.v.upper(),first.upper(),second.upper())
 if upper<0:return True,0,upper
 scores=[float((abs(center.g[i])*ar(rads[i])+sum((abs(f.h[i][j])*ar(rads[i]*rads[j])for j in range(3)),arb(0))).abs_upper())for i in range(3)];axis=max(range(3),key=lambda i:scores[i])if all(math.isfinite(v)for v in scores)else max(range(3),key=lambda i:rads[i]);return False,axis,upper

import itertools
reports=[]
for axis in range(3):
 for fixed in itertools.product([0,1],repeat=2):
  pairs=[(F(v),F(v))for v in fixed];pairs.insert(axis,(F(0),F(1)));stack=[tuple(sum((list(v)for v in pairs),[]))];nodes=leaves=0;maximum=-arb(10**20)
  while stack:
   box=stack.pop();ok,_,upper=check(box);nodes+=1
   if ok:leaves+=1;maximum=max(maximum,upper)
   else:
    i=2*axis;assert box[i+1]-box[i]>F(1,10**30);mid=(box[i]+box[i+1])/2;left=list(box);right=list(box);left[i+1]=mid;right[i]=mid;stack.extend([tuple(left),tuple(right)])
   assert nodes<100000
  assert nodes==2*leaves-1;reports.append(dict(axis=axis,fixed=fixed,nodes=nodes,leaves=leaves,maximum_upper=str(maximum)));print(reports[-1],flush=True)
(D/'edges_result.json').write_text(json.dumps(dict(status='COMPLETE12EDGES INTERVAL PASS',precision=ctx.prec,kernel_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),reports=reports),indent=2))
