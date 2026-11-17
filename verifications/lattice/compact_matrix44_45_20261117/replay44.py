import sys,math,json,time,hashlib
from pathlib import Path
from fractions import Fraction as F
import sympy as s
from flint import arb,ctx
from jet import Jet,pe2
if not __debug__:raise RuntimeError('Do not use -O')
sys.set_int_max_str_digits(0);ctx.prec=384;D=Path(__file__).parent/'input44';mode=sys.argv[1];assert mode in ['main','cap'];z=json.load(open(D/'frozen.json'));co=list(map(s.Rational,z['co']));Bs=[[[s.Rational(x)for x in row]for row in B]for B in z['B']];x=s.symbols('x');pol=[]
for i,j in [(0,0),(0,1),(1,1)]:
 P=s.Poly(sum(co[k]*Bs[k][i][j]*s.gegenbauer(2*k,21,x)/s.gegenbauer(2*k,21,1)for k in range(len(co))),x);assert all(P.nth(k)==0 for k in range(1,P.degree()+1,2));pol.append([arb(int(P.nth(2*k).p))/int(P.nth(2*k).q)for k in range(len(co))])
assert len(Bs)==len(co) and co[0]==1 and min(co)>=0
for B in Bs:assert len(B)==2 and all(len(row)==2 for row in B) and B[0][1]==B[1][0] and B[0][0]>0 and B[1][1]>0 and B[0][0]*B[1][1]>B[0][1]**2

def ar(v):v=F(v);return arb(v.numerator)/v.denominator
def bil(t,u,v):return pe2(pol[0],t)+pe2(pol[1],t)*(u+v)+pe2(pol[2],t)*u*v
def parts(x,y,z):
 e=ar(F(6,13));R=ar(F(19,13));a=1+e*x;b=a+e*(1-x)*y;ha=(1-x)**12;hb=ha*(1-y)**12
 if mode=='cap':dot=(a+b-R)*z/2;return bil(dot**2/(a*b),ha,hb)
 c=b+e*(1-x)*(1-y)*z;hc=hb*(1-z)**12;dot=(a+b-c)/2
 return bil(dot**2/(a*b),ha,hb)+(1-z)**12*bil((a-dot)**2/(a*c),ha,hc)+((1-y)*(1-z))**12*bil((b-dot)**2/(b*c),hb,hc)
def check(box):
 mids=[(box[2*i]+box[2*i+1])/2 for i in range(3)];rads=[(box[2*i+1]-box[2*i])/2 for i in range(3)];args=[Jet(arb(ar(m),ar(r)),[arb(int(i==j))for j in range(3)])for i,(m,r)in enumerate(zip(mids,rads))];f=parts(*args);center=parts(*[Jet(ar(m),[arb(int(i==j))for j in range(3)])for i,m in enumerate(mids)])
 if center.v>0:raise RuntimeError(('COUNTEREXAMPLE',list(map(str,mids)),str(center.v)))
 first=center.v+sum((abs(f.g[i])*ar(rads[i])for i in range(3)),arb(0));second=center.v+sum((abs(center.g[i])*ar(rads[i])for i in range(3)),arb(0))+sum((abs(f.h[i][j])*ar(rads[i]*rads[j]/2)for i in range(3)for j in range(3)),arb(0));upper=min(f.v.upper(),first.upper(),second.upper())
 if upper<-ar(F(1,10**7)):return True,0,upper
 scores=[float((abs(center.g[i])*ar(rads[i])+sum((abs(f.h[i][j])*ar(rads[i]*rads[j])for j in range(3)),arb(0))).abs_upper())for i in range(3)];axis=max(range(3),key=lambda i:scores[i])if all(math.isfinite(v)for v in scores)else max(range(3),key=lambda i:rads[i]);return False,axis,upper

trace=(D/(mode+'_trace.bin')).read_bytes();stack=[tuple([F(0),F(1)]*3)];leaves=0;maxupper=-arb(10**20);start=time.time()
for num,code in enumerate(trace):
 assert stack and code in (0,1,2,3)
 box=stack.pop()
 if code:
  axis=code-1;i=2*axis;mid=(box[i]+box[i+1])/2;assert box[i]<mid<box[i+1]
  left=list(box);right=list(box);left[i+1]=mid;right[i]=mid;stack.extend([tuple(left),tuple(right)])
 else:
  ok,axis,upper=check(box);assert ok,(num,box,str(upper));leaves+=1;maxupper=max(maxupper,upper)
 if num%10000==0:print(mode,num,leaves,len(stack),time.time()-start,flush=True)
assert not stack and len(trace)==2*leaves-1
result=dict(status='PASS all accepted leaves independently replayed; exact complete trace coverage',mode=mode,nodes=len(trace),leaves=leaves,precision=ctx.prec,maxupper=str(maxupper),input_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),trace_sha=hashlib.sha256(trace).hexdigest())
(Path(__file__).parent/(mode+'_independent_result.json')).write_text(json.dumps(result,indent=2));print(result,flush=True)
