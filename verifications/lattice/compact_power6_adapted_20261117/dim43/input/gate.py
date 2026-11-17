import sys,math,json,time,hashlib
from pathlib import Path
from fractions import Fraction as F
import sympy as s
from flint import arb,ctx
from jet import Jet,pe2
if not __debug__:raise RuntimeError('Do not use -O')
sys.set_int_max_str_digits(0);ctx.prec=256;D=Path(__file__).parent;mode=sys.argv[1];OD=Path(sys.argv[sys.argv.index('--out')+1]) if '--out' in sys.argv else D;OD.mkdir(exist_ok=True);assert mode in ['main','cap'];z=json.load(open(D/'frozen.json'));co=list(map(s.Rational,z['co']));Bs=[[[s.Rational(x)for x in row]for row in B]for B in z['B']];x=s.symbols('x');pol=[]
for i,j in [(i,j)for i in range(3)for j in range(i,3)]:
 P=s.Poly(sum(co[k]*Bs[k][i][j]*s.gegenbauer(2*k,s.Rational(41,2),x)/s.gegenbauer(2*k,s.Rational(41,2),1)for k in range(len(co))),x);assert all(P.nth(k)==0 for k in range(1,P.degree()+1,2));pol.append([arb(int(P.nth(2*k).p))/int(P.nth(2*k).q)for k in range(len(co))])
assert len(Bs)==len(co) and co[0]==1 and min(co)>=0
for B in Bs:
 M=s.Matrix(B);assert M==M.T and all(M[:k,:k].det()>0 for k in [1,2,3])

def ar(v):v=F(v);return arb(v.numerator)/v.denominator
def bil(t,u,v):
 out=Jet(0);idx=0
 for i in range(3):
  for j in range(i,3):
   out+=pe2(pol[idx],t)*(u**i*v**j if i==j else u**i*v**j+u**j*v**i);idx+=1
 return out
def parts(x,y,z):
 e=ar(F(6,13));R=ar(F(19,13));a=1+e*x;b=a+e*(1-x)*y;ha=(1-x)**3;hb=ha*(1-y)**3
 if mode=='cap':dot=(a+b-R)*z/2;return bil(dot**2/(a*b),ha,hb)
 c=b+e*(1-x)*(1-y)*z;hc=hb*(1-z)**3;dot=(a+b-c)/2
 return bil(dot**2/(a*b),ha,hb)+(1-z)**6*bil((a-dot)**2/(a*c),ha,hc)+((1-y)*(1-z))**6*bil((b-dot)**2/(b*c),hb,hc)
def check(box):
 mids=[(box[2*i]+box[2*i+1])/2 for i in range(3)];rads=[(box[2*i+1]-box[2*i])/2 for i in range(3)];args=[Jet(arb(ar(m),ar(r)),[arb(int(i==j))for j in range(3)])for i,(m,r)in enumerate(zip(mids,rads))];f=parts(*args);center=parts(*[Jet(ar(m),[arb(int(i==j))for j in range(3)])for i,m in enumerate(mids)])
 if center.v>0:raise RuntimeError(('COUNTEREXAMPLE',list(map(str,mids)),str(center.v)))
 first=center.v+sum((abs(f.g[i])*ar(rads[i])for i in range(3)),arb(0));second=center.v+sum((abs(center.g[i])*ar(rads[i])for i in range(3)),arb(0))+sum((abs(f.h[i][j])*ar(rads[i]*rads[j]/2)for i in range(3)for j in range(3)),arb(0));upper=min(f.v.upper(),first.upper(),second.upper())
 if upper<-ar(F(1,10**7)):return True,0,upper
 scores=[float((abs(center.g[i])*ar(rads[i])+sum((abs(f.h[i][j])*ar(rads[i]*rads[j])for j in range(3)),arb(0))).abs_upper())for i in range(3)];axis=max(range(3),key=lambda i:scores[i])if all(math.isfinite(v)for v in scores)else max(range(3),key=lambda i:rads[i]);return False,axis,upper
stack=[tuple([F(0),F(1)]*3)];trace=bytearray();nodes=leaves=0;maximum=-arb(10**20);start=time.time()
if '--root' in sys.argv:stack=[tuple(map(F,json.loads(Path(sys.argv[sys.argv.index('--root')+1]).read_text())))];print('ROOT',list(map(str,stack[0])),flush=True)
if '--resume' in sys.argv:
 ck=json.load(open(OD/(mode+'_checkpoint.json')));trace=bytearray((OD/(mode+'_trace.bin')).read_bytes());assert len(trace)==ck['nodes'];nodes=ck['nodes'];leaves=ck['leaves'];stack=[tuple(map(F,b))for b in ck['pending']];maximum=-ar(F(1,10**7));print('RESUMING',mode,nodes,leaves,len(stack),flush=True)

while stack:
 box=stack.pop();ok,axis,upper=check(box);nodes+=1
 if ok:trace.append(0);leaves+=1;maximum=max(maximum,upper)
 else:
  i=2*axis;mid=(box[i]+box[i+1])/2;assert box[i+1]-box[i]>F(1,10**30);left=list(box);right=list(box);left[i+1]=mid;right[i]=mid;stack.extend([tuple(left),tuple(right)]);trace.append(axis+1)
 if nodes%5000==0:
  print(mode,nodes,len(stack),leaves,time.time()-start,str(upper),flush=True);(OD/(mode+'_trace.bin')).write_bytes(trace);(OD/(mode+'_checkpoint.json')).write_text(json.dumps(dict(nodes=nodes,leaves=leaves,pending=[list(map(str,b))for b in stack])))
 if nodes>=2000000:raise RuntimeError('Coverage incomplete: node budget reached')
(OD/(mode+'_trace.bin')).write_bytes(trace);out=dict(status='COMPLETE INTERVAL COVERAGE PASS',mode=mode,nodes=nodes,leaves=leaves,precision_bits=ctx.prec,threshold='-1/10000000',maximum_upper=str(maximum),input_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),trace_sha=hashlib.sha256(trace).hexdigest());(OD/(mode+'_result.json')).write_text(json.dumps(out,indent=2));print(out,flush=True)
