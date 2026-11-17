import sys,math,json,time,hashlib
from pathlib import Path
from fractions import Fraction as F
import sympy as s
from flint import arb,ctx
from jet import Jet,pe2
if not __debug__:raise RuntimeError('Do not use -O')
sys.set_int_max_str_digits(0);ctx.prec=384;D=Path(__file__).parent/'input';mode=sys.argv[1];OD=Path(sys.argv[sys.argv.index('--out')+1]) if '--out' in sys.argv else D;OD.mkdir(exist_ok=True);assert mode in ['main','cap'];z=json.load(open(D/'frozen.json'));assert hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest()=='0619be767fbb2d7b15b157af8bd282126425442d1339ce24ec33f250aed61615';assert z['n']==43 and z['powers']==[6,9,12] and z['base_power']==6 and z['feature_stride']==3;co=list(map(s.Rational,z['co']));Bs=[[[s.Rational(x)for x in row]for row in B]for B in z['B']];x=s.symbols('x');pol=[]
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

tracepath=Path(sys.argv[sys.argv.index('--trace')+1])if '--trace'in sys.argv else D/(mode+'_trace.bin')
trace=tracepath.read_bytes();stack=[tuple([F(0),F(1)]*3)];leaves=0;maxupper=-arb(10**20);start=time.time();offset=0
outroot=Path(__file__).parent
if '--resume-audit'in sys.argv:
 ck=json.loads((outroot/'prefix_audit_checkpoint.json').read_text());offset=ck['nodes'];assert ck['precision']==ctx.prec and hashlib.sha256(trace[:offset]).hexdigest()==ck['prefix_sha'] and ck['input_sha']==hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest();stack=[tuple(map(F,b))for b in ck['pending']];leaves=ck['leaves'];maxupper=-ar(F(1,10**7))
for num,code in enumerate(trace[offset:],offset):
 assert stack and code in (0,1,2,3)
 box=stack.pop()
 if code:
  axis=code-1;i=2*axis;mid=(box[i]+box[i+1])/2;assert box[i]<mid<box[i+1]
  left=list(box);right=list(box);left[i+1]=mid;right[i]=mid;stack.extend([tuple(left),tuple(right)])
 else:
  ok,axis,upper=check(box);assert ok,(num,box,str(upper));leaves+=1;maxupper=max(maxupper,upper)
 if num%10000==0:print(mode,num,leaves,len(stack),time.time()-start,flush=True)
if '--partial'in sys.argv:
 ck=dict(status='PARTIAL ONLY; NOT A CERTIFICATE',nodes=len(trace),leaves=leaves,pending=[list(map(str,b))for b in stack],prefix_sha=hashlib.sha256(trace).hexdigest(),input_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),precision=ctx.prec)
 (outroot/'prefix_audit_checkpoint.json').write_text(json.dumps(ck,indent=2));print(ck['status'],len(trace),leaves,len(stack));sys.exit(0)
assert not stack and len(trace)==2*leaves-1
result=dict(status='PASS all accepted leaves independently replayed; exact complete trace coverage',mode=mode,nodes=len(trace),leaves=leaves,precision=ctx.prec,maxupper=str(maxupper),maxupper_note=('conservative -1e-7 ceiling after prefix resume; not measured worst leaf' if offset else 'maximum of replayed leaf upper bounds'),input_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),trace_sha=hashlib.sha256(trace).hexdigest())
(outroot/(mode+'_independent_result.json')).write_text(json.dumps(result,indent=2));print(result,flush=True)
