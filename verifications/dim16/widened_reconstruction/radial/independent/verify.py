"""Independent exact Fourier/export, ray and full piecewise pair audit of fixed16 witness.
No producer imports. Closed-form Laguerre and direct Fourier identities are independent.
"""
if not __debug__:raise RuntimeError('Run without optimization')
import json,sys,time,hashlib,math
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
sys.set_int_max_str_digits(0);ctx.prec=384
D=Path(__file__).resolve().parent;start=time.monotonic()
def read(f):return json.loads((D/f).read_text())
def aq(x):
 x=F(x);return arb(Q(x.numerator,x.denominator))
def poly(xs):return Poly([Q(str(F(x))) for x in xs])
def out(name,x):(D/name).write_text(json.dumps(x,indent=2)+'\n')
cert=read('candidate_exact.json');source=read('conditional_refit_exact.json');raw=read('conditional_refit.json')
assert hashlib.sha256((D/'candidate_exact.json').read_bytes()).hexdigest()=='8e47f1948dde8765ee8849c3c73dbb39eb55fd49b3ca065a1677ab04e8cb3be8'
P=poly(cert['physical']);H=poly(cert['fourier']);C=F(cert['C_exact']);lam=F(cert['lambda_count_exact'])
assert C>=0 and lam>=0 and H(0)==1 and cert['U_exact']=='7724'
assert F(cert['threshold_squared_exact'])==F(1000,499)
assert len(raw['radial'])==64 and P.degree()==H.degree()==63
# Closed form, not the exporter's three-term recurrence.
Ps=Poly();Hs=Poly()
for k,val in enumerate(raw['radial']):
 co=F(float(val));norm=math.comb(k+7,k)
 L=poly([F((-2)**j*math.comb(k+7,k-j),math.factorial(j)*norm)for j in range(k+1)])
 Ps+=Q(str(co))*L;Hs+=((-1)**k)*Q(str(co))*L
assert Ps==poly(source['physical']) and Hs==poly(source['fourier'])
assert F(source['C'])==F(float(raw['C'])) and F(source['lambda_count_exact'])==F(float(raw['lambda_count']))
N=Hs(0);eps=Q(1,10**6);E=Poly([1,Q(-1,4)])
assert N>0 and F(cert['epsilon_exact'])==F(1,10**6)
assert P==(Ps/N+eps*E)/(1-eps) and H==(Hs/N-eps*E)/(1-eps)
assert Q(str(C))==Q(source['C'])/N/(1-eps)
assert Q(str(lam))==Q(source['lambda_count_exact'])/N/(1-eps)
# Direct monomial Gaussian Fourier map: X^k -> k! L_k^7(Y).
Hdirect=Poly()
for k,co in enumerate(P.coeffs()):
 L=poly([F((-1)**j*math.factorial(k)*math.comb(k+7,k-j),math.factorial(j))for j in range(k+1)])
 Hdirect+=co*L
assert Hdirect==H
algebra={'status':'PASS','degree':63,'source_laguerre_closed_form':True,'direct_monomial_fourier_identity':True,'normalization':'1','C_exact':str(C),'lambda_exact':str(lam),'exact_repair':'epsilon=1e-6, degree1 normalized Laguerre, renormalized by1-epsilon'}
out('ALGEBRA.json',algebra);print('ALGEBRA PASS',flush=True)

# Exact rational Sturm with positive primitive rescaling. No complex root backend.
def primitive(p):
 z=p.numer();g=z.content();assert g>0;return Poly(z//g)
def variations(vals):
 z=[1 if x>0 else -1 for x in vals if x];return sum(x!=y for x,y in zip(z,z[1:]))
ray_results={}
for label,p,lo,sign in [('fourier',H,Q(0),1),('physical',P,Q(21,2),-1)]:
 chain=[primitive(p),primitive(p.derivative())]
 while chain[-1].degree()>0:
  rem=-(chain[-2]%chain[-1]);assert rem;chain.append(primitive(rem))
  if len(chain)%10==0:print(label,'STURM',len(chain),time.monotonic()-start,flush=True)
 v0=variations([r(lo)for r in chain]);vinf=variations([r[r.degree()]for r in chain])
 assert v0==vinf and sign*p(lo)>0
 ray_results[label]={'status':'PASS','ray_start':str(lo),'sign':sign,'degrees':[r.degree()for r in chain],'endpoint_variations':v0,'positive_infinity_variations':vinf}
 out(label+'_STURM.json',{'polynomials':[[str(x)for x in r.coeffs()]for r in chain],**ray_results[label]})
 print(label,'FULL RAY PASS',flush=True)
out('RAYS.json',ray_results)

# Independent pair proof: first-derivative mean-value, exact piecewise boundaries.
fit=json.loads((D/'inputs/fit.json').read_text(),parse_float=str);ang=read('inputs/angular_dim16.json');pars=list(map(F,fit['parameters']))
VK=list(map(F,['1.02','1.05','1.1','1.2','1.4','2']));LK=list(map(F,['1.2','1.4','1.5','1.6','1.65','1.75']))
assert len(pars)==14 and min(pars)>=0
assert hashlib.sha256((D/'inputs/fit.json').read_bytes()).hexdigest()=='438c5f0a7d55eca26e6b4c5adbba0c3ec542f0bb18cda15314596456ff2cd4aa'
assert hashlib.sha256((D/'inputs/angular_dim16.json').read_bytes()).hexdigest()=='bbe0a3d60b7f808c21c81631eb584bf665f38402110639eba0250bae0c4adee4'
B=F(ang['B_exact']);assert B>0
TP=2*arb.pi();CR=aq(C/B);AL=aq(lam);AP=arb_poly([arb(x)for x in P.coeffs()]);DP=AP.derivative()
threshold=F(500,499);end=F(7,4)
def parts(mid):
 a=pars[0];b=pars[7];c=-aq(pars[0])
 for w,k in zip(pars[1:7],VK):
  if mid>k:a+=w;c-=aq(w*k)
 for w,k in zip(pars[8:],LK):
  if mid>k:b+=w;c-=aq(w)*aq(k).log()
 return a,b,c

def phival(v,part):
 a,b,c=part;return aq(a)*v+aq(b)*v.log()+c

def profile_interval(lo,hi,part):return phival(aq(lo),part).union(phival(aq(hi),part))

def centered(poly,lo,hi):
 mid=(lo+hi)/2;rad=(hi-lo)/2;cs=poly(arb_poly([mid,rad]));return cs[0]+arb(0,sum((abs(x)for x in cs.coeffs()[1:]),arb(0)).upper())

def eval_point(v,part,active):
 x=TP*aq(v);h=(-phival(aq(v),part)).exp()
 return AP(x)-CR*x.exp()*h*(2-h)-(AL*x.exp()if active else arb(0))

def pair_upper(lo,hi,part,active):
 v=aq(lo).union(aq(hi));x=TP*v;h=(-profile_interval(lo,hi,part)).exp()
 a,b,_=part;dphi=aq(a)+aq(b)/v
 # d/dv [e^X(2h-h^2)]=e^X[2pi*h*(2-h)-2phi'*h*(1-h)]
 # This factorization differs from the producer's Taylor-2 formula.
 deriv=TP*centered(DP,TP*aq(lo),TP*aq(hi))-CR*x.exp()*(TP*h*(2-h)-2*dphi*h*(1-h))
 if active:deriv-=AL*TP*x.exp()
 return eval_point((lo+hi)/2,part,active)+abs(deriv)*aq((hi-lo)/2)

# Validate original partition with exact endpoints and every knot, then refine separately.
saved=read('partition.json');ordered=sorted(saved,key=lambda z:F(z['lo']))
assert F(ordered[0]['lo'])==1 and F(ordered[-1]['hi'])==end
for i,z in enumerate(ordered):
 lo,hi=F(z['lo']),F(z['hi']);assert lo<hi
 if i:assert F(ordered[i-1]['hi'])==lo
 assert not(lo<threshold<hi)
 assert bool(z['count_active'])==((lo+hi)/2<threshold)
 assert all(not(lo<k<hi)for k in VK+LK)
assert any(F(z['lo'])==threshold and not z['count_active']for z in ordered)
stack=[(F(z['lo']),F(z['hi']),0)for z in ordered];leaves=[];checks=0
while stack:
 lo,hi,depth=stack.pop();checks+=1;assert checks<200000 and depth<24
 part=parts((lo+hi)/2);active=(lo+hi)/2<threshold;upper=pair_upper(lo,hi,part,active)
 if upper<0:leaves.append({'lo':str(lo),'hi':str(hi),'active':active,'upper':str(upper),'refinement':depth});continue
 mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
leaves.sort(key=lambda x:F(x['lo']));assert all(F(a['hi'])==F(b['lo'])for a,b in zip(leaves,leaves[1:]));assert sum((F(z['hi'])-F(z['lo'])for z in leaves),F(0))==F(3,4)
assert TP*aq(end)>aq(F(21,2))
# Profile nonnegative and increasing => 0<h<=1. Beyond7/4 physical ray is negative,
# while kernel and count allowances remain nonnegative, including later kink2.
pair={'status':'PASS','method':'Independent first-derivative mean-value, centered polynomial derivative range, endpoint-monotone profile range','checks':checks,'leaves':len(leaves),'max_refinement':max(z['refinement']for z in leaves),'coverage':['1','7/4'],'total_width':'3/4','count_threshold':str(threshold),'right_limit_without_count':True,'tail':'X>=7pi/2>21/2; P<0 and both allowances>=0'}
out('PAIR.json',pair);out('partition_independent.json',leaves);print('FULL PAIR PASS',pair,flush=True)

# Rational Machin pi bounds and the exact volume/objective. No Arb pi in final rounding.
def atan_bounds(q,terms=160):
 s=F(0)
 for k in range(terms):s+=F((-1)**k,(2*k+1)*q**(2*k+1))
 nxt=F((-1)**terms,(2*terms+1)*q**(2*terms+1));return min(s,s+nxt),max(s,s+nxt)
l5,u5=atan_bounds(5);l239,u239=atan_bounds(239);pi_lo=16*l5-4*u239;pi_hi=16*u5-4*l239
assert F(31415,10000)<pi_lo<pi_hi<F(31416,10000)
objective=F(str(P(0)))+C+7724*lam;obj7725=objective+lam
vl=pi_lo**8/(2**8*math.factorial(8));vu=pi_hi**8/(2**8*math.factorial(8))
dlow=vl*objective;dhigh=vu*objective;target=F('0.024900975');comparator=F('0.02492121')
assert 0<dlow<dhigh<target<comparator
rounding={'status':'PASS','objective_exact':str(objective),'count':7724,'volume_factor':'pi^8/(256*8!)','target_exact':str(target),'density_interval_arb_for_display':str(aq(objective)*arb.pi()**8/(2**8*math.factorial(8))),'pi_machin_terms':160,'pi_lower_exact':str(pi_lo),'pi_upper_exact':str(pi_hi),'density_lower_exact':str(dlow),'density_upper_exact':str(dhigh),'published_comparator':'0.02492121','conservative7725_display':str(aq(obj7725)*arb.pi()**8/(2**8*math.factorial(8)))}
out('ROUNDING.json',rounding)
result={'status':'PASS','scope':'Complete independent frozen radial witness audit; density additionally uses independently audited angular geometry and spherical count','algebra':algebra,'rays':ray_results,'pair':pair,'strict_density_upper':'0.024900975','seconds':time.monotonic()-start,'input_sha256':{f:hashlib.sha256((D/f).read_bytes()).hexdigest()for f in ['candidate_exact.json','conditional_refit.json','conditional_refit_exact.json','partition.json','inputs/fit.json','inputs/angular_dim16.json']}}
out('VERIFIED.json',result);print('FULL RADIAL PASS',result,flush=True)
