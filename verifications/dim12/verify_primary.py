"""Complete exact 12D curved-profile cyclic certificate replay."""
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
import json,math,subprocess,sys,hashlib,time
if not __debug__:raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0)
ctx.prec=384;D=Path(__file__).resolve().parent;start=time.time()
r=json.loads((D/'certificate.json').read_text());ang=json.loads((D/'angular_exact.json').read_text())
assert r['dimension']==ang['dimension']==12 and r['degree']==63
n=12;s=Q(r['s']);aa=Q(r['a_exact']);bb=Q(r['b_exact']);cut=Q(r['cut_X']);B=Q(r['B_exact']);M=Q(r['M_exact'])
assert s==Q(12,8) and aa==10 and bb==16 and cut==12 and M>0
assert B==Q(ang['B_exact']) and aa==Q(ang['a_exact']) and bb==Q(ang['b_exact'])
pol=P(list(map(Q,r['physical'])));H=P(list(map(Q,r['fourier'])));assert pol.degree()==H.degree()==63
assert Q(r['origin_cost_exact'])==M*B and Q(r['normalization_exact'])==H(0)>0
# Deterministic export from the final saved numerical discovery checkpoint.
j=json.loads((D/'discovery/radial_cut12.json').read_text());assert j['dimension']==n and j['degree']==63 and Q(j['s_exact'])==s and Q(j['a_exact'])==aa and Q(j['b_exact'])==bb and Q(str(F(j['cut_X'])))==cut
v=[Q(str(F(x).limit_denominator(10**16)))for x in j['radial']];C=Q(str(F(j['C']).limit_denominator(10**16)))
X=P([0,1]);alpha=Q(n,2);E=[P([1]),alpha-2*X]
for k in range(1,63):E.append(((2*k+alpha-2*X)*E[-1]-(k+alpha-1)*E[-2])/(k+1))
E=[p/p(0)for p in E]
assert pol==sum((x*p for x,p in zip(v,E)),P([]))
assert H==sum(((-1)**k*x*p for k,(x,p)in enumerate(zip(v,E))),P([]))
assert M==C/B*Q(1000001,1000000)
# Second, monomial Fourier formula: F[X^k exp(-X)]=k! L_k^(n/2-1)(Y) exp(-Y).
L=[P([1]),alpha-X]
for k in range(1,63):L.append(((2*k+alpha-X)*L[-1]-(k+alpha-1)*L[-2])/(k+1))
assert H==sum((pol[k]*math.factorial(k)*L[k]for k in range(64)),P([]))
# Angular verifier has exact PSD/sign-transition tests and complete finite/tail geometry.
old_hash=hashlib.sha256((D/'angular_exact.json').read_bytes()).hexdigest()
subprocess.run([sys.executable,str(D/'verify_angular.py')],cwd=D,check=True)
assert hashlib.sha256((D/'angular_exact.json').read_bytes()).hexdigest()==old_hash
angular=json.loads((D/'angular_verified.json').read_text());assert angular['pass'] and angular['single_sign_transition'] and angular['all_u_ge4_tail_pass']
# The origin-degree-zero block [[MB,-M],[-M,M/B]] is rank-one PSD exactly.
assert M*B>0 and M/B>0 and (M*B)*(M/B)-M*M==0
# Complete sign rays. Root isolation accounts for every complex root, with multiplicity.
root_counts={}
for name,p,lo in [('fourier',H,Q(0)),('physical',pol,cut)]:
 roots=p.numer().complex_roots();assert sum(m for z,m in roots)==p.degree()
 forbidden=[(str(z),m)for z,m in roots if z.imag.contains(0) and not z.real<arb(lo)]
 assert not forbidden,(name,forbidden)
 root_counts[name]=sum(m for z,m in roots)
assert H(0)>0 and pol(cut)<0
# Pair inequality on the entire compact interval [pi*s,12]. Beyond12, P<0
# while0<h<=1, so the nonnegative pair allowance completes the infinite ray.
pi=arb.pi();DP=pol.derivative()-pol;DDP=DP.derivative()-DP
ars=[arb_poly([arb(x)for x in p.coeffs()])for p in [pol,DP,DDP]]
def tight(p,x):
 c=p(arb_poly([x.mid(),x.rad()]));err=sum((abs(y)for y in c.coeffs()[1:]),arb(0));return c[0]+arb(0,err.upper())
def residual(x,order=0,interval=False):
 pp=tight(ars[order],x)if interval else ars[order](x)
 z=x/(pi*arb(s))-1;h=(-arb(aa)*z-arb(bb)*z*z).exp()
 e1=-(arb(aa)+2*arb(bb)*z)/(pi*arb(s));e2=-2*arb(bb)/(pi*arb(s))**2
 if order==0:tail=-2*arb(M)*h+arb(M)*h*h
 elif order==1:tail=-2*arb(M)*e1*h+2*arb(M)*e1*h*h
 elif order==2:tail=-2*arb(M)*(e1*e1+e2)*h+arb(M)*(4*e1*e1+2*e2)*h*h
 else:raise AssertionError('Invalid derivative order')
 out=pp*(-x).exp()+tail;assert out.is_finite();return out
stack=[(Q(0),Q(1),0)];leaves=depth=0
while stack:
 a,b,dep=stack.pop();m=(a+b)/2;radius=(b-a)/2
 xm=pi*arb(s)+(arb(cut)-pi*arb(s))*arb(m)
 xr=(arb(cut)-pi*arb(s))*arb(radius)
 whole=pi*arb(s)+(arb(cut)-pi*arb(s))*arb(m,radius)
 f0=residual(xm);assert not f0>0,('Positive pair point',str(xm),str(f0))
 upper=f0+abs(residual(xm,1))*xr+abs(residual(whole,2,True))*xr*xr/2
 if upper<0:leaves+=1;depth=max(depth,dep);continue
 assert dep<50,('Unresolved pair box',a,b)
 stack.extend([(a,m,dep+1),(m,b,dep+1)])
obj=(pol(0)+M*B)/H(0);assert obj>0
rho=arb(obj)*(pi*arb(s)/4)**(arb(n)/2)/(arb(n)/2+1).gamma()
# Fully rational rounding from exact integer-dimensional ball volume.
def atan_series(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
pi_upper=16*atan_series(5,101)-4*atan_series(239,100)
sf=F(str(s));power=2**128;root_upper=F(math.isqrt(sf.numerator*power*power//sf.denominator)+1,power)
assert root_upper*root_upper>sf
k=n//2
if n%2:rational_upper=F(str(obj))*pi_upper**k*sf**k*root_upper/(2**k*math.prod(range(1,n+1,2)))
else:rational_upper=F(str(obj))*pi_upper**k*sf**k/(4**k*math.factorial(k))
target=F(r['strict_density_upper']);assert target==F('0.079503107') and rational_upper<target and rho<arb(Q(str(target)))
out={'pass':True,'dimension':n,'strict_density_upper':str(target),'density_interval':str(rho),'rational_rounding_upper':str(rational_upper),'exact_export':True,'exact_guard':'1000001/1000000','two_Fourier_reconstructions':True,'root_multiplicities_accounted':root_counts,'complete_Fourier_ray':True,'physical_tail_X':'12','complete_pair_intervals':leaves,'pair_depth':depth,'angular_full_cyclic_gate':angular,'origin_rank_one_PSD':True,'all_dilation_scope':'Curved monotone-difference profile and exact one-sign-transition Q; see THEOREM.md.','elapsed_seconds':time.time()-start}
(D/'full_replay_result.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items()if k not in ['rational_rounding_upper','angular_full_cyclic_gate']})
