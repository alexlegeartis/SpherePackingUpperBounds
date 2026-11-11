"""Exact radial export, rigorous roots, full scalar pair interval."""
from pathlib import Path
from fractions import Fraction as F
import json,math,argparse,time
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
ctx.prec=384;D=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('candidate');a=ap.parse_args();path=Path(a.candidate);j=json.loads(path.read_text());n=j['dimension'];d=j['degree'];s=Q(j['s_exact']);aa=Q(j['a_exact']);bb=Q(j['b_exact']);cut=Q(str(F(j['cut_X'])));ang=json.loads(Path(j['angular_source']).read_text());B=Q(ang['B_exact']);X=P([0,1]);v=[Q(str(F(x).limit_denominator(10**16))) for x in j['radial']];C=Q(str(F(j['C']).limit_denominator(10**16)));alpha=Q(n,2)
E=[P([1]),alpha-2*X]
for l in range(1,d):E.append(((2*l+alpha-2*X)*E[-1]-(l+alpha-1)*E[-2])/(l+1))
E=[p/p(0)for p in E];pol=sum((x*p for x,p in zip(v,E)),P([]));H=sum(((-1)**l*x*p for l,(x,p)in enumerate(zip(v,E))),P([]));assert H(0)>0
M=C/B*Q(1000001,1000000);assert M>0
out={'dimension':n,'degree':d,'s':str(s),'a_exact':str(aa),'b_exact':str(bb),'cut_X':str(cut),'B_exact':str(B),'M_exact':str(M),'physical':list(map(str,pol.coeffs())),'fourier':list(map(str,H.coeffs())),'origin_cost_exact':str(M*B),'normalization_exact':str(H(0)),'scope':'Radial export. A PASS still requires separate angular and cyclic gates.'}
(D/(path.stem+'_exact.json')).write_text(json.dumps(out,indent=2)+'\n')
reject={};st=time.time()
for name,p,lo in [('fourier',H,Q(0)),('physical',pol,cut)]:
 print('Root check',name,p.degree(),flush=True);roots=p.numer().complex_roots();assert sum(m for z,m in roots)==p.degree()
 bad=[(float(z.real.mid()),m)for z,m in roots if z.imag.contains(0) and not z.real<arb(lo)];reject[name]=bad
 print(name,bad,'seconds',time.time()-st,flush=True)
(D/(path.stem+'_rejection.json')).write_text(json.dumps({'global_sign_roots':reject},indent=2)+'\n')
if any(reject.values()):raise RuntimeError('Forbidden radial roots; add cuts before claiming a bound')
assert H(0)>0 and pol(cut)<0
AP=arb_poly([arb(x)for x in pol.coeffs()]);pi=arb.pi()
def tight(p,x):
 c=p(arb_poly([x.mid(),x.rad()]));err=sum((abs(y)for y in c.coeffs()[1:]),arb(0));return c[0]+arb(0,err.upper())
DP=pol.derivative()-pol;DDP=DP.derivative()-DP
ars=[AP,arb_poly([arb(x)for x in DP.coeffs()]),arb_poly([arb(x)for x in DDP.coeffs()])]
def residual(x,order=0,interval=False):
 pp=tight(ars[order],x)if interval else ars[order](x)
 z=x/(pi*arb(s))-1;h=(-arb(aa)*z-arb(bb)*z*z).exp()
 e1=-(arb(aa)+2*arb(bb)*z)/(pi*arb(s));e2=-2*arb(bb)/(pi*arb(s))**2
 if order==0: tail=-2*arb(M)*h+arb(M)*h*h
 elif order==1:tail=-2*arb(M)*e1*h+2*arb(M)*e1*h*h
 elif order==2:tail=-2*arb(M)*(e1*e1+e2)*h+arb(M)*(4*e1*e1+2*e2)*h*h
 else:raise AssertionError('order')
 return pp*(-x).exp()+tail
stack=[(Q(0),Q(1),0)];leaves=depth=0
while stack:
 a,b,dep=stack.pop();m=(a+b)/2;r=(b-a)/2;xm=pi*arb(s)+(arb(cut)-pi*arb(s))*arb(m);xr=(arb(cut)-pi*arb(s))*arb(r);whole=pi*arb(s)+(arb(cut)-pi*arb(s))*arb(m,r)
 f0=residual(xm)
 if f0>0:
  (D/(path.stem+'_pair_rejection.json')).write_text(json.dumps({'X':float(xm.mid()),'positive_residual':str(f0)},indent=2)+'\n')
  raise RuntimeError(f'Positive pair point X={xm}: {f0}')
 upper=f0+abs(residual(xm,1))*xr+abs(residual(whole,2,True))*xr*xr/2
 if upper<0:leaves+=1;depth=max(depth,dep);continue
 if dep>=50:raise RuntimeError(f'Pair unresolved {a} {b}, {upper}')
 stack.extend([(a,m,dep+1),(m,b,dep+1)])
obj=(pol(0)+M*B)/H(0);rho=arb(obj)*(pi*arb(s)/4)**(arb(n)/2)/(arb(n)/2+1).gamma()
report={'pass':True,'density_interval':str(rho),'objective_exact':str(obj),'pair_leaves':leaves,'pair_depth':depth,'elapsed':time.time()-st,'scope':'Radial gates only; not density bound until angular and cyclic proofs pass.'}
(D/(path.stem+'_radial_pass.json')).write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
