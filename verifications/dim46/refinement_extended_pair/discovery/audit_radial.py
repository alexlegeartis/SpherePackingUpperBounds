"""Exact export and radial proof attempt for sampled perturbation candidates."""
import sys,json,math,argparse
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).parent;S=D/'inputs'
ap=argparse.ArgumentParser();ap.add_argument('candidate');args=ap.parse_args();cand=json.loads(Path(args.candidate).read_text());n=cand['dimension'];d=cand['degree'];r=json.loads((S/'adaptive_cyclic40_radial111.json').read_text());ang=json.loads((D/'cyclic_angular40_exact.json').read_text())
q=lambda x:fmpq(x.numerator,x.denominator);af=lambda x:arb(q(x))
rf=lambda a,k:math.prod((a+j for j in range(k)),start=F(1))
N48=F(r['fourier_polynomial'][0]);B=F(ang['B_exact']);a=F(n,2);x=[F(float(v)).limit_denominator(10**14) for v in cand['x']];p=[x[0]*F(v)/N48 for v in r['physical_polynomial']]
for l in range(d+1):
 norm=rf(a,l)/math.factorial(l)
 for j in range(l+1):p[j]+=x[l+1]*(-2)**j*rf(a+j,l-j)/math.factorial(l-j)/math.factorial(j)/norm
h=[sum((p[k]*math.factorial(k)*(-1)**j*rf(a+j,k-j)/math.factorial(k-j)/math.factorial(j) for k in range(j,len(p))),F(0)) for j in range(len(p))]
cs=list(map(F,ang['polynomial_coefficients']));b0=sum((cs[2*j]*rf(F(1,2),j)/rf(a,j) for j in range(len(cs)//2)),F(0));Bn=B/b0;M=x[-1]/Bn*F(10001,10000)
assert h[0]>0 and M>0
P=fmpq_poly(list(map(q,p)));H=fmpq_poly(list(map(q,h)))
(D/'unverified_export.json').write_text(json.dumps({'physical_polynomial':list(map(str,p)),'fourier_polynomial':list(map(str,h)),'M_exact':str(M),'B_exact':str(B),'certified':False},indent=2))
def ab(x):
 lo=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(110));return lo,lo+x**221/221
al,au=ab(F(1,5));bl,bu=ab(F(1,239));pl=16*al-4*bu;pu=16*au-4*bl
bad={}
for name,poly,cut in [('fourier',H,F(0)),('physical',P,F(40))]:
 print('Root audit',name,flush=True);roots=poly.numer().complex_roots();assert sum(m for z,m in roots)==111
 forbidden=[(float(z.real.mid()),m) for z,m in roots if z.imag.contains(0) and not z.real<af(cut)]
 bad[name]=forbidden;print(name,forbidden,flush=True)
if any(bad.values()):
 (D/(Path(args.candidate).stem+'_rejection.json')).write_text(json.dumps({'global_sign_roots':bad},indent=2));raise RuntimeError('Forbidden sign roots')
assert H(0)>0 and P(q(F(40)))<0
AP=arb_poly([af(v) for v in p]);stack=[(6*pl,F(40),0)];passed=0;maxdepth=0
while stack:
 l,u,dep=stack.pop();mid=(l+u)/2;rad=(u-l)/2
 co=AP(arb_poly([af(mid),af(rad)]));up=co[0]+sum((abs(co[j]) for j in range(1,len(co))),arb(0));ok=up<0
 if not ok:
  cap=2-(-arb(fmpq(17,8))*(af(l)-6*af(pu))).exp()
  expupper=(arb(fmpq(9,8))*af(u)-arb(fmpq(51,4))*af(pl)).exp()
  ok=up*expupper<af(M)*cap
 if ok:passed+=1;maxdepth=max(maxdepth,dep);continue
 if dep>=32:
  print('COUPLING_FAILURE',float(l),float(u),flush=True);raise RuntimeError('coupling')
 stack.extend([(mid,u,dep+1),(l,mid,dep+1)])
obj=(p[0]+Bn*M)/h[0];density=af(obj)*(arb.pi()*arb(3)/2)**(arb(n)/2)/(arb(n)/2+1).gamma()
out={'pass':True,'dimension':n,'density_interval':str(density),'objective_exact':str(obj),'M_exact':str(M),'b0_exact':str(b0),'physical_polynomial':[str(v) for v in p],'fourier_polynomial':[str(v) for v in h],'M_guard_exact':'10001/10000','coupled_upper_X_exact':'40','pair_intervals':passed,'pair_depth':maxdepth,'source_candidate':Path(args.candidate).name,'scope':'Radial exact checks only; inherited angular and cyclic proofs required'}
(D/(Path(args.candidate).stem+'_certificate.json')).write_text(json.dumps(out,indent=2));print('RADIAL_PASS',density,passed,maxdepth,flush=True)
