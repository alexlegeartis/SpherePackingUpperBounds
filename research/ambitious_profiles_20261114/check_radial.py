"""Exact Gaussian/Laguerre export, global rays, pair interval and density gates."""
import argparse,json,sys,math,time
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_CEILING
import numpy as np
from scipy.optimize import minimize_scalar
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
from common import source,rational,pi_bounds,ROOT,volume_arb
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=768
ap=argparse.ArgumentParser();ap.add_argument('candidate');ap.add_argument('--guard',default='1000001/1000000')
args=ap.parse_args();file=Path(args.candidate).resolve();D=file.parent;c=json.loads(file.read_text())
n,d=c['dimension'],c['degree'];s,a,b,q=map(rational,[c['s_exact'],c['a_exact'],c['b_exact'],c['q_exact']])
cut=Q(c['cut_X_exact']);assert s>0 and a+b>0 and b>0 and 0<q<2
rawp,rawh,ls=source(n,d);coeff=[Q(str(F(float(z)))) for z in c['x']]
assert len(coeff)==d+3 and coeff[0]>=0 and coeff[-1]>0
physical=coeff[0]*rawp+sum((coeff[k+1]*ls[k] for k in range(d+1)),P([]))
fourier=coeff[0]*rawh+sum(((-1)**k*coeff[k+1]*ls[k] for k in range(d+1)),P([]))
assert physical.degree()==fourier.degree() and fourier(0)>0
angular=json.loads((ROOT/c['angular_file']).read_text());B=Q(angular['B_exact'])
assert angular['dimension']==n and Q(angular['gegenbauer_coefficients'][0])==1 and B==Q(c['B_exact'])
guard=rational(args.guard);C=coeff[-1]*guard;assert guard>=1
obj=(physical(0)+C)/fourier(0);assert obj>0
out={'dimension':n,'degree':physical.degree(),'correction_degree':d,'s_exact':str(s),'a_exact':str(a),'b_exact':str(b),'q_exact':str(q),
     'C_exact':str(C),'B_exact':str(B),'C_guard_exact':str(guard),'cut_X_exact':str(cut),
     'physical_polynomial':list(map(str,physical.coeffs())),'fourier_polynomial':list(map(str,fourier.coeffs())),
     'objective_exact':str(obj),'angular_file':c['angular_file'],'certified':False}
(D/'radial_export.json').write_text(json.dumps(out,indent=2)+'\n')
pl,pu=pi_bounds();assert cut>s*pu
assert b+a*pu/pl>0 if a<0 else b>0
bad={}
for name,poly,endpoint in [('fourier',fourier,Q(0)),('physical',physical,cut)]:
    roots=poly.numer().complex_roots();assert sum(m for z,m in roots)==poly.degree()
    bad[name]=[(float(z.real.mid()),m) for z,m in roots if z.imag.contains(0) and not z.real<arb(endpoint)]
    print(name,'forbidden roots',bad[name],flush=True)
(D/'root_rejection.json').write_text(json.dumps({'global_sign_roots':bad},indent=2)+'\n')
if any(bad.values()):raise RuntimeError('Global radial sign gate failed')
assert fourier(0)>0 and physical(cut)<0
AP=arb_poly([arb(z) for z in physical.coeffs()]);pi=arb.pi();M=arb(C/B)
def profile(v):return arb(q)*(-arb(a)*v.log()-arb(b)*(v-1)).exp()
assert profile(arb(pl)/arb(pu))<2
def ratio(x):
    xx=arb(float(x));hh=profile(xx/(pi*s));den=M*(2*hh-hh*hh)
    assert den>0
    return float((AP(xx)*(-xx).exp()/den).mid())
grid=np.linspace(float(s)*np.pi,float(cut),4001);vals=np.array([ratio(x) for x in grid]);peaks=[]
for i in range(1,len(grid)-1):
    if vals[i]>vals[i-1] and vals[i]>vals[i+1]:
        opt=minimize_scalar(lambda x:-ratio(x),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-12})
        if -opt.fun>1:
            xx=arb(float(opt.x));hh=profile(xx/(pi*s));gap=AP(xx)*(-xx).exp()-M*(2*hh-hh*hh)
            if gap>0:peaks.append(float(opt.x))
(D/'pair_rejection.json').write_text(json.dumps({'pair_points':peaks},indent=2)+'\n')
if peaks:
    print('PAIR COUNTEREXAMPLES',peaks,flush=True);raise RuntimeError('Pair gate failed')
stack=[(s*pl,cut,0)];passed=depthmax=0;start=time.time()
while stack:
    lo,hi,depth=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2
    terms=AP(arb_poly([arb(mid),arb(rad)]));up=terms[0]+sum((abs(z) for z in terms.coeffs()[1:]),arb(0))
    okay=up<0
    if not okay:
        hmin=profile(arb(hi)/(arb(pl)*s));hmax=profile(arb(lo)/(arb(pu)*s))
        capmin=arb(min((2*hmin-hmin*hmin).lower(),(2*hmax-hmax*hmax).lower()))
        okay=up*(-arb(lo)).exp()<M*capmin
    if okay:passed+=1;depthmax=max(depthmax,depth);continue
    if depth>=36:
        xx=arb(mid);hh=profile(xx/(pi*s));gap=AP(xx)*(-xx).exp()-M*(2*hh-hh*hh)
        if gap>0:
            (D/'pair_rejection.json').write_text(json.dumps({'pair_points':[float(mid)],'gap':str(gap)},indent=2)+'\n')
        raise RuntimeError('Pair interval proof inconclusive or failed')
    stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
density=arb(obj)*volume_arb(n,s)
with localcontext() as dc:
    dc.prec=80
    value=Decimal(str(math.nextafter(float(density.upper()),math.inf)))
    quantum=Decimal(1).scaleb(value.adjusted()-10)
    strict_decimal=value.quantize(quantum,rounding=ROUND_CEILING)
strict=Q(str(F(str(strict_decimal))))
if n%2==0:assert 0<obj*(pu*s/4)**(n//2)/math.factorial(n//2)<strict
else:
    k=(n-1)//2
    assert 0<obj**2*4*pu**(n-1)*s**n*Q(math.factorial(k+1),math.factorial(n+1))**2<strict**2
rows=json.loads((ROOT/'RESULTS.json').read_text());old=next(row for row in rows if row['dimension']==n)
out.update({'radial_checks_pass':True,'strict_density_upper_bound':str(strict_decimal),'density_interval':str(density),
            'pair_intervals':passed,'pair_depth':depthmax,'pair_seconds':time.time()-start,
            'improves_registered':strict<Q(str(F(old['strict_density_upper_bound']))),
            'scope':'Radial/PSD-origin/rounding gates; reference angular/cyclic replay required for density conclusion'})
(D/'radial_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('RADIAL ALL GATES PASS',n,str(density),str(strict_decimal),passed,depthmax,flush=True)
