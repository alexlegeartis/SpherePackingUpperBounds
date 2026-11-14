"""Rational export and rigorous global radial gates; no angular claim here."""
import argparse,json,sys,time
from fractions import Fraction
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as P, arb, arb_poly, ctx
from model import inputs,pi_bounds,density

if not __debug__:
    raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0)
ctx.prec=512
ap=argparse.ArgumentParser()
ap.add_argument('candidate')
ap.add_argument('--guard',default='10001/10000')
args=ap.parse_args()
path=Path(args.candidate).resolve();D=path.parent
can=json.loads(path.read_text());n,d=can['dimension'],can['degree']
cut=Q(str(Fraction(can['coupled_upper_X_exact'])))
p,h,ls,B=inputs(n,d)
coeff=[Q(str(Fraction(float(v)).limit_denominator(10**14))) for v in can['x']]
assert len(coeff)==d+3 and coeff[0]>=0 and coeff[-1]>0
p=coeff[0]*p+sum((coeff[k+1]*ls[k] for k in range(d+1)),P([]))
h=coeff[0]*h+sum(((-1)**k*coeff[k+1]*ls[k] for k in range(d+1)),P([]))
guard=Q(args.guard);assert guard>=1
M=coeff[-1]/B*guard
assert h(0)>0 and M>0
obj=(p(0)+B*M)/h(0)
out={'dimension':n,'degree':d,'certified':False,'physical_polynomial':list(map(str,p.coeffs())),
     'fourier_polynomial':list(map(str,h.coeffs())),'M_exact':str(M),'B_exact':str(B),
     'b0_exact':'1','M_guard_exact':str(guard),'objective_exact':str(obj),
     'coupled_upper_X_exact':str(cut),'physical_tail_X_exact':str(cut),
     'density_interval':str(density(obj,n))}
(D/'export.json').write_text(json.dumps(out,indent=2)+'\n')
bad={}
pl,pu=pi_bounds()
assert 6*pu<cut and Q(51,4)*pl>40
for name,poly,threshold in [('fourier',h,Q(0)),('physical',p,cut)]:
    roots=poly.numer().complex_roots()
    assert sum(m for z,m in roots)==poly.degree()
    bad[name]=[(float(z.real.mid()),m) for z,m in roots if z.imag.contains(0) and not z.real<arb(threshold)]
    print(name,'forbidden roots:',bad[name],flush=True)
(D/'root_result.json').write_text(json.dumps({'global_sign_roots':bad},indent=2)+'\n')
if any(bad.values()):
    raise RuntimeError('Forbidden sign roots')
assert h(0)>0 and p(cut)<0
AP=arb_poly([arb(c) for c in p.coeffs()])
stack=[(6*pl,cut,0)];passed=0;depth=0;start=time.time()
while stack:
    lo,hi,dep=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2
    terms=AP(arb_poly([arb(mid),arb(rad)]))
    up=terms[0]+sum((abs(c) for c in terms.coeffs()[1:]),arb(0))
    ok=up<0
    if not ok:
        cap=2-(-arb(Q(17,8))*(arb(lo)-6*arb(pu))).exp()
        expupper=(arb(Q(9,8))*arb(hi)-arb(Q(51,4))*arb(pl)).exp()
        ok=up*expupper<arb(M)*cap
    if ok:
        passed+=1;depth=max(depth,dep);continue
    if dep>=34:
        (D/'pair_failure.json').write_text(json.dumps({'lo':float(lo),'hi':float(hi)},indent=2)+'\n')
        raise RuntimeError('Pair coupling inconclusive or failed')
    stack.extend([(mid,hi,dep+1),(lo,mid,dep+1)])
out.update({'radial_checks_pass':True,'pair_intervals':passed,'pair_depth':depth,
            'pair_seconds':time.time()-start,'scope':'Radial gates only; angular/cyclic proof is required'})
(D/'radial_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('RADIAL PASS',n,d,out['density_interval'],passed,depth,flush=True)
