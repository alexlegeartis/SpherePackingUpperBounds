"""Exact coefficient-sum Fourier identity, complete root and pair proof."""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as P, arb, arb_poly, ctx

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).resolve().parent
r=json.loads((D/'certificate.json').read_text())
cfg=json.loads((D/'config.json').read_text())
ang=json.loads((D/'continuum_three_point/cyclic_angular40_exact.json').read_text())
n=cfg['dimension'];assert r['dimension']==ang['dimension']==n and n in [45,46,47]
p=list(map(F,r['physical_polynomial']));h=list(map(F,r['fourier_polynomial']))
M=F(r['M_exact']);B=F(ang['B_exact']);b0=F(r['b0_exact']);a=F(n,2)
assert len(p)==len(h)==112 and h[0]>0 and M>0 and b0==1
assert F(r['B_exact'])==B and F(ang['gegenbauer_coefficients'][0])==b0
rise=lambda u,k:math.prod((u+j for j in range(k)),start=F(1))
transform=[sum((p[k]*math.factorial(k)*(-1)**j*rise(a+j,k-j)/
                math.factorial(k-j)/math.factorial(j) for k in range(j,len(p))),F(0))
           for j in range(len(p))]
assert transform==h
fq=lambda x:Q(x.numerator,x.denominator)
af=lambda x:arb(fq(x))
physical=P(list(map(fq,p)));fourier=P(list(map(fq,h)))
def atan(x):
    lo=sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(110)),F(0))
    return lo,lo+x**221/221
al,au=atan(F(1,5));bl,bu=atan(F(1,239));pl,pu=16*al-4*bu,16*au-4*bl
cut=F(r['physical_tail_X_exact'])
assert cut==F(r['coupled_upper_X_exact'])==40 and F(51,4)*pl>40
for name,poly,endpoint in [('fourier',fourier,F(0)),('physical',physical,cut)]:
    roots=poly.numer().complex_roots()
    assert sum(m for z,m in roots)==poly.degree()==111
    assert all(not z.imag.contains(0) or z.real<af(endpoint) for z,m in roots)
    print(name,'global root gate PASS',flush=True)
assert fourier(0)>0 and physical(fq(cut))<0
ap=arb_poly([af(c) for c in p]);stack=[(6*pl,cut,0)];passed=maxdepth=0
while stack:
    lo,hi,depth=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2
    terms=ap(arb_poly([af(mid),af(rad)]))
    up=terms[0]+sum((abs(c) for c in terms.coeffs()[1:]),arb(0))
    okay=up<0
    if not okay:
        cap=2-(-arb(Q(17,8))*(af(lo)-6*af(pu))).exp()
        expupper=(arb(Q(9,8))*af(hi)-arb(Q(51,4))*af(pl)).exp()
        okay=up*expupper<af(M)*cap
    if okay:
        passed+=1;maxdepth=max(maxdepth,depth);continue
    assert depth<35, (str(lo),str(hi))
    stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
origin=M*B/b0;degree0=M*b0/B
assert origin>0 and degree0>0 and origin*degree0==M*M
objective=(p[0]+origin)/h[0];assert objective==F(r['objective_exact'])>0
strict=F(cfg['strict_density_upper_bound'])
if n%2==0:
    upper=objective*(F(3,2)*pu)**(n//2)/math.factorial(n//2)
    assert 0<upper<strict
else:
    m=(n-1)//2
    # V_n(sqrt(6)/2)^2 = 4*pi^(n-1)*6^n*((m+1)!)^2/((n+1)!)^2.
    upper2=objective**2*4*pu**(n-1)*6**n*math.factorial(m+1)**2/math.factorial(n+1)**2
    assert 0<upper2<strict**2
assert strict<F(cfg['previous_repository_strict_upper'])
density=af(objective)*(arb.pi()*arb(3)/2)**(arb(n)/2)/(arb(n)/2+1).gamma()
assert density<af(strict)
out={'pass':True,'dimension':n,'strict_density_upper_bound':str(strict),
     'strict_density_upper_decimal':cfg['strict_density_upper_bound'],
     'density_interval':str(density),'exact_Fourier_coefficient_sum':True,
     'complete_global_root_accounting':True,'coupled_interval_X':'[6*pi,40]',
     'physical_negative_ray_X':'[40,infinity)','pair_intervals':passed,'pair_depth':maxdepth,
     'origin_rank_one_identity_exact':True,'positive_Fourier_normalization':True,
     'rational_Machin_rounding':True,'scope':'Primary radial gates; angular and cyclic checks are separate'}
(D/'primary_radial_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('PRIMARY RADIAL PASS',n,str(density),passed,maxdepth,flush=True)
