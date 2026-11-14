"""Primary direct Fourier identity, complete rays, pair domain and strict rounding."""
import argparse,json,math,sys,time
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
from core import inputs,rational,pi_bounds,rational_rounding
if not __debug__:raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=512
ap=argparse.ArgumentParser();ap.add_argument('package');args=ap.parse_args()
D,r,cfg,ang,n,s,a,b,q=inputs(args.package)
pc=list(map(F,r['physical_polynomial']));hc=list(map(F,r['fourier_polynomial']))
assert len(pc)==len(hc) and len(pc)>1
alpha=F(n,2)
rise=lambda x,k:math.prod((x+j for j in range(k)),start=F(1))
got=[sum((pc[k]*math.factorial(k)*(-1)**j*rise(alpha+j,k-j)/math.factorial(k-j)/math.factorial(j)
          for k in range(j,len(pc))),F(0)) for j in range(len(pc))]
assert got==hc
physical=P(list(map(Q,r['physical_polynomial'])));fourier=P(list(map(Q,r['fourier_polynomial'])))
C,B,cut=map(rational,[r['C_exact'],r['B_exact'],r['cut_X_exact']])
assert C>0 and B==Q(ang['B_exact']) and Q(ang['gegenbauer_coefficients'][0])==1
assert fourier(0)>0 and physical.degree()==fourier.degree()==r['degree']
pl,pu=pi_bounds();assert cut>s*pu
if a<0:assert b+a*pu/pl>0
root_results={}
for name,poly,endpoint in [('fourier',fourier,Q(0)),('physical',physical,cut)]:
    roots=poly.numer().complex_roots();assert sum(m for z,m in roots)==poly.degree()
    assert all(not z.imag.contains(0) or z.real<arb(endpoint) for z,m in roots)
    root_results[name]={'degree':poly.degree(),'all_roots_accounted_for':True,'forbidden_real_roots':0,'cutoff_exact':str(endpoint)}
    print('PRIMARY',name,'global ray PASS',flush=True)
assert physical(cut)<0
AP=arb_poly([arb(z) for z in physical.coeffs()]);M=arb(C/B)
def h(v):return arb(q)*(-arb(a)*v.log()-arb(b)*(v-1)).exp()
assert h(arb(pl)/arb(pu))<2
stack=[(s*pl,cut,0)];accepted=depthmax=0;start=time.time()
while stack:
    lo,hi,depth=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2
    terms=AP(arb_poly([arb(mid),arb(rad)]))
    upper=terms[0]+sum((abs(z) for z in terms.coeffs()[1:]),arb(0))
    okay=upper<0
    if not okay:
        hmin=h(arb(hi)/(arb(pl)*s));hmax=h(arb(lo)/(arb(pu)*s))
        allowance_lower=arb(min((2*hmin-hmin*hmin).lower(),(2*hmax-hmax*hmax).lower()))
        okay=upper*(-arb(lo)).exp()<M*allowance_lower
    if okay:accepted+=1;depthmax=max(depthmax,depth);continue
    assert depth<38,(str(lo),str(hi))
    stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
origin=C;cross=-q*C/B;nonorigin=q*q*C/(B*B)
assert origin>0 and nonorigin>0 and origin*nonorigin==cross*cross
objective=(physical(0)+C)/fourier(0);assert objective==rational(r['objective_exact'])
strict=rational(cfg['strict_density_upper_bound']);rational_rounding(objective,n,s,pu,strict)
density=arb(objective)*(arb.pi()*arb(s)/4)**(arb(n)/2)/(arb(n)/2+1).gamma()
assert density<arb(strict)
out={'pass':True,'dimension':n,'degree':physical.degree(),'precision_bits':ctx.prec,
     'direct_Fourier_coefficient_sum_exact':True,'root_results':root_results,
     'complete_pair_domain':True,'pair_intervals':accepted,'pair_depth':depthmax,'seconds':time.time()-start,
     'origin_Gram_identity_exact':True,'positive_Fourier_normalization':True,
     'rational_Machin_rounding':True,'strict_density_upper_bound':cfg['strict_density_upper_bound'],
     'density_interval':str(density),'scope':'Primary radial gates; angular and cyclic replay also required'}
(D/'primary_radial_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('PRIMARY RADIAL PASS',n,cfg['strict_density_upper_bound'],accepted,depthmax,flush=True)
