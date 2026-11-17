"""Independent Laguerre Fourier reconstruction and signed full-interval pair proof.

Root isolation shares FLINT/Arb with the primary path; no separate backend is
claimed. Fourier reconstruction, pair enclosures and partition checks differ.
"""
import argparse,json,math,sys,time,hashlib
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=1024
ap=argparse.ArgumentParser();ap.add_argument('package');args=ap.parse_args();D=Path(args.package).resolve()
r=json.loads((D/'certificate.json').read_text());cfg=json.loads((D/'config.json').read_text());ang=json.loads((D/'angular.json').read_text())
qf=lambda z:Q(str(Fraction(str(z))))
n=r['dimension'];assert n==cfg['dimension']==ang['dimension'] and n>=3
s,a,b,q=[qf(r[k]) for k in ['s_exact','a_exact','b_exact','q_exact']]
assert s>0 and a+b>0 and b>0 and q==1
assert a==qf(ang['profile_a_exact']) and b==qf(ang['profile_b_exact'])
p=P(list(map(Q,r['physical_polynomial'])));h=P(list(map(Q,r['fourier_polynomial'])))
C,B,cut=[qf(r[k]) for k in ['C_exact','B_exact','cut_X_exact']]
assert C>0 and B==Q(ang['B_exact']) and Q(ang['gegenbauer_coefficients'][0])==1
assert p.degree()==h.degree()==r['degree']
alpha=Q(n,2);x=P([0,1]);L0,L1=P([1]),alpha-x
reconstructed=p[0]*L0+p[1]*L1;factorial=1
for k in range(1,p.degree()):
    L2=((2*k+alpha-x)*L1-(k+alpha-1)*L0)/(k+1)
    factorial*=k+1;reconstructed+=p[k+1]*factorial*L2;L0,L1=L1,L2
assert reconstructed==h and h(0)>0 and p(cut)<0
def atan(z):
    terms=160
    lo=sum(((-1)**k*z**(2*k+1)/(2*k+1) for k in range(terms)),Q(0))
    return lo,lo+z**(2*terms+1)/(2*terms+1)
al,au=atan(Q(1,5));bl,bu=atan(Q(1,239));pl,pu=16*al-4*bu,16*au-4*bl
assert s*pu<cut
if a<0:assert b+a*pu/pl>0
for poly,endpoint in [(h,Q(0)),(p,cut)]:
    roots=poly.numer().complex_roots();assert sum(m for z,m in roots)==poly.degree()
    assert not [z for z,m in roots if z.imag.contains(0) and not z.real<arb(endpoint)]
def profile(v):
    power=v**int(-a) if a.denominator==1 else v**arb(-a)
    return arb(q)*power*(-arb(b)*(v-1)).exp()
assert profile(arb(pl)/arb(pu))<2
AP=arb_poly([arb(z) for z in p.coeffs()]);M=arb(C/B);stack=[(s*pl,cut,0)]
leaves=[];depthmax=0;start=time.time()
while stack:
    lo,hi,depth=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2
    coeff=AP(arb_poly([arb(mid),arb(rad)]));error=sum((abs(z) for z in coeff.coeffs()[1:]),arb(0))
    pbox=coeff[0]+arb(0,error.upper());xx=arb(mid,rad)
    lowh=profile(arb(hi)/(arb(pl)*s));highh=profile(arb(lo)/(arb(pu)*s))
    hl,hu=lowh.lower(),highh.upper()
    rlo=2*arb(hl)-arb(hl)**2;rhi=2*arb(hu)-arb(hu)**2
    allowance_low=min(rlo.lower(),rhi.lower())
    allowance_high=arb(1) if hl<=1<=hu else max(rlo.upper(),rhi.upper())
    allowance=arb(allowance_low).union(arb(allowance_high))
    gap=pbox*(-xx).exp()-M*allowance
    assert gap.is_finite()
    if gap<0:leaves.append((lo,hi));depthmax=max(depthmax,depth);continue
    assert depth<39,(str(lo),str(hi),str(gap))
    stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
leaves.sort();assert leaves[0][0]==s*pl and leaves[-1][1]==cut
assert all(x[1]==y[0] for x,y in zip(leaves,leaves[1:]))
objective=(p(0)+C)/h(0);assert objective==qf(r['objective_exact'])>0
assert C*q*q*C/(B*B)==(-q*C/B)**2
strict=qf(cfg['strict_density_upper_bound']);assert strict==qf(r['strict_density_upper_bound'])>0
if n%2==0:assert objective*(s*pu/4)**(n//2)/math.factorial(n//2)<strict
else:
    k=(n-1)//2
    assert objective**2*4*pu**(n-1)*s**n*Q(math.factorial(k+1),math.factorial(n+1))**2<strict**2
density=arb(objective)*(arb.pi()*arb(s)/4)**(arb(n)/2)/(arb(n)/2+1).gamma();assert density<arb(strict)
partition='\n'.join(f'{lo} {hi}' for lo,hi in leaves).encode()
out={'pass':True,'dimension':n,'precision_bits':ctx.prec,'Laguerre_recurrence_Fourier_exact':True,
     'global_rays_all_roots_accounted':True,'signed_pair_interval_proof':True,
     'complete_pair_partition_checked':True,'pair_partition_sha256':hashlib.sha256(partition).hexdigest(),
     'pair_intervals':len(leaves),'pair_depth':depthmax,'seconds':time.time()-start,
     'origin_PSD_identity_exact':True,'positive_Fourier_normalization':True,'rational_Machin_rounding':True,
     'strict_density_upper_bound':cfg['strict_density_upper_bound'],'density_interval':str(density),
     'scope':'Independent radial path; independent angular/cyclic path also required'}
(D/'independent_radial_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('INDEPENDENT RADIAL PASS',n,cfg['strict_density_upper_bound'],len(leaves),depthmax,flush=True)
