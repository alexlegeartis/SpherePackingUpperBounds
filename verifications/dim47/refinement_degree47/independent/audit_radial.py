"""Independent new radial-refit audit and exact density reconstruction.

Fourier reconstruction uses the Laguerre three-term recurrence, independently
of the primary verifier's coefficient-sum formula. Pair coupling uses the
signed full interval expression. All output is confined to this audit folder.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import sys
import time
from flint import fmpq as Q, fmpq_poly as P, arb, arb_poly, ctx

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0)
ctx.prec=512
D=Path(__file__).resolve().parent
n=47
new=json.loads((D/'inputs/certificate.json').read_text())
ang=json.loads((D/'inputs/angular47.json').read_text())
source=json.loads((D/'inputs/source48_radial.json').read_text())
candidate=json.loads((D/'inputs/discovery_candidate.json').read_text())
assert new['dimension']==ang['dimension']==candidate['dimension']==n
assert candidate['degree']==47
sha=lambda b:hashlib.sha256(b).hexdigest()
physical=P(list(map(Q,new['physical_polynomial'])))
fourier=P(list(map(Q,new['fourier_polynomial'])))
M=Q(new['M_exact'])
B=Q(ang['B_exact'])
b0=Q(new['b0_exact'])
assert b0==Q(ang['gegenbauer_coefficients'][0])==1 and M>0 and fourier(0)>0
assert physical.degree()==fourier.degree()==111

# Independently reconstruct the frozen refit from normalized L_l^(alpha-1)(2x).
alpha=Q(n,2)
x=P([0,1])
coeff=[Q(str(Fraction(float(v)).limit_denominator(10**14))) for v in candidate['x']]
assert len(coeff)==50 and coeff[0]>=0 and coeff[-1]>0
source_P=P(list(map(Q,source['physical_polynomial'])))
source_N=Q(source['fourier_polynomial'][0]);assert source_N>0
base=coeff[0]*source_P/source_N
L0,L1=P([1]),alpha-2*x
correction=coeff[1]*L0+coeff[2]*L1/L1(0)
correction_F=coeff[1]*L0-coeff[2]*L1/L1(0)
for k in range(1,47):
    L2=((2*k+alpha-2*x)*L1-(k+alpha-1)*L0)/(k+1)
    correction+=coeff[k+2]*L2/L2(0)
    correction_F+=(-1)**(k+1)*coeff[k+2]*L2/L2(0)
    L0,L1=L1,L2
assert correction.degree()==47 and base+correction==physical
assert M==coeff[-1]/B*Q(1000001,1000000)
added=[k for k in range(32,48) if coeff[k+1]]
assert added and source_P.degree()==111 and correction.degree()>31
print('REFIT_RECONSTRUCTION_AND_M_GUARD_PASS',flush=True)

# Fourier[e^{-pi r^2}(pi r^2)^k]=k! L_k^(alpha-1)(pi rho^2)e^{-pi rho^2}.
x=P([0,1]);alpha=Q(47,2)
L0,L1=P([1]),alpha-x
transformed=physical[0]*L0+physical[1]*L1
source_transformed=source_P[0]*L0+source_P[1]*L1
factorial=1
for k in range(1,111):
    L2=((2*k+alpha-x)*L1-(k+alpha-1)*L0)/(k+1)
    factorial*=k+1
    transformed+=physical[k+1]*factorial*L2
    source_transformed+=source_P[k+1]*factorial*L2
    L0,L1=L1,L2
assert transformed==fourier
assert coeff[0]*source_transformed/source_N+correction_F==fourier
print('RADIAL_FOURIER_IDENTITY_PASS',flush=True)


def arctan_bounds(x):
    terms=128
    lower=sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(terms)),Q(0))
    return lower,lower+x**(2*terms+1)/(2*terms+1)


a,b=arctan_bounds(Q(1,5))
c,d=arctan_bounds(Q(1,239))
pl,pu=16*a-4*d,16*b-4*c
assert pl<pu and Q(51,4)*pl>40
root_results={}
for name,poly,cut in [('fourier',fourier,Q(0)),('physical',physical,7*pl)]:
    roots=poly.numer().complex_roots()
    assert sum(m for z,m in roots)==111
    forbidden=[str(z) for z,m in roots if z.imag.contains(0) and not z.real<arb(cut)]
    assert not forbidden,(name,forbidden)
    root_results[name]={'degree_count':sum(m for z,m in roots),'forbidden_roots':[],
                        'cutoff_exact':str(cut),'all_roots_accounted_for':True}
    print(name.upper()+'_RAY_PASS',flush=True)
assert fourier(0)>0 and physical(7*pl)<0

AP=arb_poly([arb(c) for c in physical.coeffs()])
pi=arb.pi()
stack=[(6*pl,7*pu,0)]
leaves=[]
max_depth=0
start=time.time()
while stack:
    lo,hi,depth=stack.pop()
    mid,rad=(lo+hi)/2,(hi-lo)/2
    terms=AP(arb_poly([arb(mid),arb(rad)]))
    radius=sum((abs(c) for c in terms.coeffs()[1:]),arb(0))
    pbound=terms[0]+arb(0,radius.upper())
    xx=arb(mid,rad)
    gap=pbound*(Q(9,8)*xx-Q(51,4)*pi).exp()-arb(M)*(2-(-Q(17,8)*(xx-6*pi)).exp())
    assert gap.is_finite()
    if gap<0:
        leaves.append({'lo':str(lo),'hi':str(hi),'depth':depth,'coupling_gap_upper':str(gap.upper())})
        max_depth=max(depth,max_depth)
        continue
    assert depth<35,(str(lo),str(hi),str(gap))
    stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])

# Independently check that the accepted leaves tile the enclosing interval.
ordered=sorted(leaves,key=lambda z:Q(z['lo']))
assert Q(ordered[0]['lo'])==6*pl and Q(ordered[-1]['hi'])==7*pu
assert all(Q(a['hi'])==Q(b['lo']) for a,b in zip(ordered,ordered[1:]))
(D/'pair_interval_cover.json').write_text(json.dumps(ordered,indent=2)+'\n')
print('COUPLING_PASS',len(leaves),'depth',max_depth,'seconds',time.time()-start,flush=True)

origin=M*B
degree_zero=M/B
assert origin*degree_zero==M*M and origin>0 and degree_zero>0
objective=(physical(0)+origin)/fourier(0)
assert objective==Q(new['objective_exact'])
strict=Q(5998306274,10**16)
old_strict=Q(6010083731,10**16)
comparison=Q(6490757338,10**16)
# Gamma(49/2)=48!*sqrt(pi)/(4^24*24!).  Squaring removes sqrt(6).
volume_factor=Q(4**24*math.factorial(24),math.factorial(48))
upper_squared=objective*objective*pu**46*Q(3,2)**47*volume_factor**2
lower_squared=objective*objective*pl**46*Q(3,2)**47*volume_factor**2
assert Q(3,2)**47*volume_factor**2 == Q(4*6**47*math.factorial(24)**2,math.factorial(48)**2)
assert objective>0 and 0<lower_squared<upper_squared<strict**2
assert strict<old_strict<comparison
density=arb(objective)*(Q(3,2)*pi)**arb(Q(47,2))/arb(Q(49,2)).gamma()
assert density<arb(strict)
result={'pass':True,'dimension':47,'precision_bits':ctx.prec,
        'new_radial_refit_reconstructed_exact':True,'correction_degree':correction.degree(),
        'exact_Fourier_span_membership':True,'source_t_nonnegative':True,
        'nonzero_added_Laguerre_degrees':added,'outside_degree31_span_exact':True,
        'M_guard_exact':'1000001/1000000',
        'source_radial_file_sha256':sha((D/'inputs/source48_radial.json').read_bytes()),
        'discovery_candidate_sha256':sha((D/'inputs/discovery_candidate.json').read_bytes()),
        'new_certificate_file_sha256':sha((D/'inputs/certificate.json').read_bytes()),
        'angular_file_sha256':sha((D/'inputs/angular47.json').read_bytes()),
        'Fourier_identity_by_Laguerre_recurrence':True,'root_results':root_results,
        'coupling_intervals':len(leaves),'coupling_depth':max_depth,
        'coupling_interval_union_checked':True,'coupling_full_signed_expression':True,
        'b0_exact':str(b0),'M_exact':str(M),'B_exact':str(B),
        'origin_cost_exact':str(origin),'origin_PSD_rank_one_identity':True,
        'objective_exact':str(objective),'objective_matches_primary_exact':True,
        'damping_lower_exceeds_40_exact':True,
        'strict_density_upper':'0.0000005998306274',
        'density_interval':str(density),'rational_rounding_strict':True,'odd_dimensional_squared_volume_identity_exact':True,
        'pinned_comparator':'0.0000006490757338',
        'fractional_improvement_lower':str(1-strict/comparison),
        'scope':'New radial refit, exact origin cost and rounding. Angular/cyclic result files also required.'}
(D/'radial_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('RADIAL_REFIT_PASS',json.dumps({k:v for k,v in result.items() if k not in ['objective_exact','root_results']}),flush=True)
