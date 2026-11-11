"""Independent unchanged-radial audit and exact new-angular substitution.

Fourier reconstruction uses the Laguerre three-term recurrence, independently
of the primary verifier's coefficient-sum formula. Pair coupling uses the
signed full interval expression. All output is confined to this audit folder.
"""
from pathlib import Path
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
old=json.loads((D/'inputs/radial46.json').read_text())
new=json.loads((D/'inputs/new_complete_certificate.json').read_text())
ang=json.loads((D/'inputs/angular46.json').read_text())
assert old['dimension']==new['dimension']==ang['dimension']==46
unchanged={k:old[k]==new[k] for k in ['physical_polynomial','fourier_polynomial','M_exact']}
assert all(unchanged.values())
canonical=lambda r:json.dumps({k:r[k] for k in unchanged},sort_keys=True,separators=(',',':')).encode()
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(canonical(old))==sha(canonical(new))
physical=P(list(map(Q,new['physical_polynomial'])))
fourier=P(list(map(Q,new['fourier_polynomial'])))
M=Q(new['M_exact'])
B=Q(ang['B_exact'])
b0=Q(new['b0_exact'])
assert b0==Q(ang['gegenbauer_coefficients'][0])==1 and M>0 and fourier(0)>0
assert physical.degree()==fourier.degree()==111

# Fourier[e^{-pi r^2}(pi r^2)^k]=k! L_k^(22)(pi rho^2)e^{-pi rho^2}.
x=P([0,1])
L0,L1=P([1]),23-x
transformed=physical[0]*L0+physical[1]*L1
factorial=1
for k in range(1,111):
    L2=((2*k+23-x)*L1-(k+22)*L0)/(k+1)
    factorial*=k+1
    transformed+=physical[k+1]*factorial*L2
    L0,L1=L1,L2
assert transformed==fourier
print('RADIAL_HASH_AND_FOURIER_PASS',sha(canonical(old)),flush=True)


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
strict=Q(897267583,10**15)
old_strict=Q(914882847,10**15)
comparison=Q(9295151556,10**16)
upper=objective*(Q(3,2)*pu)**23/math.factorial(23)
lower=objective*(Q(3,2)*pl)**23/math.factorial(23)
assert 0<lower<upper<strict<old_strict<comparison
density=arb(objective)*(Q(3,2)*pi)**23/math.factorial(23)
assert density<arb(strict)
result={'pass':True,'dimension':46,'precision_bits':ctx.prec,
        'unchanged_radial_fields':unchanged,'radial_canonical_sha256':sha(canonical(old)),
        'source_radial_file_sha256':sha((D/'inputs/radial46.json').read_bytes()),
        'new_certificate_file_sha256':sha((D/'inputs/new_complete_certificate.json').read_bytes()),
        'angular_file_sha256':sha((D/'inputs/angular46.json').read_bytes()),
        'Fourier_identity_by_Laguerre_recurrence':True,'root_results':root_results,
        'coupling_intervals':len(leaves),'coupling_depth':max_depth,
        'coupling_interval_union_checked':True,'coupling_full_signed_expression':True,
        'b0_exact':str(b0),'M_exact':str(M),'B_exact':str(B),
        'origin_cost_exact':str(origin),'origin_PSD_rank_one_identity':True,
        'objective_exact':str(objective),'objective_matches_primary_exact':True,
        'damping_lower_exceeds_40_exact':True,
        'strict_density_upper':'0.000000897267583',
        'density_interval':str(density),'rational_rounding_strict':True,
        'pinned_comparator':'0.0000009295151556',
        'fractional_improvement_lower':str(1-strict/comparison),
        'scope':'Radial, exact substitution and rounding. Angular/cyclic result files also required.'}
(D/'radial_substitution_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('RADIAL_SUBSTITUTION_PASS',json.dumps({k:v for k,v in result.items() if k not in ['objective_exact','root_results']}),flush=True)
