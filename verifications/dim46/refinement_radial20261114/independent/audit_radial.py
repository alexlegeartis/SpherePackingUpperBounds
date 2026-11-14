"""Independent recurrence Fourier reconstruction and signed interval coupling.

Reads only frozen exact witness data. Root isolation shares FLINT/Arb with the
primary check; no claim of a separate arithmetic backend is made.
"""
import json, math, sys, time, hashlib
from fractions import Fraction
from pathlib import Path
from flint import fmpq as Q, fmpq_poly as P, arb, arb_poly, ctx

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')
sys.set_int_max_str_digits(0);ctx.prec=768
D=Path(__file__).resolve().parent
r=json.loads((D/'inputs/certificate.json').read_text())
cfg=json.loads((D/'inputs/config.json').read_text())
ang=json.loads((D/'inputs/angular.json').read_text())
n=cfg['dimension'];assert n in [45,46,47] and n==r['dimension']==ang['dimension']
p=P(list(map(Q,r['physical_polynomial'])))
h=P(list(map(Q,r['fourier_polynomial'])))
M=Q(r['M_exact']);B=Q(ang['B_exact']);b0=Q(r['b0_exact'])
assert M>0 and b0==Q(ang['gegenbauer_coefficients'][0])==1
assert B==Q(r['B_exact']) and p.degree()==h.degree()==111
assert Q(r['coupled_upper_X_exact'])==Q(r['physical_tail_X_exact'])==40

# Fourier[e^-X X^k] = e^-Y k! L_k^(n/2-1)(Y).
a=Q(n,2);x=P([0,1]);L0=P([1]);L1=a-x
transformed=p[0]*L0+p[1]*L1
factorial=1
for k in range(1,111):
    L2=((2*k+a-x)*L1-(k+a-1)*L0)/(k+1)
    factorial*=k+1
    transformed+=p[k+1]*factorial*L2
    L0,L1=L1,L2
assert transformed==h and h(0)>0 and p(40)<0

def atan(x):
    lo=sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(128)),Q(0))
    return lo,lo+x**257/257
al,au=atan(Q(1,5));bl,bu=atan(Q(1,239));pl,pu=16*al-4*bu,16*au-4*bl
assert Q(51,4)*pl>40
root_results={}
for name,poly,cut in [('fourier',h,Q(0)),('physical',p,Q(40))]:
    roots=poly.numer().complex_roots()
    assert sum(m for z,m in roots)==111
    assert not [z for z,m in roots if z.imag.contains(0) and not z.real<arb(cut)]
    root_results[name]={'all_111_roots_accounted_for':True,'forbidden_real_roots':0,'cutoff':str(cut)}
    print('INDEPENDENT',name,'RAY PASS',flush=True)

ap=arb_poly([arb(c) for c in p.coeffs()]);pi=arb.pi()
stack=[(6*pl,Q(40),0)];leaves=[];maxdepth=0;start=time.time()
while stack:
    lo,hi,depth=stack.pop();mid=(lo+hi)/2;rad=(hi-lo)/2
    terms=ap(arb_poly([arb(mid),arb(rad)]))
    error=sum((abs(c) for c in terms.coeffs()[1:]),arb(0))
    pbox=terms[0]+arb(0,error.upper());xx=arb(mid,rad)
    # Full signed expression, after dividing by the positive h(X).
    gap=pbox*(Q(9,8)*xx-Q(51,4)*pi).exp()-arb(M)*(2-(-Q(17,8)*(xx-6*pi)).exp())
    assert gap.is_finite()
    if gap<0:
        leaves.append((lo,hi));maxdepth=max(maxdepth,depth);continue
    assert depth<36,(str(lo),str(hi),str(gap))
    stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
leaves.sort()
assert leaves[0][0]==6*pl and leaves[-1][1]==40
assert all(a[1]==b[0] for a,b in zip(leaves,leaves[1:]))
partition='\n'.join(f'{lo} {hi}' for lo,hi in leaves).encode()
origin=M*B/b0;nonorigin=M*b0/B
assert origin>0 and nonorigin>0 and origin*nonorigin==M*M
obj=(p(0)+origin)/h(0);assert obj==Q(r['objective_exact'])>0
strict=Q(str(Fraction(cfg['strict_density_upper_bound'])))
old=Q(str(Fraction(cfg['previous_repository_strict_upper'])))
assert 0<strict<old
if n%2==0:
    upper=obj*(Q(3,2)*pu)**(n//2)/math.factorial(n//2)
    assert upper<strict
else:
    m=(n-1)//2
    volume2=4*pu**(n-1)*6**n*Q(math.factorial(m+1),math.factorial(n+1))**2
    assert 0<obj**2*volume2<strict**2
density=arb(obj)*(Q(3,2)*pi)**(arb(n)/2)/(arb(n)/2+1).gamma()
assert density<arb(strict)
out={'pass':True,'dimension':n,'precision_bits':ctx.prec,
     'exact_Fourier_identity_by_recurrence':True,'root_results':root_results,
     'signed_full_pair_expression':True,'pair_interval_union_checked':True,
     'pair_intervals':len(leaves),'pair_depth':maxdepth,'pair_seconds':time.time()-start,
     'pair_partition_sha256':hashlib.sha256(partition).hexdigest(),
     'origin_rank_one_PSD_identity':True,'normalization_positive':True,
     'objective_exact':str(obj),'rational_Machin_rounding':True,
     'strict_density_upper_bound':cfg['strict_density_upper_bound'],'density_interval':str(density),
     'previous_repository_strict_upper':cfg['previous_repository_strict_upper'],
     'fractional_improvement_lower':str(1-strict/old),
     'scope':'Independent radial audit; independent angular/cyclic audit also required'}
(D/'radial_result.json').write_text(json.dumps(out,indent=2)+'\n')
print('INDEPENDENT RADIAL PASS',n,str(density),len(leaves),maxdepth,flush=True)
