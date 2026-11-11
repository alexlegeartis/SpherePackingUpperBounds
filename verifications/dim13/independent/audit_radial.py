"""Independent exact export, monomial Fourier transform, rays and curved pair.

Export uses explicit generalized Laguerre coefficients, not producer recurrence.
The pair proof uses the ratio gap and a linear-interpolation curvature bound,
not the producer's direct-function midpoint Taylor test.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,sys,time,gzip
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Do not use Python -O')
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).resolve().parent
r=json.loads((D/'inputs/certificate.json').read_text())
a=json.loads((D/'inputs/angular_exact.json').read_text())
candidate=json.loads((D/'inputs/discovery_candidate.json').read_text())
n=13
assert r['dimension']==a['dimension']==candidate['dimension']==n
assert r['degree']==candidate['degree']==63
s=Q(r['s']);aa=Q(r['a_exact']);bb=Q(r['b_exact']);cut=Q(r['cut_X'])
assert (s,aa,bb,cut)==(Q(13,8),Q(10),Q(16),Q(12))
assert (Q(candidate['s_exact']),Q(candidate['a_exact']),Q(candidate['b_exact']),Q(str(Fraction(candidate['cut_X']))))==(s,aa,bb,cut)
assert Q(a['a_exact'])==aa and Q(a['b_exact'])==bb
physical=P(list(map(Q,r['physical'])));fourier=P(list(map(Q,r['fourier'])))
assert len(r['physical'])==len(r['fourier'])==64 and physical.degree()==fourier.degree()==63
B=Q(r['B_exact']);M=Q(r['M_exact']);assert B==Q(a['B_exact']) and B>0 and M>0

# Closed coefficient formula, distinct from the producer's export recurrence.
v=[Q(str(Fraction(float(z)).limit_denominator(10**16)))for z in candidate['radial']]
C=Q(str(Fraction(float(candidate['C'])).limit_denominator(10**16)));assert len(v)==64 and C>0
alpha=Q(n,2);rise=[Q(1)]
for k in range(63):rise.append(rise[-1]*(alpha+k))
E=[]
for k in range(64):
    coeff=[Q((-2)**j)*rise[k]/rise[j]/math.factorial(k-j)/math.factorial(j)/(rise[k]/math.factorial(k))for j in range(k+1)]
    E.append(P(coeff));assert E[-1](0)==1
assert sum((c*e for c,e in zip(v,E)),P([]))==physical
assert sum(((-1)**k*c*e for k,(c,e)in enumerate(zip(v,E))),P([]))==fourier
assert M==C/B*Q(1000001,1000000)
assert Q(r['origin_cost_exact'])==M*B and Q(r['normalization_exact'])==fourier(0)>0
print('EXACT_EXPORT_AND_NORMALIZATION_PASS',flush=True)

# Fourier[X^k exp(-X)] = k! L_k^(alpha-1)(Y) exp(-Y).
X=P([0,1]);L0=P([1]);L1=alpha-X
transformed=physical[0]*L0+physical[1]*L1;factorial=1
for k in range(1,63):
    L2=((2*k+alpha-X)*L1-(k+alpha-1)*L0)/(k+1)
    factorial*=k+1;transformed+=physical[k+1]*factorial*L2;L0,L1=L1,L2
assert transformed==fourier
print('INDEPENDENT_MONOMIAL_FOURIER_TRANSFORM_PASS',flush=True)

# Check the precise angular data and PSD linkage; cyclic geometry is separate.
G=[P([1]),X]
cs=list(map(Q,a['gegenbauer_coefficients']));angular=P(list(map(Q,a['polynomial_coefficients'])))
for ell in range(1,len(cs)-1):G.append(((2*ell+n-2)*X*G[-1]-ell*G[-2])/(ell+n-2))
assert cs[0]==1 and all(c>=0 for c in cs)
assert all(g(1)==1 for g in G) and sum((c*g for c,g in zip(cs,G)),P([]))==angular
assert sum(cs)==angular(1)==B
moment=Q(1);b0=Q(0)
for j in range(angular.degree()//2+1):
    if j:moment*=Q(2*j-1,n+2*j-2)
    b0+=angular[2*j]*moment
assert b0==1
origin=M*B;degree_zero=M/B;assert origin*degree_zero==M*M and origin>0 and degree_zero>0
print('ANGULAR_IDENTITY_AND_ORIGIN_PSD_LINKAGE_PASS',flush=True)

# Rational Sturm signs at the exact finite endpoint and positive infinity.
# Positive primitive rescaling does not change a Sturm sign sequence.
def sturm(poly):
 seq=[poly,poly.derivative()]
 while seq[-1].degree()>0:
  rem=-(seq[-2]%seq[-1]);assert rem
  num=rem.numer();seq.append(P(num/abs(num.content())))
 return seq
def variation(vals):
 signs=[int(v>0)-int(v<0)for v in vals if v]
 return sum(a!=b for a,b in zip(signs,signs[1:]))
roots_result={}
for name,poly,endpoint in [('fourier',fourier,Q(0)),('physical',physical,cut)]:
 seq=sturm(poly);v0=variation([p(endpoint)for p in seq]);vinf=variation([p[p.degree()]for p in seq])
 assert v0==vinf,(name,v0,vinf)
 roots_result[name]={'method':'Primitive rational Sturm','variation_at_endpoint':v0,'variation_at_infinity':vinf,'cutoff_exact':str(endpoint),'sequence_degrees':[p.degree()for p in seq]}
 print(name.upper()+'_RATIONAL_STURM_RAY_PASS',v0,vinf,flush=True)
assert fourier(0)>0 and physical(cut)<0

# Pair ratio R=f/h-M(2-h). h>0, so its sign is the original pair gap.
polys=[arb_poly([arb(x)for x in p.coeffs()])for p in [physical,physical.derivative(),physical.derivative().derivative()]]
pi=arb.pi();denom=pi*arb(s)
assert denom>0 and arb(cut)>denom

def tight(poly,x):
    cc=poly(arb_poly([x.mid(),x.rad()]));err=sum((abs(t)for t in cc.coeffs()[1:]),arb(0))
    return cc[0]+arb(0,err.upper())

def ratio(x,curvature=False):
    z=x/denom-1;phi=arb(aa)*z+arb(bb)*z*z
    k=(arb(aa)+2*arb(bb)*z)/denom;l=2*arb(bb)/(denom*denom)
    e=(phi-x).exp();h=(-phi).exp()
    pp=tight(polys[0],x)
    if not curvature:return pp*e-arb(M)*(2-h)
    dp=tight(polys[1],x);ddp=tight(polys[2],x)
    return (ddp+2*(k-1)*dp+((k-1)**2+l)*pp)*e+arb(M)*(k*k-l)*h

stack=[(Q(0),Q(1),0)];leaves=[];maxdepth=0;start=time.time()
while stack:
    lo,hi,depth=stack.pop()
    xl=denom+(arb(cut)-denom)*arb(lo);xh=denom+(arb(cut)-denom)*arb(hi)
    whole=denom+(arb(cut)-denom)*arb((lo+hi)/2,(hi-lo)/2)
    left,right=ratio(xl),ratio(xh)
    assert left.is_finite() and right.is_finite()
    assert not left>0 and not right>0,('POSITIVE_PAIR_ENDPOINT',str(lo),str(hi),str(left),str(right))
    curvature=ratio(whole,True);assert curvature.is_finite()
    # For a twice differentiable F, F lies below the endpoint secant plus
    # max(0,-inf F'')*(width)^2/8 throughout the interval.
    width=xh-xl;assert width>0
    upper=max(left.upper(),right.upper())+max(arb(0),-curvature.lower())*width**2/8
    if upper<0:
        leaves.append({'parameter_lo':str(lo),'parameter_hi':str(hi),'depth':depth,'ratio_gap_upper':str(upper.upper())});maxdepth=max(maxdepth,depth);continue
    assert depth<36,('UNRESOLVED_PAIR',str(lo),str(hi),str(upper))
    mid=(lo+hi)/2;stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
ordered=sorted(leaves,key=lambda x:Q(x['parameter_lo']))
assert Q(ordered[0]['parameter_lo'])==0 and Q(ordered[-1]['parameter_hi'])==1
assert all(Q(a['parameter_hi'])==Q(b['parameter_lo'])for a,b in zip(ordered,ordered[1:]))
(D/'pair_cover.json').write_text(json.dumps(ordered,indent=2)+'\n')
print('COMPLETE_CURVED_PAIR_PASS',len(leaves),'depth',maxdepth,'seconds',time.time()-start,flush=True)

# Independent Machin rational pi enclosure and correct integer-dimension volume.
def atan_bounds(t):
    lo=sum(((-1)**k*t**(2*k+1)/(2*k+1)for k in range(128)),Q(0));return lo,lo+t**257/257
al,au=atan_bounds(Q(1,5));bl,bu=atan_bounds(Q(1,239));pl,pu=16*al-4*bu,16*au-4*bl
assert 0<pl<pu and arb(pl).union(arb(pu)).overlaps(pi)
objective=(physical(0)+origin)/fourier(0);assert objective>0
strict=Q(str(Fraction('0.059905026')));assert Q(str(Fraction(r['strict_density_upper'])))==strict
k=n//2
if n%2:
    odd_factorial=math.prod(range(1,n+1,2))
    vol2=lambda p:p**(2*k)*s**n/(2**k*odd_factorial)**2
else:
    vol2=lambda p:(p**k*s**k/(4**k*math.factorial(k)))**2
lower2=objective**2*vol2(pl);upper2=objective**2*vol2(pu)
assert 0<lower2<upper2<strict*strict
volume=arb(vol2(pl)).union(arb(vol2(pu))).sqrt()
density=arb(objective)*volume;assert density<arb(strict)
comparator=Q(str(Fraction('0.06012103')));assert strict<comparator
result={'pass':True,'dimension':n,'precision_bits':ctx.prec,'degree':63,
        'export_exact':True,'export_rationalization_denominator_limit':10**16,'M_guard_exact':'1000001/1000000',
        'Fourier_identity_method':'Ordinary Laguerre recurrence for Fourier transforms of every monomial X^k exp(-X).',
        'Fourier_identity_exact':True,'positive_normalization_exact':str(fourier(0)),
        'normalization_interval':str(arb(fourier(0))),
        'B_exact':str(B),'M_exact':str(M),'origin_cost_exact':str(origin),'objective_exact':str(objective),
        's_exact':str(s),'profile_a_exact':str(aa),'profile_b_exact':str(bb),'physical_cut_X_exact':str(cut),
        'root_results':roots_result,'pair_interval':f'[pi*({s}),12]','pair_intervals':len(leaves),'pair_depth':maxdepth,
        'pair_parameter_cover_complete':True,'pair_method':'Ratio gap, endpoint secant and rigorous second-derivative lower bound.',
        'angular_polynomial_degree':angular.degree(),'stored_harmonic_degree_limit':len(cs)-1,
        'angular_expansion_and_PSD_linkage_checked':True,'origin_rank_one_PSD_exact':True,
        'full_angular_cyclic_geometry_replayed_here':False,'angular_dependency':'Separate independent curved-profile angular/cyclic audit required.',
        'strict_density_upper':str(strict),'density_interval':str(density),'rational_squared_rounding':True,
        'volume_formula':'Odd: pi^k s^(n/2)/(2^k n!!); even: pi^k s^k/(4^k k!)',
        'density_squared_rational_upper':str(upper2),
        'pinned_comparator':str(comparator),'fractional_improvement_at_least':str(1-strict/comparator),
        'scope':'All radial/normalization/PSD-linkage/rounding gates pass. Density conclusion additionally requires the linked full angular/cyclic proof.'}
(D/'radial_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('RADIAL_AUDIT_PASS',json.dumps({k:v for k,v in result.items()if k not in ['positive_normalization_exact','objective_exact','root_results']}),flush=True)
