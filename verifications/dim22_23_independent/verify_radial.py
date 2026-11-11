"""Independent 22D/23D audit: Laguerre recurrence, root signs, Taylor coupling.

No producer module is imported. Run with argument 22 or 23.
Optional --sturm requests slow exact rational Sturm instead of certified roots.
"""
import hashlib, json, math, sys, time
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq, fmpq_poly, arb, arb_poly, ctx

if not __debug__:
    raise RuntimeError('Assertions must be enabled')
sys.set_int_max_str_digits(0)
ctx.prec = 512
D = Path(__file__).resolve().parent
n = int(sys.argv[1]); assert n in (22, 23)
hashes = json.loads((D/'source_hashes.json').read_text())
for name in (f'certificate{n}.json', f'angular{n}.json'):
    assert hashlib.sha256((D/name).read_bytes()).hexdigest() == hashes[name]['sha256']
d = json.loads((D/f'certificate{n}.json').read_text())
ang = json.loads((D/f'angular{n}.json').read_text())
assert d['dimension'] == ang['dimension'] == n
q = lambda x: fmpq(x.numerator, x.denominator)
a = lambda x: arb(q(x))
frac = lambda x: F(int(x.numerator), int(x.denominator))
p = list(map(F, d['physical_polynomial']))
h = list(map(F, d['fourier_polynomial']))
s = F(d['r2']); kappa = F(d['kappa']); M = F(d['M_exact'])
assert s == F(n,8) and kappa == F(102,n) and F(d['R2']) == 7
assert s*kappa == F(51,4) and M>0 and h[0]>0
P = fmpq_poly(list(map(q,p))); H = fmpq_poly(list(map(q,h)))
assert len(p) == len(h) == 112 and P.degree() == H.degree() == 111

# The generating-function identity for the Fourier convention exp(-2*pi*i*x.y)
# gives FT[X^k exp(-X)] = k! L_k^(n/2-1)(Y) exp(-Y).
# Build these polynomials by their three-term recurrence, not the producer's
# coefficientwise rising-factorial transform.
Y = fmpq_poly([0,1]); halfn = fmpq(n,2)
L = [fmpq_poly([1]), halfn-Y]
for k in range(1,111):
    L.append((2*k+halfn-Y)*L[-1]-k*(k+halfn-1)*L[-2])
assert sum((q(p[k])*L[k] for k in range(112)), fmpq_poly()) == H
assert sum((q(h[k])*L[k] for k in range(112)), fmpq_poly()) == P
print(n, 'independent Laguerre Fourier recurrence PASS', flush=True)

# Explicit Gegenbauer sum, independent of the recurrence used by verify.py.
def rise(x,k):
    return math.prod((x+i for i in range(k)), start=F(1))
poly = list(map(F,ang['polynomial_coefficients']))
assert len(poly)==20 and ang['degree']==14 and not any(poly[15:])
cs = list(map(F,ang['gegenbauer_coefficients']))
assert len(cs)==20 and min(cs)>=0 and cs[0]==1
reconstructed = [F(0)]*20
for ell, coeff in enumerate(cs):
    norm = rise(F(n-2),ell)/math.factorial(ell)
    for j in range(ell//2+1):
        reconstructed[ell-2*j] += coeff*(-1)**j*rise(F(n-2,2),ell-j)*2**(ell-2*j)/(math.factorial(j)*math.factorial(ell-2*j)*norm)
assert reconstructed == poly
B=F(ang['B_exact']); assert B==sum(poly)==sum(cs)
b0=sum((poly[2*j]*rise(F(1,2),j)/rise(F(n,2),j) for j in range(10)),F(0))
assert b0==cs[0]==F(d['b0_exact'])==1
origin=M*B/b0; degree0=M*b0/B
assert origin>0 and degree0>0 and origin*degree0==M*M

def atan_bounds(x):
    lower=sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(100)),F(0))
    return lower,lower+x**201/201
l5,u5=atan_bounds(F(1,5)); l239,u239=atan_bounds(F(1,239))
pl=16*l5-4*u239; pu=16*u5-4*l239
assert F(51,4)*pl>40 and 3<pl<pu<4

# Positive primitive normalization preserves Sturm signs. Endpoints and +infty
# suffice to rule out roots on the entire relevant half line.
def sturm(poly):
    seq=[poly,poly.derivative()]
    while seq[-1].degree()>0:
        rem=-(seq[-2]%seq[-1]); assert rem
        num=rem.numer(); seq.append(fmpq_poly(num/abs(num.content())))
    return seq
def variation(seq,point):
    values=[v[v.degree()] if point is None else v(q(point)) for v in seq]
    signs=[int(v>0)-int(v<0) for v in values if v]
    return sum(x!=y for x,y in zip(signs,signs[1:]))
sign_records={}
for label,pol,cut,sign in [('Fourier',H,F(0),1),('physical',P,7*pl,-1)]:
    t=time.monotonic()
    if '--sturm' in sys.argv:
        print(n,'exact Sturm start',label,flush=True)
        seq=sturm(pol); left=variation(seq,cut); right=variation(seq,None)
        assert left==right and sign*pol(q(cut))>0
        sign_records[label]={'method':'exact Sturm','left_variations':left,'infinity_variations':right,'sequence_length':len(seq)}
    else:
        print(n,'certified root isolation start',label,flush=True)
        roots=pol.numer().complex_roots()
        assert sum(m for z,m in roots)==pol.degree()==111
        assert all(not z.imag.contains(0) or z.real<a(cut) for z,m in roots)
        assert sign*pol(q(cut))>0
        sign_records[label]={'method':'complete FLINT complex root isolation','total_multiplicity':111,'root_disks':len(roots)}
    print(n,'global signs PASS',label,sign_records[label],'seconds',time.monotonic()-t,flush=True)

# Full divided pair function G(X)=P(X)e^((kappa-1)X-C)-2M+M e^(-kappa X+C).
# A second-order Taylor theorem enclosure bounds G on every whole interval.
# This is different from the producer's P-only enclosure and monotonic factors.
AP=arb_poly([a(x) for x in p]); k=a(kappa); beta=k-1; C=51*arb.pi()/4; AM=a(M)
firstpoly=AP.derivative()+beta*AP
secondpoly=AP.derivative().derivative()+2*beta*AP.derivative()+beta**2*AP
stack=[(s*pl,7*pu,0)]; count=depthmax=0; digest=hashlib.sha256()
while stack:
    lo,hi,depth=stack.pop(); mid=(lo+hi)/2; radius=(hi-lo)/2
    v=a(mid); r=a(radius); ea=(beta*v-C).exp(); eb=(-k*v+C).exp()
    center=AP(v)*ea-2*AM+AM*eb
    derivative=firstpoly(v)*ea-k*AM*eb
    taylor=secondpoly(arb_poly([v,r]))
    magnitude=sum((abs(x) for x in taylor.coeffs()),arb(0))
    second=magnitude*(beta*a(hi)-C).exp()+k*k*AM*(-k*a(lo)+C).exp()
    upper=center+abs(derivative)*r+second*r*r/2
    if upper<0:
        count+=1; depthmax=max(depthmax,depth)
        digest.update(f'{lo}|{hi}\n'.encode()); continue
    assert depth<35, (n,depth)
    stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
objective=(p[0]+origin)/h[0]
assert objective==F(d['objective_exact']) and objective>0
thresholds={22:'0.003653383161',23:'0.002642540089'}
threshold=F(thresholds[n])
if n%2==0:
    m=n//2
    rational_upper=objective*(pu*s/4)**m/math.factorial(m)
    assert rational_upper<threshold
else:
    # Ball volume = 2*pi^m*s^(m+1/2)*(m+1)!/(2m+2)!.
    m=(n-1)//2
    squared_upper=objective**2*4*pu**(2*m)*s**(2*m+1)*math.factorial(m+1)**2/math.factorial(2*m+2)**2
    assert squared_upper<threshold**2
density=a(objective)*(arb.pi()*a(s)/4)**(arb(n)/2)/(arb(n)/2+1).gamma()
assert density<a(threshold)
out={'pass':True,'dimension':n,'Fourier_identity_bidirectional_Laguerre_recurrence':True,
     'global_sign_proof':sign_records,'actual_angular_degree':14,
     'harmonic_PSD_and_origin_Gram':True,'pair_method':'second-derivative Taylor bound',
     'pair_intervals':count,'pair_max_depth':depthmax,'pair_partition_sha256':digest.hexdigest(),
     'exact_outward_rounding':True,'strict_upper':thresholds[n],'density_interval':str(density)}
(D/f'radial{n}.json').write_text(json.dumps(out,indent=2)+'\n')
print(out,flush=True)
