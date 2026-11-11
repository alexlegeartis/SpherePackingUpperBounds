"""Independent restriction, Gegenbauer PSD, origin Gram and rational rounding."""
import json,math,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from flint import arb,fmpq,ctx
assert __debug__
sys.set_int_max_str_digits(0);ctx.prec=384
D=Path(__file__).resolve().parent;S=D/'continuum_three_point'
r=json.loads((S/'adaptive_cyclic40_radial111.json').read_text());a=json.loads((S/'cyclic_angular40_exact.json').read_text());certificate=json.loads((D/'certificate47_transfer.json').read_text())
p=list(map(F,r['physical_polynomial']));h48=list(map(F,r['fourier_polynomial']));pol=list(map(F,a['polynomial_coefficients']));B=F(a['B_exact'])
assert len(p)==len(h48)==112 and len(pol)==20 and sum(pol)==B and B>0
assert F(r['B'])==B and F(str(r['kappa']))==F(17,8)
M=F(r['neighbor_cost'])/B*F(1000001,1000000);assert M>0
for filename,digest in certificate['source_git_blobs'].items():
 data=(S/filename).read_bytes();assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==digest

def rising(z,k):return math.prod((z+i for i in range(k)),start=F(1))
def choose(z,k):return math.prod((z-i for i in range(k)),start=F(1))/math.factorial(k)
def direct_fourier(n):
 alpha=F(n,2)-1
 return [sum((p[k]*math.factorial(k)*(-1)**j*choose(alpha+k,k-j)/math.factorial(j) for k in range(j,len(p))),F(0))for j in range(len(p))]
assert direct_fourier(48)==h48
h47=direct_fourier(47)
# Integrate the missing frequency coordinate: integral (pi*t²)^k e^(-pi*t²)dt=(1/2)_k.
integrated=[sum((h48[k]*math.comb(k,j)*rising(F(1,2),k-j)for k in range(j,len(p))),F(0))for j in range(len(p))]
assert h47==integrated
N=h47[0];assert N>0 and N==sum((p[k]*rising(F(47,2),k)for k in range(len(p))),F(0))
# Explicit normalized Gegenbauer coefficients, not the producer recurrence.
basis=[];lam=F(45,2)
for ell in range(20):
 v=[F(0)]*20;norm=rising(2*lam,ell)/math.factorial(ell)
 for j in range(ell//2+1):v[ell-2*j]=(-1)**j*rising(lam,ell-j)*2**(ell-2*j)/(math.factorial(j)*math.factorial(ell-2*j)*norm)
 assert sum(v)==1;basis.append(v)
remaining=pol[:];coeff=[F(0)]*20
for ell in range(19,-1,-1):
 coeff[ell]=remaining[ell]/basis[ell][ell]
 for j in range(20):remaining[j]-=coeff[ell]*basis[ell][j]
assert not any(remaining) and min(coeff)>=0 and coeff[0]>1 and sum(coeff)==B
b0=coeff[0];moment=sum((pol[2*j]*rising(F(1,2),j)/rising(F(47,2),j)for j in range(10)),F(0));assert moment==b0
# Degree-zero origin/neighbor Gram matrix: [[MB/b0,-M],[-M,M*b0/B]].
origin=M*B/b0;neighbor=M*b0/B;assert origin>0 and neighbor>0 and origin*neighbor==M*M
assert M*B-origin>0
objective=(p[0]+origin)/N;assert objective>0
assert certificate['dimension']==47 and certificate['source_dimension']==48
assert F(certificate['fourier_origin_exact'])==N and F(certificate['angular_constant_exact'])==b0
assert list(map(F,certificate['gegenbauer_coefficients']))==coeff and list(map(F,certificate['fourier_coefficients']))==h47
assert F(certificate['M_exact'])==M and F(certificate['origin_cost_exact'])==origin

def atan_bounds(x):
 lower=sum((-1)**j*x**(2*j+1)/(2*j+1)for j in range(80));return lower,lower+x**161/161
l5,u5=atan_bounds(F(1,5));l239,u239=atan_bounds(F(1,239));pl=16*l5-4*u239;pu=16*u5-4*l239
assert F(51,4)*pl>40
# V47(sqrt(6)/2)^2=4*pi^46*6^47*(24!)²/(48!)².
square_upper=objective**2*4*pu**46*6**47*math.factorial(24)**2/math.factorial(48)**2
threshold=F('0.000000616879891');assert square_upper<threshold**2
q=lambda z:fmpq(z.numerator,z.denominator)
volume=2*arb.pi()**23*arb(6)**23*arb(6).sqrt()*math.factorial(24)/math.factorial(48)
density=arb(q(objective))*volume
unsaturated=(p[0]+M*B)/N;assert unsaturated>0
old=F('0.0000006490757338')
assert unsaturated**2*4*pl**46*6**47*math.factorial(24)**2/math.factorial(48)**2>old**2
unsaturated_density=arb(q(unsaturated))*volume
out={'pass':True,'dimension':47,'Fourier_marginal_and_direct_identity':True,'angular_PSD_coefficients':len(coeff),'angular_constant_exact':str(b0),'Fourier_origin_exact':str(N),'origin_Gram_determinant':'0','objective_exact':str(objective),'density_interval':str(density),'strict_upper':'0.000000616879891','unsaturated_density_interval':str(unsaturated_density),'negative_control_above_pinned_comparator':True,'source_dependency':'All three source48 continuum replays must pass.'}
(D/'independent_algebra_audit.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items()if not k.endswith('_exact')})
