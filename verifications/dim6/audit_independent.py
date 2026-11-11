"""Independent rational antipodal, Fourier transform and origin arithmetic audit.
Usage: python audit_independent.py exported_directory
Does NOT certify the full far or compact cyclic domains.
"""
if not __debug__: raise RuntimeError("Verification must not run under python -O.")
import json,sys,math
from pathlib import Path
from fractions import Fraction
def fq(x): return Q(str(Fraction(str(x))))
from flint import fmpq as Q,fmpq_poly as Poly,fmpq_mat as Mat,arb,ctx
ctx.prec=256
path=Path(sys.argv[1]);s=json.load(open(path/'witness.json'));m=json.load(open(path/'witness_psd_exact.json'))['matrices'];pre=json.load(open(path/'prerequisites.json'));n=s['dimension'];p=s['p'];K=s['K'];a=Q(n,2);X=Poly([0,1]);c=[fq(v) for v in s['x']]
# Explicit Laguerre coefficients, independent of the exporter's three-term recurrence.
def binom(x,k):
 r=Q(1)
 for j in range(k):r*= (x-j)/Q(j+1)
 return r
ls=[Poly([(-2)**j*binom(k+a-1,k-j)/Q(math.factorial(j)) for j in range(k+1)]) for k in range(len(c))]
P=sum((c[k]*ls[k]/ls[k][0] for k in range(len(c))),Poly([]));H=sum(((-1)**k*c[k]*ls[k]/ls[k][0] for k in range(len(c))),Poly([]))
# Multiplication by pi|x|^2 transforms to minus the radial Laplacian/(4pi).
monos=[Poly([1])]
for j in range(P.degree()):
 q=monos[-1];monos.append(-X*q.derivative().derivative()+(2*X-a)*q.derivative()+(a-X)*q)
assert sum((P[j]*monos[j] for j in range(P.degree()+1)),Poly([]))==H
assert P[0]==sum(c,Q(0)) and H[0]==sum(((-1)**j*z for j,z in enumerate(c)),Q(0))
M=[Mat([[Q(z) for z in row] for row in block]) for block in m];B=Mat([[M[0][i+1,j+1] for j in range(p+1)] for i in range(p+1)]);b=Mat([[M[0][i+1,0]] for i in range(p+1)]);sat=(b.transpose()*B.solve(b))[0,0];alpha=M[0][0,0]-sat;assert alpha==Q(pre['origin_slack'])
g=pre['guards']['reciprocal_gaussian_guards'];gc=Q(g['c']);eta=Q(g['eta']);zeta=Q(g['zeta']);norm=H[0]+eta*2**n*(gc-2*n)+zeta/Q(2**n)*(gc-Q(n,8));assert norm==Q(g['corrected_normalization'])
obj=(P[0]+sat+sum((fq(w)*fq(v) for w,v in zip(s['count_weights'],s['count_bounds'])),Q(0))+gc*(eta+zeta))/norm;assert obj==Q(pre['saturated_objective'])
# At angle -1 normalized degree-l Gegenbauer = (-1)^l exactly.
T=[Poly([1]),2*X-1]
for k in range(1,p):T.append(2*(2*X-1)*T[-1]-T[-2])
C=Mat([[sum(((-1)**ell*M[ell][i+(ell==0),j+(ell==0)] for ell in range(K+1)),Q(0)) for j in range(p+1)] for i in range(p+1)])

def basis(lo,hi):
 pol=[t(Poly([lo,hi-lo])) for t in T]
 return Mat([[sum((pol[k][j]*Q(math.comb(i,j),math.comb(p,j)) for j in range(i+1)),Q(0)) for k in range(p+1)] for i in range(p+1)])
maximum=None;where=None
for u in range(8):
 for v in range(8):
  co=basis(Q(u,8),Q(u+1,8))*C*basis(Q(v,8),Q(v+1,8)).transpose()
  for i in range(p+1):
   for j in range(p+1):
    if maximum is None or co[i,j]>maximum:maximum=co[i,j];where=(u,v,i,j)
result={'fourier_polynomial_identity':True,'origin_and_gaussian_arithmetic':True,'antipodal_boundary_8x8_bernstein':maximum<=0,'antipodal_max_coefficient':str(maximum),'antipodal_max_coefficient_float':float(maximum),'location':where,'saturated_density':str(arb.pi()**(arb(n)/2)/(arb(n)/2+1).gamma()/arb(2)**n*arb(obj)),'scope':'Only algebra and antipodal boundary; not a complete density certificate'}
print(json.dumps(result,indent=2));(path/'independent_algebra_audit.json').write_text(json.dumps(result,indent=2)+'\n')
