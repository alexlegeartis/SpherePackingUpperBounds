"""Exact algebraic 48D-to-nD restriction and angular-origin saturation.

Requires full source48 replay; run verify_all.py for the complete proof chain.
No numerical optimization, sampled sign test, or new radial root computation.
"""
import json, math, sys, hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq, fmpq_poly, arb, ctx
if not __debug__: raise RuntimeError('Assertions must be enabled')
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).parent;S=D/'continuum_three_point'
EXPECTED={'adaptive_cyclic40_radial111.json':'fac908429ec69fdbbc20b28b1d23fa792822cf0c','cyclic_angular40_exact.json':'33e54a15e6fcc9684847efb97235a0cafd91a457'}
for name,want in EXPECTED.items():
 data=(S/name).read_bytes();assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==want
r=json.loads((S/'adaptive_cyclic40_radial111.json').read_text());q=json.loads((S/'cyclic_angular40_exact.json').read_text())
p=list(map(F,r['physical_polynomial']));h48=list(map(F,r['fourier_polynomial']));Q=fmpq_poly(list(map(fmpq,q['polynomial_coefficients'])));B=F(q['B_exact']);M=F(r['neighbor_cost'])/B*F(1000001,1000000)
assert M>0 and F(r['B'])==B and F(str(r['kappa']))==F(17,8) and len(p)==112
qq=lambda x:fmpq(x.numerator,x.denominator)
ff=lambda x:F(int(x.numerator),int(x.denominator))
aa=lambda x:arb(qq(x))
def rising(a,k):return math.prod((a+j for j in range(k)),start=F(1))
def basis(n,d):
 t=fmpq_poly([0,1]);G=[fmpq_poly([1]),t]
 for j in range(1,d):G.append(((2*j+n-2)*t*G[-1]-j*G[-2])/(j+n-2))
 assert all(g(1)==1 for g in G)
 return G
results=[]
for n in range(25,49):
 # Triangular conversion reconstructs every exact angular coefficient.
 G=basis(n,Q.degree());rem=Q;cs=[F(0)]*len(G)
 for j in range(len(G)-1,-1,-1):
  cs[j]=ff(rem[j]/G[j][j]);rem-=qq(cs[j])*G[j]
 assert not rem and min(cs)>=0 and cs[0]>0 and sum(cs)==B
 b0=cs[0];moment=sum((ff(Q[2*j])*rising(F(1,2),j)/rising(F(n,2),j) for j in range(Q.degree()//2+1)),F(0));assert b0==moment
 N=sum((pk*rising(F(n,2),k) for k,pk in enumerate(p)),F(0));assert N>0
 if n==47:
  # Independent full Fourier transformation and codimension-one integral identity.
  hn=[sum((p[k]*math.factorial(k)*(-1)**j*rising(F(n,2)+j,k-j)/math.factorial(k-j)/math.factorial(j) for k in range(j,len(p))),F(0)) for j in range(len(p))]
  hi=[sum((h48[k]*math.comb(k,j)*rising(F(48-n,2),k-j) for k in range(j,len(p))),F(0)) for j in range(len(p))]
  assert hn==hi and hn[0]==N
  # Origin Schur saturation is an exact rank-one Gram identity.
  origin=M*B/b0;cross=-M;diag=M*b0/B
  assert origin>0 and diag>0 and origin*diag-cross*cross==0
  assert M*B-origin>0
  (D/'certificate47_transfer.json').write_text(json.dumps({'dimension':n,'source_dimension':48,'exclusion_squared':'6','kappa':'17/8','source_git_blobs':EXPECTED,'fourier_origin_exact':str(N),'angular_constant_exact':str(b0),'gegenbauer_coefficients':[str(x) for x in cs],'M_exact':str(M),'origin_cost_exact':str(origin),'fourier_coefficients':[str(x) for x in hn]},indent=2))
 obj=(p[0]+M*B/b0)/N
 a=arb(n)/2;density=aa(obj)*(arb.pi()*arb(3)/2)**a/(a+1).gamma()
 row={'dimension':n,'fourier_origin_exact':str(N),'angular_constant_exact':str(b0),'objective_exact':str(obj),'density_interval':str(density),'density_approx':float(density.mid()),'angular_PSD_verified':True}
 if n==47:
  upper=F('0.000000616879891');assert density<aa(upper)
  direct=aa((p[0]+M*B)/N)*(arb.pi()*arb(3)/2)**a/(a+1).gamma()
  assert direct>aa(F('0.0000006490757338'))
  row['unsaturated_restriction_density_interval']=str(direct)
  row['origin_cost_decrease_exact']=str(M*B-M*B/b0)
  baseline=F('0.0000006490757338');assert upper<baseline
  row.update(strict_upper=str(upper),strict_upper_decimal='0.000000616879891',pinned_comparator=str(baseline),fractional_improvement_lower=str(1-upper/baseline),full_Fourier_identity_verified=True,origin_rank_one_PSD_verified=True)
 results.append(row)
 print(n,density,'b0',float(b0),flush=True)
(D/'transfer_exact_audit.json').write_text(json.dumps({'pass':True,'scope':'Exact algebra; full source48 replay is required separately.','precision_bits':ctx.prec,'source_commit':'e756587e8fdf7827f085d6d0372cac267b110bd1','results':results},indent=2))
print('TRANSFER_ALGEBRA_PASS')
