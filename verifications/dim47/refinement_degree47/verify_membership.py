"""Exact degree47 span/export membership, separate from the density proof."""
from pathlib import Path
from fractions import Fraction as F
import json,math,sys
sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;r=json.loads((D/'certificate.json').read_text());n=r['dimension'];a=F(n,2);disc=D/'discovery';candidate=json.loads((disc/f'radial_perturb{n}_degree47.json').read_text());seed=json.loads((disc/'inputs/adaptive_cyclic40_radial111.json').read_text());angular=json.loads((disc/f'cyclic_angular{n}_exact.json').read_text());x=[F(float(v)).limit_denominator(10**14)for v in candidate['x']];assert len(x)==50 and x[0]>=0 and x[-1]>0
rise=lambda b,k:math.prod((b+j for j in range(k)),start=F(1));N48=F(seed['fourier_polynomial'][0]);ps=[F(v)/N48 for v in seed['physical_polynomial']];hs=[F(v)/N48 for v in seed['fourier_polynomial']];hn=[sum((hs[k]*math.comb(k,j)*rise(F(48-n,2),k-j)for k in range(j,len(hs))),F(0))for j in range(len(hs))];p=[x[0]*z for z in ps];h=[x[0]*z for z in hn]
for ell in range(48):
 norm=rise(a,ell)/math.factorial(ell)
 for j in range(ell+1):
  coeff=x[ell+1]*(-2)**j*rise(a+j,ell-j)/math.factorial(ell-j)/math.factorial(j)/norm;p[j]+=coeff;h[j]+=(-1)**ell*coeff
assert p==list(map(F,r['physical_polynomial'])) and h==list(map(F,r['fourier_polynomial'])) and h[0]>0
cs=list(map(F,angular['polynomial_coefficients']));b0=sum((cs[2*j]*rise(F(1,2),j)/rise(a,j)for j in range(len(cs)//2)),F(0));B=F(angular['B_exact'])/b0;M=x[-1]/B*F(1000001,1000000);assert M==F(r['M_exact']) and b0==F(r['b0_exact'])
assert (p[0]+B*M)/h[0]==F(r['objective_exact'])
added=[ell for ell in range(32,48)if x[ell+1]];assert added
out={'pass':True,'dimension':n,'source48_degree':111,'correction_degree':47,'source_t_nonnegative':True,'exact_physical_export_match':True,'exact_Fourier_span_match':True,'exact_M_export_match':True,'nonzero_added_Laguerre_degrees':added,'degree31_span_embedding':'Pad sixteen zero coefficients; basis polynomials do not depend on the maximum included degree.','outside_degree31_span':True,'reason':'The source has degree111 and every normalized Laguerre polynomial has its distinct degree. Nonzero degrees32..47 cannot be absorbed by changing source t or lower coefficients.'};(D/'membership_result.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
