"""Export exact witnesses, without asserting their feasibility."""
import json,sys,math
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as Poly
rows=json.load(open(sys.argv[1]));dest=Path(sys.argv[2]);dest.mkdir(parents=True,exist_ok=True)
for r in rows:
 if 'coefficients' not in r:continue
 n=r['dimension'];D=r['radial_degree'];x=Poly([0,1]);a=Q(n,2);ls=[Poly([1]),a-2*x]
 for k in range(1,D):ls.append(((2*k+a-2*x)*ls[-1]-(k+a-1)*ls[-2])/Q(k+1))
 norms=[Q(1)]
 for k in range(1,D+1):norms.append(norms[-1]*(a+k-1)/k)
 co=[Q(str(F(str(v)))) for v in r['coefficients']];P=sum((co[k]*ls[k]/norms[k] for k in range(D+1)),Poly([]));H=sum(((-1)**k*co[k]*ls[k]/norms[k] for k in range(D+1)),Poly([]))
 mats=[[[str(F(float(v)).limit_denominator(10**22)) for v in row] for row in B]for B in r['matrices']]
 obj=(P[0]+Q(mats[0][0][0]))/H[0]
 out={k:r[k]for k in ['dimension','radial_degree','radial_nodes','angular_degree','exclusion_squared','truncation_squared']};out.update(matrices=mats,physical_polynomial=list(map(str,P.coeffs())),fourier_polynomial=list(map(str,H.coeffs())),normalization=str(H[0]),exact_objective=str(obj),subdivisions_per_radial_cell=8,certified_density_upper=str(r['density_candidate']*1.000001),certification_status='UNVERIFIED candidate; must pass both full-domain verifiers')
 (dest/f'dim{n}.json').write_text(json.dumps(out,indent=2))
