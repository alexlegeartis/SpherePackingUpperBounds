if not __debug__: raise RuntimeError('Run without -O.')
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import numpy as np,json,hashlib,math,sys
from flint import fmpz_mat,fmpq_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
P=Path(__file__).parent
source=(P/'verify.py').read_text();exec(source[source.index('def parent():'):source.index('def rogers_comparison_intervals():')])
G,pr=parent();file=P/sys.argv[1];a=json.load(open(file));B=fmpz_mat(a['normals']);S=fmpz_mat(a['labels']);k=B.nrows();n=92-k
H,T=B.transpose().hnf(transform=True);assert H.tolist()==[[int(i==j) for j in range(k)] for i in range(92)] and abs(T.det())==1
N=fmpz_mat(T.tolist()[k:]);lifts=fmpz_mat(T.tolist()[:k]);assert N*B.transpose()==fmpz_mat(n,k);assert lifts*B.transpose()==fmpz_mat([[int(i==j) for j in range(k)] for i in range(k)])
K=B*G.inv()*B.transpose();period=N*G*N.transpose();assert period.det()==G.det()*K.det();be=F(0)
shifts=fmpq_mat(S*lifts)-fmpq_mat(S)*K.inv()*fmpq_mat(B)*G.inv();assert shifts*fmpq_mat(B.transpose())==fmpq_mat(S.nrows(),k)
coordinates=shifts*fmpq_mat(G*N.transpose())*period.inv();assert coordinates*fmpq_mat(N)==shifts
for i in range(S.nrows()):
 for j in range(i):
  d=fmpq_mat([[int(S[i,t]-S[j,t]) for t in range(k)]]);be=max(be,F(str((d*K.inv()*d.transpose())[0,0])))
assert be<16 and str(be)==a['beta'];assert len(set(tuple(r) for r in S.tolist()))==S.nrows()
d2=F(S.nrows()**2)*((16-be)/4)**n/F(int(period.det()));assert str(d2)==a['density_squared']
baseline_file='published_adaptive_results.json' if n==85 else ('published_mw92_results.json' if n==90 else 'published_simplex_results.json')
base=next(F(r['density_squared']) for r in json.load(open(P/baseline_file))['rows'] if r['dimension']==n)
def atan(q,N):
 s=sum((F((-1)**i,(2*i+1)*q**(2*i+1)) for i in range(N)),F(0));e=F((-1)**N,(2*N+1)*q**(2*N+1));return min(s,s+e),max(s,s+e)
l5,h5=atan(5,100);l239,h239=atan(239,40);pl,ph=16*l5-4*h239,16*h5-4*l239
vol=pl**(n//2)/math.factorial(n//2) if n%2==0 else F(2**(n//2+1),math.prod(range(1,n+1,2)))*pl**(n//2)
ordinary2=d2*vol**2
with localcontext() as c:
 c.prec=100;cd=(Decimal(d2.numerator)/d2.denominator).sqrt();dd=(Decimal(ordinary2.numerator)/ordinary2.denominator).sqrt();floor=F(int(dd*10**90),10**90);assert floor**2<=ordinary2;gain=100*((Decimal((d2/base).numerator)/(d2/base).denominator).sqrt()-1)
 report={'dimension':n,'source_parent_checks':pr,'normal_determinant':str(K.det()),'period_determinant':str(period.det()),'beta':str(be),'density_squared':str(d2),'center_density':str(cd),'gain_percent_vs_current_published_simplex':str(gain),'ordinary_density_lower_rational':str(floor),'ordinary_density_lower_decimal':str(Decimal(floor.numerator)/floor.denominator),'sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'all_period_lift_shift_checks':True,'all_cosets_distinct_by_parent_minimum':True}
print(json.dumps(report,indent=2));(P/(file.stem+'_verified.json')).write_text(json.dumps(report,indent=2))
