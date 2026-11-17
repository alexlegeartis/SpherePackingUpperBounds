"""Independent exact algebra, provenance, scalar correction and normalization."""
from pathlib import Path
import sys,json,math,itertools,hashlib
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,ctx
import sympy as sy
sys.set_int_max_str_digits(0);ctx.prec=512
D=Path(__file__).parent;parent=D.parent;K=D/'geometry/dim51/frozen.json';raw=K.read_bytes();ker=json.loads(raw);assert ker['dimension']==51 and ker['R']=='19/13' and ker['powers']==[6,10,14]
assert raw==(D/'scalar/kernel.json').read_bytes();co=list(map(Q,ker['co']));q=Q(ker['q']);assert sum(co)==q and co[0]==1 and min(co)>0
B=[sy.Matrix(b)for b in ker['B']];z=sy.Matrix(ker['z']);ss=sy.Rational(ker['s']);assert len(B)==len(co)==30
for b in B:
 assert b==b.T and b.shape==(3,3)
 for size in [1,2,3]:
  for ix in itertools.combinations(range(3),size):assert b.extract(ix,ix).det()>0
schur=(z.T*B[0].inv()*z)[0];assert ss>schur
x=P([0,1]);beta=P([Q(19,6),-Q(13,132)]);H=[beta**p for p in [6,10,14]];correction=P([0])
for i in range(3):
 correction+=Q(str(2*z[i]))/q*H[i]
 for j in range(3):correction+=sum((co[k]*Q(str(B[k][i,j]))for k in range(30)),Q(0))/(q*q)*H[i]*H[j]
ann=json.load(open(D/'scalar/annular_certificate.json'));assert ann['dimension']==51 and ann['d']=='22' and ann['cut']=='180';assert ann['kernel_sha256']==hashlib.sha256(raw).hexdigest();assert list(map(str,correction))==ann['kernel_polynomial'];assert correction.degree()==28 and correction(Q(418,13))==0;assert ann['s']==ker['s']
f=json.load(open(D/'scalar/repaired.json'));fp=P(list(map(Q,f['P'])));fq=P(list(map(Q,f['Q'])));assert f['dimension']==51 and f['contact_x']=='22' and fp.degree()==fq.degree()==71 and fq(0)==1;assert f['P']==ann['P']
factor=P([1])
for a in f['factors']:
 if a['type']=='negative_linear':
  v=Q(a['a']);assert v>0;factor*=P([1,1/v])
 elif a['type']=='positive_double':
  v=Q(a['mid']);assert v>0;factor*=P([1,-2/v,1/(v*v)])
 elif a['type']=='complex_pair':
  r,i=Q(a['re']),Q(a['im']);assert i>0;factor*=P([1,-2*r/(r*r+i*i),1/(r*r+i*i)])
 else:raise AssertionError(a)
assert factor==fq
alpha=Q(49,2)
def binom(a,k):
 v=Q(1)
 for i in range(k):v*= (a-i)/(i+1)
 return v
def transform(poly):
 out=P([0])
 for j in range(72):
  val=sum((poly[k]*math.factorial(k)*binom(k+alpha,k-j)for k in range(j,72)),Q(0))*Q((-1)**j,math.factorial(j));out+=val*x**j
 return out
assert transform(fq)==fp and transform(fp)==fq
# Independent exact reconstruction of the numerical-to-rational input polynomial.
cand=json.load(open(D/'source_candidate.json'));base=json.load(open(D/'source_scalar.json'));initial=P(list(map(Q,base['Q'])))
for k,d in enumerate(cand['radial_delta']):
 L=P([(-2)**j*binom(k+alpha,k-j)/math.factorial(j)for j in range(k+1)])/binom(k+alpha,k)
 initial+=Q(str(F(format(d,'.16g'))))*(-1)**k*L
initial/=initial(0);assert list(map(str,initial))==json.load(open(D/'scalar/input_roots.json'))['Q']
# An independent sufficient tail test, avoiding complex root classification.
shift=[sum((fp[k]*math.comb(k,j)*Q(180)**(k-j)for k in range(j,72)),Q(0))for j in range(72)];tail_negative=all(c<=0 for c in shift)and shift[-1]<0
caps=json.load(open(D/'scalar/caps/result.json'));assert len(caps['rows'])==len(ann['rows'])==633;edge=Q(22);budget=Q(0);checked_caps=0
for row,orig in zip(caps['rows'],ann['rows']):
 assert Q(row['l'])==edge and row['l']==orig['l'] and row['r']==orig['r'] and row['height']==orig['height'];edge=Q(row['r']);assert 0<edge-Q(row['l'])<=22;h=Q(row['height']);assert h>=0
 if h:
  c=json.load(open(D/'scalar/caps'/row['cap_certificate']));t=Q(c['s']);assert t==1-Q(11)/edge and -1<=t<1;cK=P(list(map(Q,c['K'])));poly=(x-t)*cK*cK
  if c['b']==1:poly*=1+x
  else:assert c['b']==0
  g=[P([1]),x]
  for k in range(2,poly.degree()+1):g.append(((2*k+47)*x*g[-1]-(k-1)*g[-2])/(k+48))
  cG=list(map(Q,c['gegenbauer']));assert min(cG)>=0 and cG[0]>0;assert sum((a*b for a,b in zip(cG,g)),P([0]))==poly
  ratio=poly(1)/cG[0];count=ratio.numerator//ratio.denominator;assert count==c['cap']==row['cap'];checked_caps+=1;budget+=h*count
assert edge==180 and budget==Q(caps['budget']);objective=fp[0]+Q(ker['s'])+budget;assert objective>0
value=arb(objective)*arb(Q(11,2))**(arb(51)/2)/(arb(53)/2).gamma();density=json.load(open(D/'scalar/caps/density_result.json'));assert Q(density['safe_upper'])==Q('68506041/500000000000000')==Q(str(F(density['safe_decimal'])));assert value<arb(Q('68506041/500000000000000'))
report={'status':'EXACT_ALGEBRA_AND_NORMALIZATION_PASS_GEOMETRY_PENDING','frozen_sha256':hashlib.sha256(raw).hexdigest(),'scalar_sha256':hashlib.sha256((D/'scalar/repaired.json').read_bytes()).hexdigest(),'direct_tail_shift_all_nonpositive':tail_negative,'positive_annuli_caps_checked':checked_caps,'radial_correction_degree':28,'conditional_density':str(value),'conditional_safe_decimal':density['safe_decimal'],'origin_slack':str(ss-schur),'source_reconstruction':'explicit normalized Laguerre coefficients matched exact input','fourier':'independent closed coefficient transform and positive-factor identity'}
(D/'exact_result.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if tail_negative:(D/'tail_shift_certificate.json').write_text(json.dumps({'cut':180,'coefficients':list(map(str,shift))},indent=2))
