if not __debug__: raise RuntimeError('Do not disable assertions')
import json,sys,hashlib,math
from pathlib import Path
from fractions import Fraction as F
import sympy as S
from flint import fmpq as Q,fmpq_poly as P
sys.set_int_max_str_digits(0)
D=Path(__file__).parent;  z=json.load(open(D/'replay/input_candidate.json'));src=json.load(open(D/'source37/input/frozen.json'));r=json.load(open(D/'replay/repaired.json'));ann=json.load(open(D/'replay/annular_final_certificate.json'))
assert hashlib.sha256((D/'source37/input/frozen.json').read_bytes()).hexdigest()=='78f99a960215ea1d0c44e2e93faf0cb1e06c4be1dba99cb14495433c2b09c550'
assert z['n']==36 and z['d']==17 and z['powers']==[6,9,12]
assert z['co']==src['co'] and z['q']==src['q']
co=list(map(S.Rational,z['co'])); B=[S.Matrix(b) .applyfunc(S.Rational) for b in z['B']]; Bs=[S.Matrix(b).applyfunc(S.Rational) for b in src['B']]; a=S.Rational(z['a']);assert a>0
x=S.symbols('x');Gs=[S.gegenbauer(2*k,S.Rational(35,2),x)/S.gegenbauer(2*k,S.Rational(35,2),1) for k in range(len(co))];Gt=[S.gegenbauer(2*k,17,x)/S.gegenbauer(2*k,17,1) for k in range(len(co))]
for i in range(3):
 for j in range(3):assert S.expand(sum(co[k]*(B[k][i,j]*Gt[k]-a*Bs[k][i,j]*Gs[k]) for k in range(len(co))))==0
for b in B:assert b==b.T and all(b[:k,:k].det()>0 for k in [1,2,3])
Z=S.Matrix(z['z']).applyfunc(S.Rational); ss=S.Rational(z['s']);assert ss-(Z.T*B[0].inv()*Z)[0]==S.Rational(1,10**12)
p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert ft(0)==1 and len(p)==82
prod=P([1])
for f in r['factors']:
 if f['type']=='negative_linear':v=Q(f['a']);assert v>0;prod*=P([1,1/v])
 elif f['type']=='positive_double':v=Q(f['mid']);assert v>0;prod*=P([1,-2/v,1/v**2])
 else:
  assert f['type']=='complex_pair';u,v=Q(f['re']),Q(f['im']);assert v>0;prod*=P([1,-2*u/(u*u+v*v),1/(u*u+v*v)])
assert prod==ft
X=P([0,1]);L=[P([1]),P([18,-1])]
for k in range(2,82):L.append(((2*k+16-X)*L[-1]-(k+16)*L[-2])/k)
assert sum((ft[k]*math.factorial(k)*L[k] for k in range(82)),P([0]))==p
assert sum((p[k]*math.factorial(k)*L[k] for k in range(82)),P([0]))==ft
assert all(v<0 for v in p(P([180,1])))
for key in ['B','co','z','s','q']:assert ann[key]==z[key]
assert ann['P']==r['P']
last=F(17);budget=F(0)
for row in ann['rows']:
 l,u=F(row['l']),F(row['r']);assert l==last and 0<u-l<=68;last=u;assert not l<F(323,13)<u
 t=u/17
 if t<2:cap=1
 else:
  s=1-2/t;moment=F(1);values=[]
  for k in range(501):
   ratio=F(2*k+1,36+2*k)
   if ratio>s*s:values.append((1-s*s)/(moment*(ratio-s*s)))
   moment*=ratio
  best=min(values);cap=best.numerator//best.denominator
 assert cap==row['cap'];h=F(row['height']);assert h>=0;budget+=2**37*cap*h
assert last==180 and len(ann['rows'])==653 and budget==F(ann['budget'])
density=(F(str(p[0]))+F(z['s'])+budget)*F(17,4)**18/math.factorial(18);assert density<F('2.727705e-5')
print('PASS exact linkage to published37, scaled angular restriction, origin Schur, Fourier factors and involution, full tail, cap algebra, annular budget and even normalization',float(density))
(D/'exact_result.json').write_text(json.dumps({'status':'EXACT PASS; fresh geometry and radial replay separate','density':str(density),'safe':'2.727705e-5'},indent=2))
