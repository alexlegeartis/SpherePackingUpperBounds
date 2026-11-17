import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib,sympy as S
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent;I=D/'input';R=I/'root_repair';load=lambda p:json.loads(p.read_text())
z=load(I/'frozen.json');a=load(R/'annular_final_certificate.json');r=load(R/'repaired.json');raw=load(I/'numerical_source.json')
assert hashlib.sha256((I/'frozen.json').read_bytes()).hexdigest()=='4d67f38f02a67dfba251f088e471865d2919f4b95e864aa8786e293aa657ca01'
assert z==load(R/'input_candidate.json');assert z['n']==44 and z['degree']==71 and z['d']==22 and z['powers']==[6,9,12] and z['base_power']==6 and z['feature_stride']==3
for k in ['B','co','q','z','s']:assert z[k]==a[k]
co=list(map(F,z['co']));assert co[0]==1 and min(co)>=0 and sum(co)==F(z['q'])
B=[S.Matrix([[S.Rational(t)for t in row]for row in b])for b in z['B']];zz=S.Matrix(list(map(S.Rational,z['z'])));ss=S.Rational(z['s']);T=S.Matrix([[1,-2,1],[0,2,-2],[0,0,1]])
for k,b in enumerate(B):
 assert b==b.T and all(b[:j,:j].det()>0 for j in [1,2,3])
 src=S.Matrix([[S.Rational(str((raw['B'][k][i][j]+raw['B'][k][j][i])/2))for j in range(3)]for i in range(3)]);expected=T.T*src*T
 if k==0:expected[0,0]-=S.Rational(1,10000)
 assert b==expected
zsrc=T.T*S.Matrix([S.Rational(str(v))for v in raw['z']]);zsrc[0]-=S.Rational(1,200000);assert zz==zsrc
assert ss-(zz.T*B[0].inv()*zz)[0]==S.Rational(1,10**12)
p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert p.degree()==71 and ft.degree()==71 and ft(0)==1 and a['P']==r['P'];assert r['dimension']==44 and r['contact_x']=='22'
x=P([0,1]);L=[P([1]),P([Q(22),-1])]
for k in range(2,72):L.append(((Q(2*k)+Q(20)-x)*L[-1]-(Q(k)+Q(20))*L[-2])/k)
assert sum((p[k]*factorial(k)*L[k]for k in range(72)),P([0]))==ft
assert sum((ft[k]*factorial(k)*L[k]for k in range(72)),P([0]))==p
prod=P([1])
for f in r['factors']:
 if f['type']=='negative_linear':t=Q(f['a']);assert t>0;prod*=P([1,1/t])
 elif f['type']=='positive_double':t=Q(f['mid']);assert t>0;prod*=P([1,-2/t,1/t**2])
 else:
  assert f['type']=='complex_pair';u,v=Q(f['re']),Q(f['im']);assert v>0;prod*=P([1,-2*u/(u*u+v*v),1/(u*u+v*v)])
assert prod==ft
last=F(22);budget=F(0)
for row in a['rows']:
 l,u=F(row['l']),F(row['r']);assert l==last and 0<u-l<=88;last=u;assert not l<F(418,13)<u
 t=u/22
 if t<2:cap=1
 else:
  c=1-2/t;mu=F(1);vals=[]
  for k in range(501):
   rat=F(2*k+1,44+2*k)
   if rat>c*c:vals.append((1-c*c)/(mu*(rat-c*c)))
   mu*=rat
  best=min(vals);cap=best.numerator//best.denominator
 assert cap==row['cap'];height=F(row['height']);assert height>=0;budget+=2**45*cap*height
assert last==180 and len(a['rows'])==633 and budget==F(a['budget'])
obj=F(str(p[0]))+F(ss)+budget;assert obj>0
# In dimension44 the volume multiplier is exactly (22/4)^22 / 22!.
rho=obj*F(22,4)**22/factorial(22);safe=F('1.7676809433e-6');assert rho<safe
stored=load(R/'density_result.json');assert F(stored['safe_upper'])==safe and F(stored['density_exact'])==rho
out=dict(status='EXACT44D ALGEBRA, PSD, FOURIER, FACTORS, CAP COUNTS AND NORMALIZATION PASS; GEOMETRY AND RADIAL REPLAY SEPARATE',safe=str(safe),annuli=633,budget=str(budget));(D/'exact_result.json').write_text(json.dumps(out,indent=2));print(out)
