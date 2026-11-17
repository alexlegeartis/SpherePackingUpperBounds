import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib,sympy as S
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent;I=D/'input';R=I/'root_repair';load=lambda p:json.loads(p.read_text())
z=load(I/'frozen.json');a=load(R/'annular_final_certificate.json');r=load(R/'repaired.json');raw=load(I/'numerical_source.json')
assert hashlib.sha256((I/'frozen.json').read_bytes()).hexdigest()=='fb7b8b04c23006dd56282013aa38a867f0d2552aace230cbffe8f55bc4336f5d'
assert z==load(R/'input_candidate.json');assert z['n']==40 and z['degree']==81 and z['d']==17 and z['powers']==[6,9,12] and z['base_power']==6 and z['feature_stride']==3
for k in ['B','co','q','z','s']:assert z[k]==a[k]
co=list(map(F,z['co']));assert co[0]==1 and min(co)>=0 and sum(co)==F(z['q'])
B=[S.Matrix([[S.Rational(t)for t in row]for row in b])for b in z['B']];zz=S.Matrix(list(map(S.Rational,z['z'])));ss=S.Rational(z['s']);T=S.Matrix([[1,-2,1],[0,2,-2],[0,0,1]])
for k,b in enumerate(B):
 assert b==b.T and all(b[:j,:j].det()>0 for j in [1,2,3])
 src=S.Matrix([[S.Rational(str((raw['B'][k][i][j]+raw['B'][k][j][i])/2))for j in range(3)]for i in range(3)]);expected=T.T*src*T
 if k==0:expected[0,0]-=S.Rational(1,100)
 assert b==expected
zsrc=T.T*S.Matrix([S.Rational(str(v))for v in raw['z']]);zsrc[0]-=S.Rational(1,1000);assert zz==zsrc
assert ss-(zz.T*B[0].inv()*zz)[0]==S.Rational(1,10**12)
p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert p.degree()==81 and ft.degree()==81 and ft(0)==1 and a['P']==r['P'];assert r['dimension']==40 and r['contact_x']=='17'
x=P([0,1]);L=[P([1]),P([Q(40,2),-1])]
for k in range(2,82):L.append(((Q(2*k)+Q(36,2)-x)*L[-1]-(Q(k)+Q(36,2))*L[-2])/k)
assert sum((p[k]*factorial(k)*L[k]for k in range(82)),P([0]))==ft
assert sum((ft[k]*factorial(k)*L[k]for k in range(82)),P([0]))==p
prod=P([1])
for f in r['factors']:
 if f['type']=='negative_linear':t=Q(f['a']);assert t>0;prod*=P([1,1/t])
 elif f['type']=='positive_double':t=Q(f['mid']);assert t>0;prod*=P([1,-2/t,1/t**2])
 else:
  assert f['type']=='complex_pair';u,v=Q(f['re']),Q(f['im']);assert v>0;prod*=P([1,-2*u/(u*u+v*v),1/(u*u+v*v)])
assert prod==ft
last=F(17);budget=F(0)
for row in a['rows']:
 l,u=F(row['l']),F(row['r']);assert l==last and 0<u-l<=68;last=u;assert not l<F(323,13)<u
 t=u/17
 if t<2:cap=1
 else:
  c=1-2/t;mu=F(1);vals=[]
  for k in range(501):
   rat=F(2*k+1,40+2*k)
   if rat>c*c:vals.append((1-c*c)/(mu*(rat-c*c)))
   mu*=rat
  best=min(vals);cap=best.numerator//best.denominator
 assert cap==row['cap'];height=F(row['height']);assert height>=0;budget+=2**41*cap*height
assert last==180 and len(a['rows'])==653 and budget==F(a['budget'])
obj=F(str(p[0]))+F(ss)+budget;assert obj>0
# Independent dimension40 volume normalization.
n=40;safe=F('5.9777359372e-6');stored=load(R/'density_result.json')
if n%2==0:
 rho=obj*F(17,4)**(n//2)/factorial(n//2);assert rho<safe and F(stored['density_exact'])==rho
else:
 m=(n+1)//2;num=obj**2*F(17,4)**n*16**m*factorial(m)**2/factorial(2*m)**2
 def atan(q):
  lo=sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(200)),F(0));return lo,lo+F(1,401*q**401)
 l5,u5=atan(5);l239,u239=atan(239);piL=16*l5-4*u239;assert piL>0 and safe**2*piL>num
 assert F(stored['pi_lower'])==piL and F(stored['density_squared_times_pi'])==num
assert F(stored['safe_upper'])==safe
out=dict(status='EXACT40D ALGEBRA, PSD, FOURIER, FACTORS, CAP COUNTS AND NORMALIZATION PASS; GEOMETRY AND RADIAL REPLAY SEPARATE',safe=str(safe),annuli=653,budget=str(budget));(D/'exact_result.json').write_text(json.dumps(out,indent=2));print(out)
