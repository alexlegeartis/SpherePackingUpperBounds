from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib,sympy as s,sys
from flint import fmpq as Q,fmpq_poly as P
sys.set_int_max_str_digits(0);D=Path(__file__).parent;source=D.parent/'source48/input/frozen.json';src=json.loads(source.read_text());sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha=='169ea8d0e2ec662e7d450641894d1e3ba5f280623d2d6a659d1308401d45a728'
x=s.symbols('x');co=list(map(s.Rational,src['co']));Bs=[s.Matrix([[s.Rational(v)for v in row]for row in b])for b in src['B']];basis48=[s.Poly(s.gegenbauer(2*k,23,x)/s.gegenbauer(2*k,23,1),x)for k in range(len(co))]
for n,safe in ((45,F('1.1728551765e-6')),(46,F('7.721596228e-7')),(47,F('5.130033534e-7'))):
 I=D/f'dim{n}';R=I/'root_repair';z=json.loads((I/'frozen.json').read_text());c=json.loads((R/'input_candidate.json').read_text());an=json.loads((R/'annular_final_certificate.json').read_text());r=json.loads((R/'repaired.json').read_text());rr=json.loads((I/'restriction.json').read_text())
 assert z==c and z['source_sha']==rr['source_sha']==sha and z['n']==n and z['powers']==[6,9,12] and z['d']==22
 for key in ['co','q','B','z','s']:assert an[key]==z[key]
 assert z['co']==src['co'] and z['q']==src['q'];a=s.Rational(z['a']);assert a>0
 B=[s.Matrix([[s.Rational(v)for v in row]for row in b])for b in z['B']];alpha=s.Rational(n,2)-1;basis=[s.Poly(s.gegenbauer(2*k,alpha,x)/s.gegenbauer(2*k,alpha,1),x)for k in range(len(co))]
 for i in range(3):
  for j in range(3):
   assert sum((basis[k]*co[k]*B[k][i,j]for k in range(len(co))),s.Poly(0,x))==a*sum((basis48[k]*co[k]*Bs[k][i,j]for k in range(len(co))),s.Poly(0,x))
 for k,b in enumerate(B):
  assert b==b.T and all(b[:i,:i].det()>0 for i in(1,2,3));assert b==a*s.Matrix([[s.Rational(v)for v in row]for row in rr['C'][k]])/co[k]
 Z=s.Matrix(list(map(s.Rational,z['z'])));origin=F(z['s']);assert s.Rational(z['s'])-(Z.T*B[0].inv()*Z)[0]==s.Rational(1,10**12)
 p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert p.degree()==71 and ft(0)==1 and an['P']==r['P']
 xx=P([0,1]);al=Q(str(alpha));L=[P([1]),P([al+1,-1])]
 for k in range(2,72):L.append(((2*k+al-1-xx)*L[-1]-(k+al-1)*L[-2])/k)
 assert sum((p[k]*factorial(k)*L[k]for k in range(72)),P([0]))==ft
 fact=P([1])
 for v in r['factors']:
  if v['type']=='negative_linear':aa=Q(v['a']);assert aa>0;fact*=P([1,1/aa])
  elif v['type']=='positive_double':aa=Q(v['mid']);assert aa>0;fact*=P([1,-2/aa,1/aa**2])
  else:
   assert v['type']=='complex_pair';aa,bb=Q(v['re']),Q(v['im']);assert bb>0;fact*=P([1,-2*aa/(aa*aa+bb*bb),1/(aa*aa+bb*bb)])
 assert fact==ft
 last=F(22);budget=F(0)
 for row in an['rows']:
  l,u=F(row['l']),F(row['r']);assert l==last and 0<u-l<=88;last=u;assert not l<F(418,13)<u
  rad=u/22
  if rad<2:cap=1
  else:
   t=1-2/rad;mu=F(1);vals=[]
   for k in range(501):
    rat=F(2*k+1,n+2*k)
    if rat>t*t:vals.append((1-t*t)/(mu*(rat-t*t)))
    mu*=rat
   cap=min(vals).__floor__()
  assert cap==row['cap'];height=F(row['height']);assert height>=0;budget+=2**(n+1)*cap*height
 assert last==180 and budget==F(an['budget']) and len(an['rows'])==633
 obj=F(str(p[0]))+origin+budget;assert obj>0
 if n%2==0:val=obj*F(22,4)**(n//2)/factorial(n//2);assert val<safe
 else:
  m=(n+1)//2;num=obj*obj*F(22,4)**n*16**m*factorial(m)**2/factorial(2*m)**2
  def atan(q):
   lo=sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(200)),F(0));return lo,lo+F(1,401*q**401)
  lo5,hi5=atan(5);lo239,hi239=atan(239);pilower=16*lo5-4*hi239;assert safe**2*pilower>num
 out={'status':'EXACT TRANSFER AND SCALAR ALGEBRA PASS; SOURCE GEOMETRY STILL REQUIRED','n':n,'safe':str(safe),'source_sha':sha,'frozen_sha':hashlib.sha256((I/'frozen.json').read_bytes()).hexdigest(),'annuli':633};(D/f'exact_result{n}.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
