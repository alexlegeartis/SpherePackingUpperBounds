import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib,sympy as sy
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent;I=D/'input';R=I/'root_repair'
z=json.loads((I/'frozen.json').read_text());c=json.loads((R/'input_candidate.json').read_text());a=json.loads((R/'annular_final_certificate.json').read_text());r=json.loads((R/'repaired.json').read_text())
assert hashlib.sha256((I/'frozen.json').read_bytes()).hexdigest()=='169ea8d0e2ec662e7d450641894d1e3ba5f280623d2d6a659d1308401d45a728'
assert z==c and z['n']==48 and z['powers']==[6,9,12] and z['base_power']==6 and z['feature_stride']==3 and z['d']==22
for key in ['co','B','q','z','s']:assert a[key]==z[key],key
co=list(map(F,z['co']));q=F(z['q']);assert co[0]==1 and min(co)>=0 and sum(co)==q
Bs=[sy.Matrix([[sy.Rational(v)for v in row]for row in b])for b in z['B']];zz=sy.Matrix(list(map(sy.Rational,z['z'])));origin=sy.Rational(z['s'])
for b in Bs:assert b==b.T and all(b[:i,:i].det()>0 for i in (1,2,3))
assert origin-(zz.T*Bs[0].inv()*zz)[0]==sy.Rational(1,10**12)
# Exact chosen repairs, relative to preserved discovery decimals (not a numerical premise).
raw=json.loads((I/'numerical_source.json').read_text())
for k,b in enumerate(Bs):
 for i in range(3):
  for j in range(3):assert b[i,j]==sy.Rational(str((raw['B'][k][i][j]+raw['B'][k][j][i])/2))-(sy.Rational(1,10000)if(k,i,j)==(0,0,0)else 0)
for i in range(3):assert zz[i]==sy.Rational(str(raw['z'][i]))-(sy.Rational(1,200000)if i==0 else 0)
p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert p.degree()==71 and ft(0)==1 and a['P']==r['P']
x=P([0,1]);lag=[P([1]),P([24,-1])]
for k in range(2,72):lag.append(((2*k+22-x)*lag[-1]-(k+22)*lag[-2])/k)
assert sum((p[k]*factorial(k)*lag[k]for k in range(72)),P([0]))==ft
fact=P([1])
for v in r['factors']:
 if v['type']=='negative_linear':aa=Q(v['a']);assert aa>0;fact*=P([1,1/aa])
 elif v['type']=='positive_double':aa=Q(v['mid']);assert aa>0;fact*=P([1,-2/aa,1/aa**2])
 else:
  assert v['type']=='complex_pair';aa,bb=Q(v['re']),Q(v['im']);assert bb>0;fact*=P([1,-2*aa/(aa*aa+bb*bb),1/(aa*aa+bb*bb)])
assert fact==ft
last=F(22);budget=F(0)
for row in a['rows']:
 l,u=F(row['l']),F(row['r']);assert l==last and 0<u-l<=88;last=u
 assert not l<F(418,13)<u
 rad=u/22
 if rad<2:cap=1
 else:
  t=1-2/rad;mu=F(1);vals=[]
  for k in range(501):
   rat=F(2*k+1,48+2*k)
   if rat>t*t:vals.append((1-t*t)/(mu*(rat-t*t)))
   mu*=rat
  v=min(vals);cap=v.numerator//v.denominator
 assert cap==row['cap'];h=F(row['height']);assert h>=0;budget+=2**49*cap*h
assert last==180 and budget==F(a['budget']) and len(a['rows'])==633
rho=(F(str(p[0]))+F(origin)+budget)*F(22,4)**24/factorial(24);safe=F('3.4559190288e-7');assert rho<safe
out=dict(status='EXACT PSD, FOURIER, CAP, COVERAGE, NORMALIZATION AND ROUNDING PASS; MAIN GEOMETRY STILL REQUIRED',safe_bound=str(safe),density_approx=float(rho),budget_approx=float(budget),annuli=len(a['rows']),radial_nodes=a['nodes'])
(D/'exact_result.json').write_text(json.dumps(out,indent=2));print(out)
