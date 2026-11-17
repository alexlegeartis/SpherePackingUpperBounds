import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent;I=D/'input';R=I/'root_repair46'
z=json.loads((I/'frozen.json').read_text());c=json.loads((R/'input_candidate.json').read_text());a=json.loads((R/'annular_final_certificate.json').read_text());r=json.loads((R/'repaired.json').read_text())
assert hashlib.sha256((I/'frozen.json').read_bytes()).hexdigest()=='fe60310a72412cc57e47dc0f0825a0aa7806dd80d3c7c0c045888c6925957dab'
for key in ['co','B','q']:assert z[key]==c[key],key
for key in ['C','a','z','s','q']:assert a[key]==c[key],key
assert F(c['a'])==1 and F(a['L'])==1
assert F(c['z'][0])==F(z['z'][0])-F(1,200000) and F(c['z'][1])==F(z['z'][1])
co=list(map(F,z['co']));q=F(z['q']);assert co[0]==1 and min(co)>=0 and sum(co)==q
Bs=[[[F(v)for v in row]for row in b]for b in z['B']];zz=list(map(F,c['z']));origin=F(c['s'])
assert [[[F(v)*co[k]for v in row]for row in b]for k,b in enumerate(z['B'])]==[[[F(v)for v in row]for row in b]for b in c['C']]
for b in Bs:assert b[0][1]==b[1][0] and b[0][0]>0 and b[0][0]*b[1][1]>b[0][1]**2
b=Bs[0];det=b[0][0]*b[1][1]-b[0][1]**2;assert origin*det>=b[1][1]*zz[0]**2-2*b[0][1]*zz[0]*zz[1]+b[0][0]*zz[1]**2
p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert p.degree()==71 and ft(0)==1 and a['P']==r['P']
x=P([0,1]);lag=[P([1]),P([23,-1])]
for k in range(2,72):lag.append(((2*k+21-x)*lag[-1]-(k+21)*lag[-2])/k)
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
   rat=F(2*k+1,46+2*k)
   if rat>t*t:vals.append((1-t*t)/(mu*(rat-t*t)))
   mu*=rat
  v=min(vals);cap=v.numerator//v.denominator
 assert cap==row['cap'];h=F(row['height']);assert h>=0;budget+=2**47*cap*h
assert last==180 and budget==F(a['budget']) and len(a['rows'])==633
rho=(F(str(p[0]))+origin+budget)*F(22,4)**23/factorial(23);safe=F('7.766101459e-7');assert rho<safe
out=dict(status='EXACT PSD, FOURIER, CAP, COVERAGE, NORMALIZATION AND ROUNDING PASS',safe_bound=str(safe),density_approx=float(rho),budget_approx=float(budget),annuli=len(a['rows']),radial_nodes=a['nodes'])
(D/'exact_result.json').write_text(json.dumps(out,indent=2));print(out)
