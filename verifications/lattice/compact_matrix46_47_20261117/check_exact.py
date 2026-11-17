import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,hashlib
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent;source=json.loads((D/'source48.json').read_text());assert hashlib.sha256((D/'source48.json').read_bytes()).hexdigest()=='aca9b48134ba9456995cd6a0cfc0aafcc419f8acdf6fe10f42af5477ca59cc86'
x=P([0,1])
def gegen(n,N):
 lam=Q(n-2,2);gs=[P([1]),2*lam*x]
 for k in range(2,N+1):gs.append((2*(k+lam-1)*x*gs[-1]-(k+2*lam-2)*gs[-2])/k)
 return [v/v(1)for v in gs]
g48=gegen(48,20);co=list(map(Q,source['co']));bs=[[[Q(v)for v in r]for r in b]for b in source['B']];original=[[sum((co[k]*bs[k][i][j]*g48[2*k]for k in range(len(co))),P([0]))for j in range(2)]for i in range(2)]
results=[]
for n,safe in [(46,'7.790880151e-7'),(47,'5.186745468e-7')]:
 S=D/f'dim{n}';z=json.loads((S/'input_candidate.json').read_text());r=json.loads((S/'repaired.json').read_text());a=json.loads((S/'annular_final_certificate.json').read_text());rest=json.loads((D/f'restriction{n}.json').read_text());C=[[[Q(v)for v in row]for row in b]for b in z['C']];g=gegen(n,2*(len(C)-1));assert z['C']==rest['C'] and z['q']==source['q'];q=F(z['q']);amp=F(z['a']);assert amp>0
 for i in range(2):
  for j in range(2):assert sum((C[k][i][j]*g[2*k]for k in range(len(C))),P([0]))==original[i][j]
 for b in C:assert b[0][1]==b[1][0] and b[0][0]>=0 and b[1][1]>=0 and b[0][0]*b[1][1]>=b[0][1]**2
 zz=list(map(F,z['z']));origin=F(z['s']);b=[[amp*F(str(v))for v in row]for row in C[0]];det=b[0][0]*b[1][1]-b[0][1]**2;assert det>0 and origin*det>=b[1][1]*zz[0]**2-2*b[0][1]*zz[0]*zz[1]+b[0][0]*zz[1]**2
 for k in ['C','a','q','z','s']:assert a[k]==z[k],k
 p=P(list(map(Q,r['P'])));ft=P(list(map(Q,r['Q'])));assert p.degree()==71 and ft(0)==1 and a['P']==r['P'];alpha=Q(n,2)-1;lag=[P([1]),P([alpha+1,-1])]
 for k in range(2,72):lag.append(((2*k-1+alpha-x)*lag[-1]-(k-1+alpha)*lag[-2])/k)
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
  l,u=F(row['l']),F(row['r']);assert l==last and 0<u-l<=88;last=u;assert not l<F(418,13)<u;rad=u/22
  if rad<2:cap=1
  else:
   t=1-2/rad;mu=F(1);vals=[]
   for k in range(501):
    rat=F(2*k+1,n+2*k)
    if rat>t*t:vals.append((1-t*t)/(mu*(rat-t*t)))
    mu*=rat
   v=min(vals);cap=v.numerator//v.denominator
  assert cap==row['cap'];h=F(row['height']);assert h>=0;budget+=2**(n+1)*cap*h
 assert last==180 and budget==F(a['budget']) and len(a['rows'])==633
 obj=F(str(p[0]))+origin+budget;assert obj>0
 if n==46:assert obj*F(22,4)**23/factorial(23)<F(safe)
 else:
  def atan(q,N):return sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
  piL=16*atan(5,201)-4*atan(239,200);num=obj*obj*F(22,4)**47*16**24*factorial(24)**2/factorial(48)**2;assert num<F(safe)**2*piL
 results.append(dict(n=n,status='PASS separate exact restriction, PSD, Fourier, caps, coverage and volume',safe_upper=safe,radial_nodes=a['nodes'],budget_approx=float(budget)))
(D/'exact_result.json').write_text(json.dumps(results,indent=2));print(results)
