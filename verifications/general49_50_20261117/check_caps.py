from pathlib import Path
from fractions import Fraction as F
import json,math,sys
from flint import fmpq as Q,fmpq_poly as Poly
sys.set_int_max_str_digits(0);P=Path(__file__).resolve().parent;X=Poly([0,1]);results=[]
def atan(q,N):return sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
piL=16*atan(5,201)-4*atan(239,200)
for n,safe in [(49,'2.85041752e-7'),(50,'2.04993038e-7')]:
 D=P/f'caps{n}'
 S=P/f'dim{n}'
 original=json.loads((S/'annular_final_certificate.json').read_text());new=json.loads((D/'result.json').read_text());w=json.loads((S/'repaired.json').read_text())
 assert len(original['rows'])==len(new['rows'])==632
 last=Q(22);budget=Q(0);checked={}
 for v,old in zip(new['rows'],original['rows']):
  l,r=Q(v['l']),Q(v['r']);assert l==last and r-l==Q(1,4);last=r
  assert v['l']==old['l'] and v['r']==old['r'];h=Q(old['height']);assert h>=0
  # Use the independently replayed height, not an unchecked new radial value.
  cap=old['cap'];assert cap>0
  file=v['cap_certificate']
  if file is not None:
   c=json.loads((D/file).read_text());s=Q(c['s']);assert s==1-Q(11)/r and 0<s<1 and c['b'] in [0,1]
   K=Poly([Q(a)for a in c['K']]);f=(X-s)*K*K
   if c['b']==1:f*=X+1
   # Independent normalized Gegenbauer recurrence, not producer closed formula.
   lam=Q(n-2,2);gs=[Poly([1]),X]
   for k in range(1,f.degree()):gs.append((2*(k+lam)*X*gs[k]-k*gs[k-1])/(k+2*lam))
   rem=f;co=[Q(0)]*(f.degree()+1)
   for k in range(f.degree(),-1,-1):co[k]=rem[k]/gs[k][k];rem-=co[k]*gs[k]
   assert rem==0 and co[0]>0 and all(a>=0 for a in co)
   assert co==list(map(Q,c['gegenbauer']))
   ratio=f(1)/co[0];C=int(ratio.numerator//ratio.denominator);assert C==c['cap']
   cap=min(cap,C);checked[file]=f.degree()
  assert cap==v['cap'];budget+=h*cap
 assert last==180
 a,z,s=map(Q,[original['a'],original['z'],original['s']]);assert a>0 and s-z*z/a==Q(1,10**12)
 assert list(map(Q,original['P']))==list(map(Q,w['P']))
 objective=Q(w['P'][0])+s+budget;assert objective>0
 if n%2==0:
  density=F(str(objective))*F(22,4)**(n//2)/math.factorial(n//2);assert density<F(safe);desc={'exact_density':str(density)}
 else:
  m=(n+1)//2;square=F(str(objective))**2*4*22**n*math.factorial(m)**2/math.factorial(2*m)**2
  assert square<F(safe)**2*piL;desc={'exact_pi_times_density_squared':str(square)}
 results.append({'dimension':n,'safe_upper':safe,'status':'PASS independent cap recurrence, complete coverage, independent heights and exact safe rounding','distinct_polynomials':len(checked),'budget':str(budget),**desc})
 print(n,'PASS',safe,'caps',len(checked),'budget',float(budget),flush=True)
(P/'result.json').write_text(json.dumps(results,indent=2))
