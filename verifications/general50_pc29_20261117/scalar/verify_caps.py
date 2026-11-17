import json,sys,math,hashlib
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);n=int(sys.argv[1]);S=Path(__file__).parent;D=S/'caps';original=json.loads((S/'annular_certificate.json').read_text());r=json.loads((D/'result.json').read_text());assert len(r['rows'])==len(original['rows']);x=P([0,1]);budget=Q(0)
def rise(a,k):return math.prod((a+j for j in range(k)),start=Q(1))
def expand(co):
 poly=P([0]);lam=Q(n-2,2)
 for k,c in enumerate(co):
  g=[Q(0)]*(k+1);norm=rise(2*lam,k)/math.factorial(k)
  for j in range(k//2+1):g[k-2*j]=(-1)**j*rise(lam,k-j)*2**(k-2*j)/(math.factorial(j)*math.factorial(k-2*j)*norm)
  poly+=c*P(g)
 return poly
for j,(v,old) in enumerate(zip(r['rows'],original['rows'])):
 assert v['l']==old['l'] and v['r']==old['r'] and v['height']==old['height'] and v['old_cap']==old['cap']
 if v['cap_certificate'] is not None:
  c=json.loads((D/v['cap_certificate']).read_text());s=Q(c['s']);assert s==1-Q(original['d'])/(2*Q(v['r'])) and -1<=s<1 and c['b']in[0,1]
  K=P(list(map(Q,c['K'])));f=(x-s)*K*K
  if c['b']:f*=x+1
  co=list(map(Q,c['gegenbauer']));assert min(co)>=0 and co[0]>0 and expand(co)==f
  ratio=f(1)/co[0];assert ratio==Q(c['ratio']);C=ratio.numerator//ratio.denominator;assert C==c['cap'] and v['cap']==min(old['cap'],C)
 else:assert v['cap']==old['cap']
 budget+=Q(v['height'])*v['cap']
assert budget==Q(r['budget']);s=F(original['s']);w=json.loads((S/'repaired.json').read_text());objective=F(w['P'][0])+s+F(str(budget));assert objective>0
if n%2==0:density=(F(original['d'])/4)**(n//2)/math.factorial(n//2)*objective;bound={'density_exact':str(density)}
else:
 def atan(v,N):return sum((F((-1)**k,(2*k+1)*v**(2*k+1))for k in range(N)),F(0))
 pi_lower=16*atan(5,200)-4*atan(239,201);m=(n+1)//2;sq=objective**2*4*math.factorial(m)**2*F(original['d'])**n/math.factorial(2*m)**2/pi_lower;scale=10**100;t=sq*scale**2;density=F(math.isqrt(t.numerator//t.denominator)+1,scale);assert density*density>sq;bound={'squared_density_upper':str(sq),'pi_lower':str(pi_lower),'density_rational_upper':str(density)}
safe=F(math.ceil(density*10**15),10**15);assert density<=safe
out=dict(status='EXACT RECONSTRUCTED CAPS AND DENSITY PASS; INDEPENDENT AUDIT PENDING',dimension=n,safe_upper=str(safe),safe_decimal=f'{float(safe):.15g}',new_budget=str(budget),source_annular_sha256=hashlib.sha256((S/'annular_certificate.json').read_bytes()).hexdigest(),**bound);(D/'density_result.json').write_text(json.dumps(out,indent=2));print(n,'safe',out['safe_decimal'],'density',float(density),'budget',float(budget),flush=True)
