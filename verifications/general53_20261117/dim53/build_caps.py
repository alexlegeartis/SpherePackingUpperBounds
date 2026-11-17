import sys,json,math,time
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from flint import fmpq as Q,fmpq_poly as P
from scipy.special import roots_jacobi
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);n=int(sys.argv[1]);base=Path(__file__).parent;out=base/'caps';out.mkdir(exist_ok=True);data=json.loads((base/'annular_certificate.json').read_text());x=P([0,1]);G=[P([1]),x]
def gegen(poly):
 while len(G)<=poly.degree():
  k=len(G);G.append(((2*k+n-4)*x*G[-1]-(k-1)*G[-2])/(k+n-3))
 rem=poly;co=[Q(0)]*(poly.degree()+1)
 for k in range(poly.degree(),-1,-1):co[k]=rem[k]/G[k][k];rem-=co[k]*G[k]
 assert rem==0
 return co
@lru_cache(None)
def jacobi(a,b,K):
 alpha,beta=Q(n-3,2)+a,Q(n-3,2)+b;J=[P([1])];h=[Q(1)if b==0 else Q(n-1,n)]
 for k in range(K):
  t=2*k+alpha+beta;A=(beta*beta-alpha*alpha)/(t*(t+2))
  if k==0:J.append(x-A)
  else:
   B=4*k*(k+alpha)*(k+beta)*(k+alpha+beta)/(t*t*(t*t-1));J.append((x-A)*J[-1]-B*J[-2]);h.append(B*h[-1])
 # h[k] for monic degree k uses recurrence B_k.
 if K:
  k=K;t=2*k+alpha+beta;B=4*k*(k+alpha)*(k+beta)*(k+alpha+beta)/(t*t*(t*t-1));h.append(B*h[-1])
 assert len(J)==len(h)==K+1
 return J,h

def candidate(s,k,b):
 J,h=jacobi(1,b,k-1);K=sum((Ji*Ji(s)/hi for Ji,hi in zip(J,h)),P([0]));f=(x-s)*K*K
 if b:f*=x+1
 co=gegen(f)
 if co[0]<=0 or min(co)<0:return None
 bound=f(1)/co[0];C=int(bound.numerator//bound.denominator)
 assert C>0 and sum(co)==f(1)
 return dict(s=str(s),k=k,b=b,K=list(map(str,K)),gegenbauer=list(map(str,co)),ratio=str(bound),cap=C)

newrows=[];budget=Q(0)
for j,row in enumerate(data['rows']):
 h=Q(row['height']);old=row['cap'];rec=None
 if h>0:
  s=1-Q(data['d'])/(2*Q(row['r']));sf=float(s);chosen=None
  for k in range(1,201):
   r10=roots_jacobi(k,(n-1)/2,(n-3)/2)[0][-1];r11=roots_jacobi(k,(n-1)/2,(n-1)/2)[0][-1]
   if sf<r10:chosen=(k,0);break
   if sf<r11:chosen=(k,1);break
  assert chosen
  for k,b in [chosen,(chosen[0]+1,chosen[1])]:
   rec=candidate(s,k,b)
   if rec is not None:break
  assert rec is not None,(j,chosen)
  C=min(old,rec['cap']);(out/f'cap_{j}.json').write_text(json.dumps(rec,indent=2))
  print(n,j,'k',rec['k'],'b',rec['b'],'improvement',old/C,'oldcost',float(h*old),'newcost',float(h*C),flush=True)
 else:C=old
 budget+=h*C;newrows.append({**row,'old_cap':old,'cap':C,'cap_certificate':None if rec is None else f'cap_{j}.json'})
result=dict(status='EXACT LEVENSHTEIN POLYNOMIAL CAP CERTIFICATES PASS',dimension=n,budget=str(budget),old_budget=data['budget'],rows=newrows)
(out/'result.json').write_text(json.dumps(result,indent=2));print('PASS',n,'oldbudget',float(Q(data['budget'])),'newbudget',float(budget),flush=True)
