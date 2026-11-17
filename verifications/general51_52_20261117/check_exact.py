from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,sys
from flint import fmpq as Q,fmpq_poly as Poly
sys.set_int_max_str_digits(0);P=Path(__file__).resolve().parent;results=[]
def atan(q,N):return sum((F((-1)**j,(2*j+1)*q**(2*j+1))for j in range(N+1)),F(0))
piL=16*atan(5,201)-4*atan(239,200)
for n,safe in [(51,'1.52983244e-7'),(52,'1.06516626e-7')]:
 S=P/f'dim{n}';r=json.loads((S/'repaired.json').read_text());ann=json.loads((S/'annular_certificate.json').read_text());contact=json.loads((S/'contact_cap_certificate.json').read_text());p=Poly(list(map(Q,r['P'])));q=Poly(list(map(Q,r['Q'])));d=Q(r['contact_x']);x=Poly([0,1]);alpha=Q(n,2)-1;Ls=[Poly([1]),Poly([alpha+1,-1])]
 for k in range(2,p.degree()+1):Ls.append(((2*k-1+alpha-x)*Ls[-1]-(k-1+alpha)*Ls[-2])/k)
 assert sum((p[k]*factorial(k)*Ls[k]for k in range(p.degree()+1)),Poly([0]))==q and q(0)==1
 fac=Poly([1])
 for v in r['factors']:
  if v['type']=='negative_linear':a=Q(v['a']);assert a>0;fac*=Poly([1,1/a])
  elif v['type']=='positive_double':a=Q(v['mid']);assert a>0;fac*=Poly([1,-2/a,1/a**2])
  else:
   assert v['type']=='complex_pair';a,b=Q(v['re']),Q(v['im']);assert b>0;fac*=Poly([1,-2*a/(a*a+b*b),1/(a*a+b*b)])
 assert fac==q
 a,z,s=map(Q,[ann['a'],ann['z'],ann['s']]);assert list(map(Q,ann['P']))==list(p)
 if n==51:
  assert a>0 and a+2*z<0 and s-z*z/a==Q(1,10**12)
  import hashlib;assert hashlib.sha256((S/'angular.json').read_bytes()).hexdigest()=='ec0b022471f3eb27a44dc75b1fc6f8ceed486098fa04911e346495b42b681243'
 else:assert a==z==s==0
 Cdir=P/f'caps{n}';caps=json.loads((Cdir/'result.json').read_text());assert len(caps['rows'])==len(ann['rows']);last=d;budget=Q(0);num=0
 for i,(row,crow)in enumerate(zip(ann['rows'],caps['rows'])):
  l,rr=Q(row['l']),Q(row['r']);assert l==last and 0<rr-l<=d;last=rr;assert row['l']==crow['l']and row['r']==crow['r'];height=Q(row['height']);assert height>=0
  if i==0:assert rr==d*(1+Q(contact['epsilon']));cap=contact['cap']
  else:
   k=row['cap_degree'];ss=1-d/(2*rr);gap=Q(k,k+n-1)-ss;assert gap>0
   # Exact factorial/rising-product reconstruction independent of floating selection.
   ratio=Q(1)
   for j in range(k):ratio*=Q(2*(n-1+j),n-1+2*j)
   v=(1-ss)*ratio/gap;cap=int(v.numerator//v.denominator)
  assert cap==row['cap']
  if crow['cap_certificate']is not None:
   cc=json.loads((Cdir/crow['cap_certificate']).read_text());ss=Q(cc['s']);assert ss==1-d/(2*rr)and -1<=ss<1
   K=Poly(list(map(Q,cc['K'])));f=(x-ss)*K*K;assert cc['b']in[0,1]
   if cc['b']:f*=x+1
   lam=Q(n-2,2);gs=[Poly([1]),x]
   for k in range(1,f.degree()):gs.append((2*(k+lam)*x*gs[k]-k*gs[k-1])/(k+2*lam))
   rem=f;co=[Q(0)]*(f.degree()+1)
   for k in range(f.degree(),-1,-1):co[k]=rem[k]/gs[k][k];rem-=co[k]*gs[k]
   assert rem==0 and min(co)>=0 and co[0]>0 and co==list(map(Q,cc['gegenbauer']))
   v=f(1)/co[0];cap=min(cap,int(v.numerator//v.denominator));num+=1
  assert cap==crow['cap'];budget+=cap*height
 assert last==(180 if n==51 else 250)
 objective=p[0]+s+budget;assert objective>0
 if n%2:
  m=(n+1)//2;v=F(str(objective))**2*4*F(str(d))**n*factorial(m)**2/factorial(2*m)**2;assert v<F(safe)**2*piL
 else:v=F(str(objective))*(F(str(d))/4)**(n//2)/factorial(n//2);assert v<F(safe)
 results.append({'dimension':n,'safe_upper':safe,'status':'PASS independent exact Fourier/caps/coverage/PSD/rounding','budget':str(budget),'annuli':len(ann['rows']),'nodes':ann['nodes']});print(n,'PASS',safe,float(budget),flush=True)
(P/'exact_result.json').write_text(json.dumps(results,indent=2))
