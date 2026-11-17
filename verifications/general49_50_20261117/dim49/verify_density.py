import json,math,sys
from pathlib import Path
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);D=Path(__file__).parent
r=json.loads((D/'repaired.json').read_text());c=json.loads((D/'annular_final_certificate.json').read_text());n=r['dimension'];d=F(r['contact_x']);assert d==22 and len(c['rows'])==632
budget=F(0)
for j,row in enumerate(c['rows']):
 assert F(row['l'])==22+F(j,4) and F(row['r'])==22+F(j+1,4)
 ss=1-F(11)/F(row['r']);k=row['cap_degree'];gap=F(k,k+n-1)-ss;assert isinstance(k,int) and k>=0 and gap>0
 rat=F(1)
 for h in range(k):rat*=F(n-1+h)/F(n-1+2*h,2)
 cap=(1-ss)*rat/gap;assert row['cap']==cap.numerator//cap.denominator
 assert F(row['height'])>=0;budget+=row['cap']*F(row['height'])
assert budget==F(c['budget']);a,z,s=map(F,[c['a'],c['z'],c['s']]);assert a>0 and s==z*z/a+F(1,10**12) and s*a-z*z>0 and a+2*z<0
assert F(c['q'])==F(json.loads((D/'angular.json').read_text())['B_exact']) and c['P']==r['P']
obj=F(r['P'][0])+s+budget;assert obj>0
out=dict(status='PRODUCER COMPLETE UNRESTRICTED UPPER; INDEPENDENT AUDIT PENDING',dimension=n,budget=str(budget),objective=str(obj),annuli=632,nodes=c['nodes'])
if n%2==0:
 density=(d/4)**(n//2)/math.factorial(n//2)*obj;out['density_exact']=str(density)
else:
 def atan(inv,N):return sum((F((-1)**k,(2*k+1)*inv**(2*k+1)) for k in range(N)),F(0))
 lo=16*atan(5,200)-4*atan(239,201);hi=16*atan(5,201)-4*atan(239,200);assert F(314159,100000)<lo<hi<F(314160,100000)
 m=(n+1)//2;sq=obj**2*4*math.factorial(m)**2*d**n/math.factorial(2*m)**2/lo
 scale=10**100;v=sq*scale**2;density=F(math.isqrt(v.numerator//v.denominator)+1,scale);assert density**2>sq
 out.update(pi_lower=str(lo),pi_upper=str(hi),squared_density_upper=str(sq),rational_density_upper=str(density))
safe=F(math.ceil(density*10**15),10**15);assert safe>=density
out.update(safe_upper=str(safe),safe_decimal=f'{float(safe):.15g}');(D/'density_result.json').write_text(json.dumps(out,indent=2));print(n,'density',float(density),'safe',float(safe),'budget',float(budget),flush=True)
