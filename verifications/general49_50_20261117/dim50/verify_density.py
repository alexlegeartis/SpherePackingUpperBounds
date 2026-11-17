import json,math,sys
from pathlib import Path
from fractions import Fraction as F
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);D=Path(__file__).parent
r=json.loads((D/'repaired.json').read_text());c=json.loads((D/'annular_final_certificate.json').read_text());assert len(c['rows'])==632
budget=F(0)
for j,row in enumerate(c['rows']):
 assert F(row['l'])==22+F(j,4) and F(row['r'])==22+F(j+1,4)
 ss=1-F(11)/F(row['r']);k=row['cap_degree'];gap=F(k,k+49)-ss;assert isinstance(k,int) and k>=0 and gap>0
 rat=F(1)
 for h in range(k):rat*=F(49+h)/F(49+2*h,2)
 cap=(1-ss)*rat/gap;assert row['cap']==cap.numerator//cap.denominator
 assert F(row['height'])>=0;budget+=row['cap']*F(row['height'])
assert budget==F(c['budget']);a,z,s=map(F,[c['a'],c['z'],c['s']]);assert a>0 and s==z*z/a+F(1,10**12) and s*a-z*z>0 and a+2*z<0
assert F(c['q'])==F(json.loads((D/'angular.json').read_text())['B_exact'])
assert c['P']==r['P'] and r['dimension']==50 and F(r['contact_x'])==22
obj=F(r['P'][0])+s+budget;density=F(22,4)**25/math.factorial(25)*obj;safe=F('2.06573053e-7');assert density<safe
out=dict(status='PRODUCER COMPLETE UNRESTRICTED50 UPPER PASS; INDEPENDENT AUDIT PENDING',density_exact=str(density),safe_upper=str(safe),safe_decimal='2.06573053e-7',budget=str(budget),objective=str(obj),annuli=632,nodes=c['nodes'])
(D/'density_result.json').write_text(json.dumps(out,indent=2));print('EXACT DENSITY PASS',float(density),'safe',float(safe),flush=True)
