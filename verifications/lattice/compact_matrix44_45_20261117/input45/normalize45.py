import sys;sys.set_int_max_str_digits(0)
from fractions import Fraction as F
from pathlib import Path
import json,math,decimal
from flint import arb,ctx
ctx.prec=1024
for n in [45]:
 p=Path(f'research/lattice4445_independent_20261117/input45/root_repair{n}');z=json.load(open(p/'annular_final_certificate.json'));obj=F(z['P'][0])+F(z['s'])+F(z['budget']);assert obj>0
 if n%2==0:
  val=obj*F(22,4)**(n//2)/math.factorial(n//2);flt=float(val);safe=F(math.ceil(flt*10**16),10**16);assert val<safe;proof=dict(density_exact=str(val))
 else:
  m=(n+1)//2;num=obj*obj*F(22,4)**n*16**m*math.factorial(m)**2/math.factorial(2*m)**2
  def atan_bounds(q,N=200):
   lower=sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0));return lower,lower+F(1,(2*N+1)*q**(2*N+1))
  l5,u5=atan_bounds(5);l239,u239=atan_bounds(239);piL=16*l5-4*u239;flt=math.sqrt(float(num)/math.pi);safe=F(math.ceil(flt*10**16),10**16);assert safe*safe*piL>num;proof=dict(density_squared_times_pi=str(num),pi_lower=str(piL))
 assert safe>0
 print(n,flt,'safe',float(safe),'budget',float(F(z['budget'])),flush=True);json.dump(dict(n=n,safe_upper=str(safe),exact_comparison_pass=True,status='PRODUCER PASS; independent audit pending',**proof),open(p/'density_result.json','w'),indent=2)
