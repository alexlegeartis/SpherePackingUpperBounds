import sys,json,math
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
ctx.prec=384;sys.set_int_max_str_digits(0)
D=Path(__file__).parent;n=int(sys.argv[1]);j=json.loads((D/f'dim{n}_radial_comparison.json').read_text());raw=j['outcomes'][-1]
def rat(x):return Q(str(F(float(x))))
X=Poly([0,1]);a=Q(n,2);gs=[Poly([1]),a-2*X]
for k in range(1,len(raw['radial'])-1):gs.append(((2*k+a-2*X)*gs[-1]-(k+a-1)*gs[-2])/(k+1))
gs=[p/p(0) for p in gs]
P=sum((rat(c)*g for c,g in zip(raw['radial'],gs)),Poly());H=sum(((-1)**k*rat(c)*g for k,(c,g) in enumerate(zip(raw['radial'],gs))),Poly());C=rat(raw['C']);N=H(0);P/=N;H/=N;C/=N
for epsilon in [Q(1,10**6),Q(1,10**5),Q(1,10**4)]:
 p=(P+epsilon*(1-4*X/n))/(1-epsilon);h=(H-epsilon*(1-4*X/n))/(1-epsilon);c=C/(1-epsilon)*Q(1000001,1000000)
 roots={}
 for name,pol in [('physical',p),('fourier',h)]:
  rr=pol.numer().complex_roots();assert sum(m for z,m in rr)==pol.degree();reals=[]
  for z,m in rr:
   if z.imag.contains(0):
    assert z.imag.is_exact() and z.imag==0
    if not z.real<0:reals.append(str(z.real))
  roots[name]=reals
 obj=p(0)+c;vol=arb.pi()**arb(a)*(arb(Q(n,32)))**arb(a)/arb(a+1).gamma();density=arb(obj)*vol
 print('EPS',str(epsilon),'ROOTS',roots,'DENSITY',str(density),flush=True)
 if not roots['fourier'] and h(0)>0:
  out=dict(dimension=n,physical=list(map(str,p.coeffs())),fourier=list(map(str,h.coeffs())),C_exact=str(c),epsilon_exact=str(epsilon),s_exact=str(Q(n,8)),density_expression=str(density),certified=False,roots=roots)
  (D/f'radial_exact{n}.json').write_text(json.dumps(out,indent=2));break
