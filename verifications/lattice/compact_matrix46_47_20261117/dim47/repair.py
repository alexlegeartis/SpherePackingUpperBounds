import sys
if not __debug__:raise RuntimeError("Assertions must remain enabled")
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
import json,math,time
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
ctx.prec=768;D=Path(__file__).parent;source=json.load(open(D/'input_roots.json'));old=Poly(list(map(Q,source['Q'])));roots=old.complex_roots();neg=[];pos=[];complexes=[]
for z,m in roots:
 assert m==1
 if z.imag.contains(0):
  assert z.real<0 or z.real>0
  (neg if z.real<0 else pos).append(z.real)
 elif z.imag>0:complexes.append(z)
assert len(neg)>=1 and len(pos)%2==0 and len(neg)+len(pos)+2*len(complexes)==old.degree();pos.sort(key=lambda x:float(x.mid()))
def rat(x):
 # Rational approximation is an exact chosen coefficient, not an interval claim.
 return Q(str(F(x.mid().str(110, radius=False))))
q=Poly([1]);factors=[]
for nr in neg:
 a=-rat(nr);assert a>0;q*=Poly([1,1/a]);factors.append(dict(type='negative_linear',a=str(a)))
for i in range(0,len(pos),2):
 mid=rat((pos[i]+pos[i+1])/2);assert mid>0;factor=Poly([1,-2/mid,1/(mid*mid)]);q*=factor;factors.append(dict(type='positive_double',mid=str(mid),old_width=str(pos[i+1]-pos[i])))
for z in complexes:
 re,im=rat(z.real),rat(z.imag);assert im>0;norm=re*re+im*im;factor=Poly([1,-2*re/norm,1/norm]);q*=factor;factors.append(dict(type='complex_pair',re=str(re),im=str(im)))
assert q.degree()==71 and q(0)==1
# Exact Fourier involution on polynomial times exp(-pi*r^2).
x=Poly([0,1]);Ls=[Poly([1]),Poly([Q("47/2"),-1])]
for k in range(2,72):Ls.append(((2*k+Q("43/2")-x)*Ls[-1]-(k+Q("43/2"))*Ls[-2])/k)
p=Poly([0])
for k in range(72):p+=q[k]*math.factorial(k)*Ls[k]
back=Poly([0])
for k in range(72):back+=p[k]*math.factorial(k)*Ls[k]
assert back==q and p[p.degree()]<0
out=dict(dimension=47,degree=71,x_convention='pi*r^2',contact_x='22',Q=list(map(str,q)),P=list(map(str,p)),factors=factors,fourier_nonnegative_by_positive_factors=True,exact_fourier_involution=True,spatial_certified=False,source='input_candidate.json')
(D/'repaired.json').write_text(json.dumps(out,indent=2));print('POSITIVE FACTORS AND EXACT TRANSFORM PASS','positive_real_pairs',len(pos)//2,'complex_pairs',len(complexes),'newf0',float(p[0]),flush=True)
proots=p.complex_roots();possible=[(z,m)for z,m in proots if z.imag.contains(0)and not z.real<0];print('PHYSICAL REAL ROOTS',[(str(z.real),m)for z,m in possible],flush=True);(D/'physical_roots.json').write_text(json.dumps([(str(z),int(m))for z,m in proots],indent=2))
