import json,math,sys
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as Poly
sys.set_int_max_str_digits(0)
if not __debug__:raise RuntimeError('Assertions required')
D=Path(__file__).parent;z=json.load(open(D/'root_repair45/repaired.json'));assert z['dimension']==45 and z['contact_x']=='22' and z['degree']==71
q=Poly([1])
for f in z['factors']:
 t=f['type']
 if t=='negative_linear':
  a=Q(f['a']);assert a>0;q*=Poly([1,1/a])
 elif t=='positive_double':
  m=Q(f['mid']);assert m>0;q*=Poly([1,-2/m,1/(m*m)])
 elif t=='complex_pair':
  r=Q(f['re']);i=Q(f['im']);assert i>0;nn=r*r+i*i;q*=Poly([1,-2*r/nn,1/nn])
 else:raise AssertionError(t)
assert q==Poly(list(map(Q,z['Q']))) and q[0]==1 and q.degree()==71
x=Poly([0,1]);L=[Poly([1]),Poly([Q(45,2),-1])]
for k in range(2,72):L.append(((2*k+Q(41,2)-x)*L[-1]-(k+Q(41,2))*L[-2])/k)
p=sum((q[k]*math.factorial(k)*L[k]for k in range(72)),Poly([0]));assert p==Poly(list(map(Q,z['P'])))
assert sum((p[k]*math.factorial(k)*L[k]for k in range(72)),Poly([0]))==q
print('EXACT FACTORS, FOURIER TRANSFORM AND ZERO NORMALIZATION PASS')
