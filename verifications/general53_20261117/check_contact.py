from pathlib import Path
from fractions import Fraction as F
import json,math,sys
from flint import fmpq as Q,fmpq_poly as Poly
sys.set_int_max_str_digits(0);P=Path(__file__).resolve().parent
for n in [53]:
 D=P/f'dim{n}';c=json.loads((D/'contact_cap_certificate.json').read_text());co=list(map(Q,c['gegenbauer']));assert co[0]==1 and all(a>=0 for a in co)
 lam=Q(n-2,2)
 def rise(a,k):
  z=Q(1)
  for i in range(k):z*=a+i
  return z
 f=Poly([0])
 for k,a in enumerate(co):
  g=[Q(0)]*(k+1)
  for j in range(k//2+1):g[k-2*j]=(-1)**j*rise(lam,k-j)*2**(k-2*j)*math.factorial(k)/(math.factorial(j)*math.factorial(k-2*j)*rise(2*lam,k))
  f+=a*Poly(g)
 assert f==Poly(list(map(Q,c['polynomial'])))
 s=Q(c['s']);assert Q(c['epsilon'])==Q(1,1000) and s==Q(501,1001)
 st=[f,f.derivative()]
 while st[-1]!=0:
  rem=-(st[-2]%st[-1])
  if rem==0:break
  st.append(rem)
 def variations(t):
  signs=[1 if v(t)>0 else -1 for v in st if v(t)!=0]
  return sum(a!=b for a,b in zip(signs,signs[1:]))
 assert f(-1)<0 and f(s)<0
 assert variations(Q(-1))-variations(s)==0
 assert f(1)==sum(co) and c['cap']==int(f(1).numerator//f(1).denominator)
 print(n,'PASS independent closed Gegenbauer reconstruction, exact Sturm nonpositivity and neighbor cap',c['cap'],flush=True)
