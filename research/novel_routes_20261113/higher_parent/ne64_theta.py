"""Exact formal series for the standard 3-modular theta basis.
The use of this theta-basis theorem is an external mathematical dependency.
"""
from pathlib import Path
import math,json,re
import sympy as sp
P=Path(__file__).resolve().parent;N=16
# Source-format, discriminant, and exponent checks are independent of the
# externally supplied 3-modularity and minimum theorem.
source=(P/'CQ64b.html').read_text()
raw=source[source.index('<a NAME="GRAM"'):source.index('<a NAME="DIVISORS"')]
vals=[int(x) for x in re.findall(r'-?\d+',re.sub('<[^>]*>',' ',raw))]
assert vals[:2]==[64,0] and len(vals)==2082
G=sp.zeros(64);cursor=2
for i in range(64):
 for j in range(i+1):G[i,j]=G[j,i]=vals[cursor];cursor+=1
L,D=G.LDLdecomposition(hermitian=False)
assert all(D[i,i]>0 for i in range(64))
assert sp.prod(D[i,i] for i in range(64))==3**32
assert all(G[i,i]%2==0 for i in range(64))
assert all(x.q==1 for x in 3*G.inv())

def mul(a,b):
 c=[0]*(N+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b[:N+1-i]):c[i+j]+=x*y
 return c

def power(a,n):
 b=[1]+[0]*N
 for _ in range(n):b=mul(b,a)
 return b
A=[0]*(N+1)
for x in range(-10,11):
 for y in range(-10,11):
  q=x*x-x*y+y*y
  if q<=N:A[q]+=1
F=[1]+[0]*N
for k in range(1,N+1):
 for d in (k,3*k):
  f=[0]*(N+1)
  for j in range(7):
   if j*d<=N:f[j*d]=(-1)**j*math.comb(6,j)
  F=mul(F,f)
F=[0]+F[:-1]
basis=[mul(power(A,32-6*j),power(F,j)) for j in range(6)]
coefs=[];theta=[0]*(N+1)
for j,b in enumerate(basis):
 c=1 if j==0 else-theta[j];assert b[j]==1;coefs.append(c)
 theta=[x+c*y for x,y in zip(theta,b)]
assert theta[:6]==[1,0,0,0,0,0];assert theta[6]==138458880
iso=3**31-2*3**15-1;short=theta[6]+theta[9]
assert short<iso
report={'source_gram_dimension':64,'source_gram_determinant':str(3**32),'positive_definite_exact':True,'discriminant_exponent3_exact':True,'basis':'sum c_j theta_A2^(32-6j) (eta(tau)eta(3tau))^(6j), j=0,...,5',
'coefficients':coefs,'theta_coeff_q_norm_half':theta,'isotropic_nonzero_discriminant_classes_lower_bound':iso,
'dual_norm4_vectors':theta[6],'dual_norm6_vectors':theta[9],
'isotropic_classes_occupied_by_short_vectors_at_most':short,
'isotropic_classes_with_coset_min_at_least8_at_least':iso-short,
'centerdensity_65_squared':str(3**33)+'/4','centerdensity_65_approx':3**16*math.sqrt(3)/2,
'dependency':'The standard strongly3-modular theta-series ring theorem; published Ne64 minimum12,det3^32,and3-modularity.'}
(P/'ne64_theta.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
