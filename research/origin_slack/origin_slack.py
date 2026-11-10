"""Exact origin-coordinate Schur tightening. Source continuum replay is separate."""
import json,math
from pathlib import Path
import sympy as s
P=Path(__file__).resolve().parent
assert __debug__
def read(name):return json.loads((P/name).read_text(),parse_float=str)
def matrix(a):return s.Matrix([[s.Rational(x) for x in row] for row in a])
def positive(A):
 L,D=A.LDLdecomposition(hermitian=False)
 assert L*D*L.T==A and all(x>0 for x in D.diagonal())
def sqrt_lower(q):
 q=s.Rational(q); scale=10**30
 a=math.isqrt(int(q.p)*scale**2//int(q.q))
 r=s.Rational(a,scale)
 assert r*r<=q and (r+1/s.Integer(scale))**2>q
 return r
out=[]
for n in [4,5,19,20,21]:
 if n<6:
  psd=read('dim4_psd.json' if n==4 else 'witness_psd_exact.json')
  raw=read('dim4_witness.json' if n==4 else 'witness.json')
  Q=matrix(psd['matrices'][0]); c=list(map(s.Rational,raw['x'][:32])); H=sum((-1)**k*t for k,t in enumerate(c))
  old=s.Rational('0.633430327' if n==4 else '0.510522525'); r2=s.Integer(1)
 else:
  raw=read(f'dim{n}_certificate.json'); Q=matrix(raw['matrices'][0]); H=s.Rational(raw['normalization']); old=s.Rational(raw['certified_density_upper']);r2=s.Rational(raw['exclusion_squared'])
 assert H>0
 positive(Q)
 B=Q[1:,1:];b=Q[1:,0];positive(B)
 alpha=Q[0,0]-(b.T*B.inv()*b)[0]
 assert alpha>0 and alpha==1/Q.inv()[0,0]
 e=s.zeros(Q.rows,1);e[0]=1
 W=(B.inv()*b).T.col_join(s.eye(B.rows))
 assert W*B*W.T==Q-alpha*e*e.T
 # Volume v_n (r/2)^n, with rational lower enclosure.
 pi=s.Rational(223,71)
 if n%2==0:
  volume=pi**(n//2)/s.factorial(n//2)*(r2/4)**(n//2)
 else:
  k=(n-1)//2
  volume=2**(2*k+1)*s.factorial(k)*pi**k/s.factorial(2*k+1)*(r2/4)**k*sqrt_lower(r2/4)
 decrease_lower=alpha/H*volume
 exact_volume=s.pi**s.Rational(n,2)/s.gamma(s.Rational(n,2)+1)*(r2/4)**s.Rational(n,2)
 digits=12
 guaranteed_decrease=s.floor(decrease_lower*10**digits)/s.Integer(10**digits)
 assert guaranteed_decrease>0 and guaranteed_decrease<=decrease_lower
 new=old-guaranteed_decrease
 result={'dimension':n,'alpha_exact':str(alpha),'alpha_approx':str(s.N(alpha,20)),'normalization_exact':str(H),'density_decrease_approx':str(s.N(alpha/H*exact_volume,25)),'density_decrease_lower_exact':str(decrease_lower),'old_certified_decimal':str(old),'conditional_outward_upper':format(float(new),'.12f'),'new_upper_exact':str(new),'PSD_factorization_exact':True,'continuum_source_replayed_here':False}
 out.append(result)
 print(n,result['alpha_approx'],result['density_decrease_approx'],result['conditional_outward_upper'])
(P/'origin_slack_results.json').write_text(json.dumps(out,indent=2)+'\n')
