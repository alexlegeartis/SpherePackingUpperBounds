from fractions import Fraction as F
from pathlib import Path
import json, math
from flint import arb,ctx
assert __debug__
ctx.prec=256
P=Path(__file__).resolve().parent
sources=P.parent
import hashlib
for name,digest in json.loads((P/'source_inputs.json').read_text()).items():
 assert hashlib.sha256((P.parent.parent/name).read_bytes()).hexdigest()==digest, name
registry={r['dimension']:r for r in json.loads((P.parent.parent/'RESULTS.json').read_text())}
def ar(q):return arb(q.numerator)/arb(q.denominator)
def tighten(Q):
 B=[r[1:] for r in Q[1:]];b=[r[0] for r in Q[1:]];n=len(B)
 L=[[F(i==j) for j in range(n)] for i in range(n)];d=[]
 for j in range(n):
  d.append(B[j][j]-sum(L[j][k]**2*d[k] for k in range(j)))
  assert d[-1]>0
  for i in range(j+1,n):
   L[i][j]=(B[i][j]-sum(L[i][k]*L[j][k]*d[k] for k in range(j)))/d[j]
 assert all(B[i][j]==sum(L[i][k]*d[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
 y=[]
 for i in range(n): y.append(b[i]-sum(L[i][j]*y[j] for j in range(i)))
 z=[y[i]/d[i] for i in range(n)];x=[F(0)]*n
 for i in reversed(range(n)):x[i]=z[i]-sum(L[j][i]*x[j] for j in range(i+1,n))
 assert all(sum(B[i][j]*x[j] for j in range(n))==b[i] for i in range(n))
 alpha=Q[0][0]-sum(b[i]*x[i] for i in range(n));assert alpha>0
 W=[x]+[[F(i==j) for j in range(n)]for i in range(n)]
 assert all(Q[i][j]-(alpha if i==j==0 else 0)==sum(W[i][a]*B[a][b]*W[j][b] for a in range(n) for b in range(n)) for i in range(n+1) for j in range(n+1))
 return alpha
out=[]
producer=json.loads((P/'schur_records.json').read_text())
for n in [4,5,19,20,21]:
 if n<6:
  src=sources/('dim4/refinement_p6' if n==4 else 'dim5/refinement_p4')
  raw=json.loads((src/'witness.json').read_text(),parse_float=str)
  Q=[[F(x) for x in r]for r in json.loads((src/'witness_psd_exact.json').read_text())['matrices'][0]]
  c=list(map(F,raw['x'][:32]));H=sum((-1)**k*x for k,x in enumerate(c));obj=(sum(c)+Q[0][0]+sum(F(a)*F(b) for a,b in zip(raw['count_weights'],raw['count_bounds'])))/H;r2=F(1)
 else:
  raw=json.loads((sources/f'dim{n}'/'certificate.json').read_text(),parse_float=str)
  Q=[[F(x)for x in r]for r in raw['matrices'][0]];H=F(raw['normalization']);r2=F(raw['exclusion_squared'])
  assert H==F(raw['fourier_polynomial'][0])
  obj=(F(raw['physical_polynomial'][0])+Q[0][0])/H
  assert obj==F(raw['exact_objective'])
 alpha=tighten(Q); record=next(x for x in producer if x['dimension']==n)
 assert alpha==F(record['alpha_exact']) and H==F(record['normalization_exact'])
 if n%2==0:
  k=n//2;volume=arb.pi()**k/math.factorial(k)*ar(r2/4)**k
 else:
  k=(n-1)//2;volume=2**(2*k+1)*math.factorial(k)*arb.pi()**k/math.factorial(2*k+1)*ar(r2/4)**k*ar(r2/4).sqrt()
 density=volume*ar(obj-alpha/H);decrease=volume*ar(alpha/H)
 assert density<ar(F(record['new_upper_exact']))
 sharp={4:'0.633428734445',5:'0.510501094478',19:'0.009618343139',20:'0.006979497272',21:'0.005149549312'}[n]
 assert density<ar(F(sharp))
 assert json.loads((P/'certified_bounds.json').read_text())[str(n)]==sharp
 assert registry[n]['strict_density_upper_bound']==sharp and registry[n]['package']=='verifications/origin_slack'
 assert decrease>ar(F(record['density_decrease_lower_exact']))
 item={'dimension':n,'density_interval':str(density),'decrease_interval':str(decrease),'producer_upper_verified':record['conditional_outward_upper'],'sharper_upper_verified':sharp};out.append(item);print(json.dumps(item),flush=True)
(P/'independent_all_results.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS independent Fraction algebra and Arb enclosures in all five dimensions')
