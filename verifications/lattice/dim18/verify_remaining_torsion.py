if not __debug__: raise RuntimeError('Do not run with -O')
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import json, math
out={};exceptional=[];other=F(0)
for v in combinations_with_replacement(range(3),18):
 if not any(v):continue
 costs=[]
 for k in [1,2]:
  rr=[min(k*x%5,5-k*x%5) for x in v]
  cost=F(sum(x*x for x in rr),50)
  if sum(rr)%2:cost+=F(25-10*max(rr),50)
  costs.append(cost)
 cost=min(costs)
 if cost>F(22,25):exceptional.append((v,cost))
 else:other=max(other,cost)
assert len(exceptional)==1 and exceptional[0][0]==(1,)*9+(2,)*9 and exceptional[0][1]==1
assert other==F(22,25)
rows=[]
for p in range(13,131):
 if any(p%d==0 for d in range(2,math.isqrt(p)+1)):continue
 residues=lambda a:[min(k*a%p,p-k*a%p) for k in range(1,p)]
 m=(p-1);r1=residues(1)
 avmax=min(F(sum(max(a,b) for a,b in zip(r1,residues(j))),m*p) for j in range(2,(p+1)//2))
 corr=F(sum(p-2*r for r in r1 if r%2),m*p)
 assert corr==F(p+1,4*p)
 bound=max(F(18*(p+1),12*p),F(17*(p+1),12*p)+corr,F(18*(p+1),12*p)+1-2*avmax)/2
 assert bound<=F(37,39)
 rows.append({'p':p,'physical_witness':str(bound)})
# Unit ratios for order25, both signs included independently.
un=[k for k in range(1,25) if math.gcd(k,25)==1]
r=lambda k:min(k%25,25-k%25)
mean=F(sum(r(k)**2 for k in un),len(un)*625)
assert mean==F(21,250)
maxbound=min(F(sum(max(r(k),r(a*k)) for k in un),len(un)*25) for a in un if a not in [1,24])
assert maxbound==F(42,125)
assert (18*mean+1-2*maxbound)/2==F(23,25)
# Single absolute unit-type case scales to18 copies of±1/25, parity even.
assert F(18,2*625)<F(23,25)
assert 24*18>18*20
assert F(18*6,8)<14
assert max(F(w+(3 if w%2 else 0),18) for w in range(1,14))==F(8,9)
S=F(1,25);dist=F(421,400)
assert F(17,4)**2>18 and dist==1+F(21,16)*S
for witness in [F(43,49),F(22,27),F(95,121),F(37,39),F(22,25),F(23,25),F(8,9)]:assert dist*witness<1
assert S<F(8,17)
assert F(33100,256)*(1+S/36)**9<131
allowed=[i for i in range(1,131) if i%9 and i%25 and all(i%r['p'] for r in rows) and i%7 and i%11]
assert max(allowed)==128 and 128*256==32768
def atan(q,N):return sum(((-1)**j*F(1,(2*j+1)*q**(2*j+1))for j in range(N)),F(0))
pi_lower=16*atan(5,20)-4*atan(239,21)
assert pi_lower>F('3.14159')
prior=json.loads(Path(__file__).with_name('prior_certificate.json').read_text())
assert F(prior['claims']['density'])<=F('0.0103704834')
assert F('0.0103704834')*2**18*math.factorial(9)/pi_lower**9<33100
out={'passed':True,'large_primes':rows,'p5_exception':'only (9,9) fullsupport','order25_witness':'23/25','index_upper':131,'allowed_index_max':128,'conditional_rho':32768,'gap':'1/25'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2));print('INDEPENDENT18D REMAINING CHECKS PASSED')
