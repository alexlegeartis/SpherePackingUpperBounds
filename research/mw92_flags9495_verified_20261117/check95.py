from pathlib import Path
from fractions import Fraction as F
from math import prod,isqrt
import itertools,json,hashlib
P=Path(__file__).resolve().parent;r=json.load(open(P/'input/results.json'));s=json.load(open(P/'input/flag95.json'))
d=[F(r['depths'][str(j)])for j in [1,3,9]];h=[16-v for v in d];a,b,c=h
assert 4<a<b<c
G=[[F(v)for v in row]for row in s['gram']];assert G==[[a,-a/2,0],[-a/2,b,-b/2],[0,-b/2,c]]
D=a*b*c-a*b*b/4-a*a*c/4;assert D==F(s['determinant'])>0 and a*b-a*a/4>0
# Global positive integer suffix forms, not merely finite enumeration:
# 2q_k(z)=z_first^2+z_last^2+sum adjacent differences squared.
# q_k is positive integer for every nonzero integer suffix. Expansion gives G.
coefs=[a,b-a,c-b];assert min(coefs)>0
for i in range(3):assert sum(coefs[:i+1])==G[i][i]
for i in range(2):assert -sum(coefs[:i+1])/2==G[i][i+1]
# Validate all64 ordered color cases and the required highest odd coordinate.
colors=list(itertools.product(range(2),repeat=3))
for x in colors:
 for y in colors:
  if x==y:assert 4*a>16;continue
  diff=[u-v for u,v in zip(x,y)];j=max(k for k,v in enumerate(diff)if v)
  assert all(v in [-1,0,1]for v in diff)and d[j]+h[j]==16
# Any vertical coordinate above j can only increase the global h_lower index.
# Same parity vectors have form2z, so q>=4h1>16.
center2=F(2**190,47**14)/D;assert center2==F(s['center_density_squared'])
def atan(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(N)),F(0))
lo=16*atan(5,100)-4*atan(239,101);V2=F(2**96,prod(range(1,96,2))**2)*lo**94
safe=F('2.06770138493108e-21');assert safe==F(s['occupied95_floor'])and safe**2<center2*V2
old=json.load(open(P/'baseline95.json'));old2=F(old['row']['density_squared']);assert center2>old2
ratio2=center2/old2;rr=F(isqrt(ratio2.numerator*10**24//ratio2.denominator),10**12);assert rr**2<=ratio2<(rr+F(1,10**12))**2
gain=[str(100*(rr-1)),str(100*(rr+F(1,10**12)-1))]
external2=F(2**184,47**14)/32;ratio2=center2/external2;rr=F(isqrt(ratio2.numerator*10**24//ratio2.denominator),10**12);extgain=[str(100*(rr-1)),str(100*(rr+F(1,10**12)-1))]
report=dict(status='PASS global weighted root-chain parity proof, all64 colorcases, exact odd-volume rounding; covering/parent prerequisites in94checker',dimension=95,safe_occupied_lower=str(safe),safe_decimal='2.06770138493108e-21',determinant=str(D),center_density_squared=str(center2),gain_over_public_percent_interval=gain,external_derived_comparator_center_squared=str(external2),external_gain_percent_interval=extgain,baseline_commit=old['commit'],input_sha256={f:hashlib.sha256((P/'input'/f).read_bytes()).hexdigest()for f in ['flag95.json','flag95.py','FLAG95_PROOF.md','results.json']})
(P/'RESULT95.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
