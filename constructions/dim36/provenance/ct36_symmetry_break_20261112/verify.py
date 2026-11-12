from pathlib import Path
from fractions import Fraction as F
import json,numpy as np,math
P=Path(__file__).resolve().parent
def K(n,k,x):return sum((-1)**j*3**(k-j)*math.comb(x,j)*math.comb(n-x,k-j)for j in range(max(0,k-(n-x)),min(k,x)+1))
for r in json.load(open(P/'code_sweep.json')):
 n=r['length'];f=list(map(F,r['coefficients']));assert all(x>=0 for x in f)
 assert all(1+sum(f[k-1]*K(n,k,i)for k in range(1,n+1))<=0 for i in range(8,n+1))
 assert 1+sum(f[k-1]*K(n,k,0)for k in range(1,n+1))==F(r['code_upper_rational'])
N=np.load(P/'reduced.npz')['N'];J=json.load(open(P/'holes.json'));ref={(r['p'],tuple(r['a'])):r for r in json.load(open(P/'refined_holes.json'))}
for r in J['rows']:
 p=r['p'];a=np.array(r['a'],dtype=np.int64);s=ref.get((p,tuple(r['a'])),r);z=np.array(s.get('refined_witness',s['witness']),dtype=np.int64);y=a-p*z;q=F(int(y@N@y),4*p*p)
 assert q==F(s.get('refined_q',s['q']));assert F(24576,15625)*(min(q,F(8))/8)**18<=1
print('PASS12 rational LP certificates and900 exact two-coset rejection witnesses')
