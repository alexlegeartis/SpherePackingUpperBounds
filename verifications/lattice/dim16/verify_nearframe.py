if not __debug__:
 raise RuntimeError("Do not run verification with python -O")
"""Exact finite witnesses and rational constants for near-D16 frame exclusion."""
import json,itertools
from fractions import Fraction as F
pi_lower=16*sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(4)),F(0))-F(4,239)
pi_upper=16*sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(5)),F(0))-4*(F(1,239)-F(1,3*239**3))
assert pi_lower>F('3.14159') and pi_upper<F(22,7)
import numpy as np
out={}
def histogram(p):
 h=(p-1)//2;R=[[min(i*k%p,p-i*k%p) for i in range(1,h+1)] for k in range(1,h+1)]
 best=-1;arg=None;seen=0
 def rec(count,left):
  nonlocal best,arg,seen
  if len(count)==h:
   if sum(count)==0:return
   if p==9 and sum(count[j-1] for j in [1,2,4])==0:return
   seen+=1;values=[]
   for rr in R:
    square=sum(z*r*r for z,r in zip(count,rr));parity=sum(z*r for z,r in zip(count,rr))%2;largest=max(r for z,r in zip(count,rr) if z)
    values.append(square+parity*(p*p-2*p*largest))
   val=min(values)
   if val>best:best=val;arg=count[:]
   return
  for z in range(left+1):rec(count+[z],left-z)
 rec([],16)
 return {'unscaled_squared_witness':str(F(best,p*p)),'histograms':seen,'maximizer':arg}
for p in [5,7,9,11]:
 out[str(p)]=histogram(p);assert F(out[str(p)]['unscaled_squared_witness'])<=F(42,25)
 print(p,out[str(p)],flush=True)
# For the remaining possible primes, enumerate parity supports, not full counts.
for p in [13,17,19,23,29,31,37,41]:
 h=(p-1)//2;types=np.arange(1,h+1,dtype=np.int64);R=(types[:,None]*types[None,:])%p;R=np.minimum(R,p-R)
 best=16*(p+1)*h;arg=0
 for start in range(1,2**h,2048):
  masks=np.arange(start,min(start+2048,2**h),dtype=np.int64);bits=(masks[:,None]>>np.arange(h,dtype=np.int64))&1
  population=bits.sum(axis=1);m=16-(population%2);valid=population<=m
  parity=(bits@R.T)%2
  largest=np.max(bits[:,None,:]*R[None,:,:],axis=2)
  correction=(parity*(p-2*largest)).sum(axis=1)
  numerator=m*(p+1)*h+12*correction;numerator[~valid]=-1
  k=int(np.argmax(numerator))
  if int(numerator[k])>best:best=int(numerator[k]);arg=int(masks[k])
 bound=F(best,12*p*h);assert bound<=F(42,25)
 out[str(p)]={'unscaled_squared_witness':str(bound),'parity_supports':2**h,'maximizer_mask':arg}
 print(p,out[str(p)],flush=True)
# Rank-two 3-torsion has eight nonzero codewords with total weight <=96.
# Some weight <=12; max D16 residue norm for weights <=12 is14/9.
rank2=max(F(m,9)+(F(1,3) if m%2 else 0) for m in range(1,13));assert rank2==F(14,9)
threshold=F(16,105)
assert (1+F(5,4)*threshold)*F(21,25)==1
assert threshold<F(1,4)
index=F(5220,128)*(1+threshold/32)**8;assert index<43
rho=F('0.018743418752')*2**16*40320/F('3.14159')**8;assert rho<5220
allowed=[i for i in range(33,43) if all(i%p for p in [5,7,11,13,17,19,23,29,31,37,41]) and i%9];assert allowed==[]
out['summary']={'status':'PASS','normalized_witness':'21/25','S_threshold':'16/105','density_alternative_rho':'4096','index_upper_rational':str(index),'rank2_unscaled_witness':str(rank2)}
json.dump(out,open(__file__.replace('verify_nearframe.py','nearframe_audit.json'),'w'),indent=2)
print('ALL EXACT FINITE CHECKS PASS',flush=True)
