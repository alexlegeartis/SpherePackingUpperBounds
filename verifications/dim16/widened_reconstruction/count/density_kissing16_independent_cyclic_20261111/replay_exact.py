"""Verify every leaf of the saved cyclic Bernstein proof in exact integer arithmetic.
The polynomial input is separately checked against the original factors by audit_cyclic.py.
This replay never uses floating point, Arb, or sampled points for acceptance.
"""
if not __debug__:
    raise RuntimeError("Certificate verification requires assertions; do not run Python with -O or -OO.")

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,math,time
from fractions import Fraction
import numpy as np
from flint import fmpq,fmpq_mat,fmpz,fmpz_mat
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'density_kissing16_reweight_frozen_20261111';start=time.monotonic()
proof=json.loads((ROOT/'proof_events.json').read_text());polyb=(ROOT/'polynomial.json').read_bytes();assert hashlib.sha256(polyb).hexdigest()==proof['polynomial_sha256']
for name,sha in proof['input_sha256'].items():assert hashlib.sha256((SRC/name).read_bytes()).hexdigest()==sha
poly=json.loads(polyb);assert poly['source_sha256']==proof['input_sha256'];D={tuple(e):fmpq(z) for e,z in poly['terms']};d=max(max(e) for e in D);n=d+1;assert d==16
cap=fmpq(501,1000);a=fmpq(-1);delta=cap+1
BT=fmpq_mat([[sum((fmpq(math.comb(e,j)*math.comb(i,j),math.comb(d,j))*a**(e-j)*delta**j for j in range(min(i,e)+1)),fmpq(0))for e in range(n)]for i in range(n)])
A=np.empty((n,n,n),dtype=object);A.fill(fmpq(0))
for e,z in D.items():A[e]=z
for axis in range(3):
 V=np.moveaxis(A,axis,0).reshape(n,n*n);M=BT*fmpq_mat(V.tolist());V=np.array([[M[i,j]for j in range(n*n)]for i in range(n)],dtype=object);A=np.moveaxis(V.reshape(n,n,n),0,axis)
P=[-1+delta*fmpq(i,d)for i in range(n)];S=[1-2*delta*fmpq(i,d)+delta**2*fmpq(i*(i-1),d*(d-1))for i in range(n)]
G=np.empty_like(A)
for i,j,k in np.ndindex(G.shape):G[i,j,k]=1+2*P[i]*P[j]*P[k]-S[i]-S[j]-S[k]
common=fmpz(1)
for v in list(A.flat)+list(G.flat):common=common.lcm(v.denominator)
for e in np.ndindex(A.shape):A[e]=A[e].numerator*(common//A[e].denominator);G[e]=G[e].numerator*(common//G[e].denominator)
# Both child coefficient tensors have common denominator old_denominator*2^d.
# These exact integer matrices implement the de Casteljau halves.
LEFT=fmpz_mat([[math.comb(i,j)*2**(d-i) if j<=i else 0 for j in range(n)]for i in range(n)])
RIGHT=fmpz_mat([[math.comb(d-i,j-i)*2**i if j>=i else 0 for j in range(n)]for i in range(n)])
def split_exact(X,axis):
 V=np.moveaxis(X,axis,0).reshape(n,n*n);M=fmpz_mat(V.tolist());out=[]
 for T in (LEFT,RIGHT):
  Z=T*M;Y=np.array([[Z[i,j]for j in range(n*n)]for i in range(n)],dtype=object);out.append(np.moveaxis(Y.reshape(n,n,n),0,axis))
 return out
stack=[(A,G,((Fraction(0),Fraction(1)),)*3,0)];counts={'H':0,'M':0,'G':0,'S':0,'split':0};volume=Fraction(0);maxdepth=0
for idx,event in enumerate(proof['events']):
 assert stack,'extra events after stack exhausted'
 H,G,box,depth=stack.pop();maxdepth=max(maxdepth,depth);kind=event['kind'];counts[kind]+=1
 if kind=='split':
  axis=event['axis'];assert axis in (0,1,2);hl,hr=split_exact(H,axis);gl,gr=split_exact(G,axis)
  a,b=box[axis];m=(a+b)/2;bx1=list(box);bx2=list(box);bx1[axis]=(a,m);bx2[axis]=(m,b)
  stack.extend([(hr,gr,tuple(bx2),depth+1),(hl,gl,tuple(bx1),depth+1)])
 else:
  if kind=='H':assert all(v<0 for v in H.flat),'nonnegative H Bernstein coefficient'
  elif kind=='G':assert all(v<0 for v in G.flat),'nonnegative Gram Bernstein coefficient'
  elif kind=='S':assert box[0][0]>box[1][1] or box[1][0]>box[2][1],'not outside sorted chamber'
  elif kind=='M':
   num,den=map(fmpz,event['lambda_rational']);assert num>=0 and den>0
   assert all(den*h+num*g<0 for h,g in zip(H.flat,G.flat)),'nonnegative H+lambda G coefficient'
  else:raise AssertionError('unknown event')
  volume+=math.prod(b-a for a,b in box)
 if (idx+1)%200==0:print('EXACT_REPLAY',idx+1,'seconds',time.monotonic()-start,flush=True)
assert not stack,'unfinished boxes';assert volume==1
out={'verified':True,'arithmetic':'Exact FLINT rational initialization and arbitrary-size integer dyadic subdivision. Every leaf inequality checked by integer signs. No floating point or interval arithmetic in replay acceptance.','events':len(proof['events']),'counts':counts,'max_depth':maxdepth,'covered_unit_cube_volume_exact':str(volume),'common_initial_denominator_decimal_digits':len(str(common)),'proof_events_sha256':hashlib.sha256((ROOT/'proof_events.json').read_bytes()).hexdigest(),'polynomial_sha256':proof['polynomial_sha256'],'input_sha256':proof['input_sha256'],'seconds':time.monotonic()-start}
(ROOT/'exact_replay_result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
