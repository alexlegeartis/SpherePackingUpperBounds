#!/usr/bin/env python3
"""Independent exact scalar verification; Python standard library only.
Uses a sign-preserving primitive pseudo-remainder Sturm sequence.
"""
from fractions import Fraction as Q
from pathlib import Path
from functools import reduce
import json,math,time
ROOT=Path(__file__).resolve().parent

def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a

def primitive(a):
 a=trim(a);g=reduce(math.gcd,(abs(x) for x in a),0)
 return [x//g for x in a] if g else [0]

def integers(a):
 den=math.lcm(*(x.denominator for x in a))
 return primitive([x.numerator*(den//x.denominator) for x in a])

def negative_remainder(a,b):
 r=a.copy();n=len(b)-1;lead=b[-1];s=1 if lead>0 else -1
 while len(r)-1>=n and r!=[0]:
  shift=len(r)-len(b);c=r[-1];r=[abs(lead)*x for x in r]
  for j,v in enumerate(b):r[j+shift]-=s*c*v
  r=primitive(r)
 return primitive([-x for x in r])

def sturm(a):
 a=primitive(a);b=primitive([i*a[i] for i in range(1,len(a))]);seq=[a,b]
 while seq[-1]!=[0] and len(seq[-1])>1:
  r=negative_remainder(seq[-2],seq[-1])
  if r==[0]:break
  seq.append(r)
 return seq

def sign_at(p,x):
 if x is None:v=p[-1]
 else:
  u,q=x.numerator,x.denominator;n=len(p)-1
  v=sum(c*u**i*q**(n-i) for i,c in enumerate(p))
 return (v>0)-(v<0)

def variations(signs):
 a=[s for s in signs if s]
 return sum(x!=y for x,y in zip(a,a[1:]))

def verify():
 start=time.time();D=json.loads((ROOT/'scalar.json').read_text())
 f=[Q(a) for a in D['F']];g=[Q(a) for a in D['G']];n=len(f)-1
 transformed=[sum((f[k]*(-1)**j*math.factorial(k)*math.comb(k+1,k-j)/math.factorial(j) for k in range(j,n+1)),Q(0)) for j in range(n+1)]
 if transformed!=g:raise ValueError('The Fourier polynomial identity fails.')
 if not (g[0]==1 and g[-1]>0 and f[-1]<0):raise ValueError('Scalar normalization/sign test failed.')
 alpha=Q(D['alpha']);boundary=alpha*Q(D["cutoff"])**2
 sf=sturm(integers(f));sg=sturm(integers(g))
 results={}
 for label,seq,left in [('spatial_tail',sf,boundary),('Fourier_half_line',sg,Q(0))]:
  a=[sign_at(p,left) for p in seq];b=[sign_at(p,None) for p in seq]
  vl,vr=variations(a),variations(b)
  if vl-vr!=0:raise ValueError('A scalar polynomial has roots on a forbidden half-line.')
  results[label]={'signs_at_left':a,'signs_at_infinity':b,'variation_at_left':vl,'variation_at_infinity':vr,'root_count':vl-vr}
 if sign_at(integers(f),boundary)>=0:raise ValueError('Spatial tail sign failed.')
 print('PASS exact Gaussian Fourier-polynomial identity, G(0)=1.')
 print('PASS scalar spatial tail is strictly negative; Fourier polynomial is strictly positive.')
 for k,v in results.items():print(k,v)
 print('Standard-library verification seconds:',round(time.time()-start,3))
 (ROOT/'scalar_standard_verification.json').write_text(json.dumps(results,indent=2))
 return results
if __name__=='__main__':verify()
