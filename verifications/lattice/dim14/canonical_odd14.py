"""Exact finite arithmetic for the conditional 14D frame gap; no optimizer."""
import json, math
from fractions import Fraction as F

def histogram(n,p):
 h=(p-1)//2
 R=[[min(i*k%p,p-i*k%p) for i in range(1,h+1)] for k in range(1,h+1)]
 best=-1;arg=None;seen=0
 def rec(c,left):
  nonlocal best,arg,seen
  if len(c)==h:
   if not sum(c):return
   if p==9 and not sum(c[j-1] for j in [1,2,4]):return
   seen+=1;values=[]
   for rr in R:
    square=sum(z*r*r for z,r in zip(c,rr));parity=sum(z*r for z,r in zip(c,rr))%2
    mx=max(r for z,r in zip(c,rr) if z)
    values.append(square+parity*(p*p-2*p*mx))
   value=min(values)
   if value>best:best=value;arg=c[:]
   return
  for z in range(left+1):rec(c+[z],left-z)
 rec([],n)
 return {'physical_squared_witness':str(F(best,2*p*p)),'cases':seen,'histogram':arg}
