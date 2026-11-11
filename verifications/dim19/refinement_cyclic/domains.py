"""Explicit harmonic formula and fully rational Bernstein domain audit."""
from pathlib import Path
import json,math,itertools
import numpy as np
from flint import fmpq as Q,fmpq_mat as M
assert __debug__
P=Path(__file__).resolve().parent
d=json.loads((P/'certificate.json').read_text()); n=19;K=5;L=16
assert (d['dimension'],d['radial_nodes'],d['angular_degree'])==(n,K,L)
from fractions import Fraction
r=Q(str(Fraction(d['exclusion_squared'])));R=Q(str(Fraction(d['truncation_squared'])));assert r==Q(19,8) and R==Q(57,20)
As=[[[Q(x)for x in row]for row in A]for A in d['matrices']]
cells=[(r+(R-r)*i/4,r+(R-r)*(i+1)/4)for i in range(4)]
def rising(a,k):
 z=Q(1)
 for i in range(k):z*=a+i
 return z
def accumulate(out,e,v):out[e]=out.get(e,Q(0))+v
harm=[]
for l in range(L+1):
 out={};norm=rising(Q(n-2),l)/math.factorial(l)
 for j in range(l//2+1):
  k=l-2*j;c=(-1)**j*rising(Q(n-2,2),l-j)/(math.factorial(j)*math.factorial(k)*norm*r**l)
  for a in range(k+1):
   for b in range(k-a+1):
    z=k-a-b
    accumulate(out,(a+j,b+j,z),c*math.factorial(k)/(math.factorial(a)*math.factorial(b)*math.factorial(z))*(-1)**z)
 harm.append({e:v for e,v in out.items()if v})
def kernel(i,j):
 a,b=cells[i];c,e=cells[j]
 wi={i:(b/(b-a),-1/(b-a)),i+1:(-a/(b-a),1/(b-a))};wj={j:(e/(e-c),-1/(e-c)),j+1:(-c/(e-c),1/(e-c))}
 out={}
 for l,A in enumerate(As):
  off=int(l==0)
  for ii,ui in wi.items():
   for jj,vj in wj.items():
    for aa,bb in itertools.product(range(2),repeat=2):
     coef=A[ii+off][jj+off]*ui[aa]*vj[bb]
     for ex,value in harm[l].items():accumulate(out,(ex[0]+aa,ex[1]+bb,ex[2]),coef*value)
 return {e:v for e,v in out.items()if v}
def transform_axis(A,T,axis):
 V=np.moveaxis(A,axis,0);shape=V.shape;rows=shape[0];cols=math.prod(shape[1:]);prod=T*M(V.reshape(rows,cols).tolist())
 return np.moveaxis(np.array([[prod[i,j]for j in range(cols)]for i in range(rows)],dtype=object).reshape(shape),0,axis)
def coefficients(poly,box):
 degrees=[max(e[k]for e in poly)for k in range(3)]; A=np.full(tuple(d+1 for d in degrees),Q(0),dtype=object)
 for e,v in poly.items():A[e]=v
 for axis,(lo,hi)in enumerate(box):
  d=degrees[axis]
  # Separate affine monomial substitution, then power-to-Bernstein conversion.
  T=M([[Q(math.comb(k,j))*lo**(k-j)*(hi-lo)**j if k>=j else Q(0)for k in range(d+1)]for j in range(d+1)])
  A=transform_axis(A,T,axis)
  T=M([[Q(math.comb(i,j),math.comb(d,j))if i>=j else Q(0)for j in range(d+1)]for i in range(d+1)])
  A=transform_axis(A,T,axis)
 return A
def split(A,axis):
 a=np.moveaxis(A,axis,0).copy();d=a.shape[0]-1;left=np.empty_like(a);right=np.empty_like(a);left[0]=a[0];right[d]=a[-1]
 for j in range(1,d+1):
  a=(a[:-1]+a[1:])/Q(2);left[j]=a[0];right[d-j]=a[-1]
 return np.moveaxis(left,0,axis),np.moveaxis(right,0,axis)
def check(poly,box):
 stack=[(coefficients(poly,box),0)];nodes=0;depth=0
 while stack:
  A,k=stack.pop();nodes+=1;depth=max(depth,k)
  if all(x<0 for x in A.flat):continue
  assert k<30 and nodes<20000
  # Float only chooses subdivision axis, never certifies a sign or excludes a region.
  axis=max(range(3),key=lambda j:max((abs(float(x))for x in np.diff(A,axis=j).flat),default=0))
  a,b=split(A,axis);stack.extend([(a,k+1),(b,k+1)])
 return {'nodes':nodes,'depth':depth}
cache={(i,j):kernel(i,j)for i in range(4)for j in range(i,4)};results=[]
for i,j,k in itertools.combinations_with_replacement(range(4),3):
 out=cache[i,j].copy()
 for e,v in cache[i,k].items():accumulate(out,(e[0],e[2],e[1]),v)
 for e,v in cache[j,k].items():accumulate(out,(e[2],e[0],e[1]),v)
 result=check(out,[cells[i],cells[j],cells[k]]);results.append(['cyclic',i,j,k,result]);print(results[-1],flush=True)
for i,j in itertools.combinations_with_replacement(range(4),2):
 result=check(cache[i,j],[cells[i],cells[j],(R,4*R)]);results.append(['far',i,j,result]);print(results[-1],flush=True)
(P/'domains_audit.json').write_text(json.dumps(results,indent=2)+'\n');print('PASS all cyclic and far domains, fully exact rational arithmetic, no geometry exclusions')
