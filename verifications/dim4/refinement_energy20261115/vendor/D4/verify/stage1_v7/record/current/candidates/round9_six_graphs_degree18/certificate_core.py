"""Exact D4 two-annulus certificate checker: all decisions use integers/rationals.
Positive integer Gram rays ensure PSD. Bernstein subdivision is exhaustive,
not sampling. UNKNOWN and VIOLATION never constitute a proof.
"""
import json,math,time,itertools
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
import numpy as np
ROOT=Path(__file__).resolve().parent
BREAKS=tuple(map(Q,json.loads((ROOT/'kernel.json').read_text())['breaks']))
NS=len(BREAKS)-1
DEG=int(json.loads((ROOT/"kernel.json").read_text())["local_degree"])
NC=DEG+1
if not 1<=DEG<=14: raise ValueError("Unsupported radial degree.")
if BREAKS[0]!=1 or BREAKS[-1]>2 or any(a>=b for a,b in zip(BREAKS,BREAKS[1:])):raise ValueError('Unsupported annulus intervals.')

def conv(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):c[i+j]+=u*v
 return c

def gram():
 x=json.loads((ROOT/'kernel.json').read_text());L=x['angular_degree'];den=int(x['common_denominator'])
 if x['dimension']!=4 or x['local_degree']!=DEG or list(map(Q,x['breaks']))!=list(BREAKS):raise ValueError('Wrong kernel parameters.')
 if len(x['blocks'])!=L+1 or den<=0:raise ValueError('Invalid Gram description.')
 mats=[]
 for ell,block in enumerate(x['blocks']):
  n=NC*NS+1 if ell==0 else NC*NS
  if block['dimension']!=n:raise ValueError('Wrong block size.')
  a=np.zeros((n,n),dtype=object)
  for ray in block['rays']:
   w=ray['weight'];v=ray['vector']
   if not isinstance(w,int) or w<=0 or len(v)!=n or any(not isinstance(z,int) for z in v):raise ValueError('Invalid Gram ray.')
   for i,u in enumerate(v):
    for j,z in enumerate(v):a[i,j]+=w*u*z
  mats.append(a)
 return mats,den,L

def radial(s):
 spec=json.loads((ROOT/'kernel.json').read_text()).get('radial_basis_intervals')
 left,right=tuple(map(Q,spec[s])) if spec is not None else BREAKS[s:s+2]
 assert left<=BREAKS[s]<BREAKS[s+1]<=right
 width=right-left;start=[-left/width,1/width];end=[right/width,-1/width]
 if any(x.denominator!=1 for x in start+end):raise ValueError('Nonintegral radial coefficients need explicit denominators.')
 start=list(map(int,start));end=list(map(int,end));rows=[]
 for i in range(NC):
  p=[1]
  for _ in range(i):p=conv(p,start)
  for _ in range(DEG-i):p=conv(p,end)
  rows.append([math.comb(DEG,i)*v for v in p])
 return np.array(rows,dtype=object)

def gegen(L):
 den=math.lcm(*range(1,L+2));rows=[]
 for ell in range(L+1):
  p=[0]*(ell+1)
  for j in range(ell//2+1):p[ell-2*j]=(-1)**j*math.comb(ell-j,j)*2**(ell-2*j)*(den//(ell+1))
  if sum(p)!=den:raise ValueError('Gegenbauer normalization failed.')
  rows.append(p)
 return rows,den

def _build_uncached(kind,regions,eta=Q(0)):
 As,aden,L=gram();P,pden=gegen(L);Bs=[radial(s) for s in range(NS)]
 def matrices(i,j):
  out=[]
  for ell,A in enumerate(As):
   M=A[1:,1:] if ell==0 else A
   out.append(Bs[i].T@M[NC*i:NC*i+NC,NC*j:NC*j+NC]@Bs[j])
  return out
 if kind=='full':
  C=np.zeros((3*L+1,)*3,dtype=object);den=aden*pden
  for ai,aj,ak in [(0,1,2),(0,2,1),(1,2,0)]:
   R=matrices(regions[ai],regions[aj])
   for ell,M in enumerate(R):
    for k,pk in enumerate(P[ell]):
     if not pk:continue
     fac=pk*2**(L-k)
     for i,j in itertools.product(range(NC),repeat=2):
      if not M[i,j]:continue
      pre=fac*M[i,j];a0=L-k+i;b0=L-k+j
      for u in range(k+1):
       for v in range(k-u+1):
        w=k-u-v;coef=pre*math.comb(k,u)*math.comb(k-u,v)*(-1)**w
        idx=[0,0,0];idx[ai]=a0+2*u;idx[aj]=b0+2*v;idx[ak]=L+2*w;C[tuple(idx)]+=coef
  idx=(L,L,L)
  bounds=[(BREAKS[s],BREAKS[s+1]) for s in regions]
 elif kind=='outside':
  rsq=BREAKS[-1]**2;RD=rsq.denominator
  C=np.zeros((2*L+NC,2*L+NC,L+1),dtype=object);den=aden*pden*RD**L
  W=[{(0,0,0):1}];terms=[((1,1,0),-2*RD),((2,0,1),RD),((1,1,1),2*RD),((0,2,1),RD),((0,0,1),-rsq.numerator)]
  for _ in range(L):
   d={}
   for a,v in W[-1].items():
    for b,z in terms:
     t=tuple(x+y for x,y in zip(a,b));d[t]=d.get(t,0)+v*z
   W.append({k:v for k,v in d.items() if v})
  for ell,M in enumerate(matrices(*regions)):
   for k,pk in enumerate(P[ell]):
    if not pk:continue
    fac=pk*(2*RD)**(L-k)
    for i,j in itertools.product(range(NC),repeat=2):
     if not M[i,j]:continue
     pre=fac*M[i,j]
     for (u,v,w),z in W[k].items():C[L-k+i+u,L-k+j+v,w]+=pre*z
  idx=(L,L,0);bounds=[(BREAKS[s],BREAKS[s+1]) for s in regions]+[(Q(0),Q(1))]
 else:raise ValueError(kind)
 dd=math.lcm(den,eta.denominator)
 if dd!=den:C*=dd//den
 C[idx]-=eta.numerator*(dd//eta.denominator)*2**L
 return C,dd,bounds

@lru_cache(maxsize=8)
def _base_build(kind,regions):
 return _build_uncached(kind,regions,Q(0))

def build(kind,regions,eta=Q(0)):
 # Cache only exact unmodified numerator data; a tolerance shift always copies.
 C,den,bounds=_base_build(kind,tuple(regions))
 if eta==0:return C,den,bounds
 _,_,L=gram();dd=math.lcm(den,eta.denominator)
 shifted=C.copy()*(dd//den)
 idx=(L,L,L) if kind=='full' else (L,L,0)
 shifted[idx]-=eta.numerator*(dd//eta.denominator)*2**L
 return shifted,dd,bounds

# The polynomial construction and subdivision remain unchanged. This local
# exact-integer backend is loaded by path so a previously imported helper
# from another candidate cannot be silently reused. Its source is included
# in the global verifier's beginning/end hash contract.
import importlib.util as _exact_importlib
_exact_spec = _exact_importlib.spec_from_file_location(
    '_candidate_exact_bernstein_backend', ROOT/'exact_bernstein_fast.py')
_exact_backend = _exact_importlib.module_from_spec(_exact_spec)
_exact_spec.loader.exec_module(_exact_backend)
transform = _exact_backend.transform
to_bernstein = _exact_backend.to_bernstein

def split(A,axis,upper):
 w=np.moveaxis(A,axis,0).copy();n=w.shape[0]-1;left=np.empty_like(w);right=np.empty_like(w);left[0]=w[0];right[n]=w[n]
 for j in range(1,n+1):
  w=(w[:-1]+w[1:]+int(upper))//2;left[j]=w[0];right[n-j]=w[-1]
 return np.moveaxis(left,0,axis),np.moveaxis(right,0,axis)

def certify(C,den,bounds,collect=False,node_limit=120000,depth_limit=54):
 start=time.time();C,den=to_bernstein(C,den,bounds);big=max(abs(int(v)) for v in C.flat)
 p=min(50,58-(big.bit_length()-den.bit_length()+1))
 if p<0:return {'status':'UNKNOWN','reason':'coefficient dynamic range'}
 low=np.array([int(v)*2**p//den for v in C.flat],dtype=np.int64).reshape(C.shape)
 high=np.array([-((-int(v)*2**p)//den) for v in C.flat],dtype=np.int64).reshape(C.shape)
 if max(abs(low).max(),abs(high).max())>=2**59:raise ValueError('Integer overflow guard failed.')
 dims=C.ndim;del C;stack=[(low,high,((0,1),)*dims,(0,)*dims)];nodes=leaves=0;maxdepth=0;fails=[];unknown=False
 while stack:
  lo,hi,box,lev=stack.pop();nodes+=1;maxdepth=max(maxdepth,sum(lev))
  if hi.max()<=0:leaves+=1;continue
  fail=None
  for corner in itertools.product((0,1),repeat=dims):
   ix=tuple(-1 if c else 0 for c in corner)
   if lo[ix]>0:
    z=[Q(box[i][corner[i]],2**lev[i]) for i in range(dims)]
    point=[bounds[i][0]+z[i]*(bounds[i][1]-bounds[i][0]) for i in range(dims)]
    fail={'point':list(map(str,point)),'integer_lower':int(lo[ix]),'binary_scale':p};break
  if fail:
   fails.append(fail)
   if not collect or len(fails)>=100:break
   continue
  if nodes>=node_limit or sum(lev)>=depth_limit:unknown=True;break
  axis=min(range(dims),key=lambda i:lev[i])
  ll,lr=split(lo,axis,False);ul,ur=split(hi,axis,True);v=list(lev);v[axis]+=1;v=tuple(v)
  bl=list(box);br=list(box);a,b=box[axis];bl[axis]=(2*a,a+b);br[axis]=(a+b,2*b)
  stack.extend([(lr,ur,tuple(br),v),(ll,ul,tuple(bl),v)])
 status='VIOLATION' if fails else ('UNKNOWN' if unknown or stack else 'CERTIFIED')
 return dict(status=status,nodes=nodes,leaves=leaves,maxdepth=maxdepth,failures=fails,binary_scale=p,seconds=time.time()-start)

def check(kind,regions,eta=Q(0),collect=False):
 C,den,bounds=build(kind,regions,eta);res=certify(C,den,bounds,collect=collect);res.update(kind=kind,regions=list(regions),eta=str(eta))
 print('EXACT',kind,regions,{k:v for k,v in res.items() if k!='failures'},'violations',len(res['failures']),flush=True)
 return res

def pair_checks(origin_vector=None):
 As,den,L=gram();S=json.loads((ROOT/'scalar.json').read_text());fs=list(map(Q,S['F']));alpha=Q(S['alpha']);out=[]
 def add(a,b):
  c=[Q(0)]*max(len(a),len(b))
  for i,v in enumerate(a):c[i]+=v
  for i,v in enumerate(b):c[i]+=v
  return c
 for region in range(NS):
  B=radial(region);j=NC*region
  v=([Q(int(z),den) for z in As[0][0,1+j:1+j+NC]] if origin_vector is None else list(origin_vector[j:j+NC]))
  if len(v)!=NC: raise ValueError('Origin vector size mismatch.')
  coupling=[sum((v[i]*int(B[i,k]) for i in range(NC)),Q(0)) for k in range(DEG+1)]
  diag=[Q(0)]*(2*DEG+1);bound=2*max(abs(z) for z in v)
  for ell,A in enumerate(As):
   M=A[1:,1:] if ell==0 else A;M=M[j:j+NC,j:j+NC];bound+=max(Q(int(M[i,i]),den) for i in range(NC));Z=B.T@M@B
   for i,k in itertools.product(range(NC),repeat=2):diag[i+k]+=Q(int(Z[i,k]),den)
  D=add(diag,[2*x for x in coupling]);scalar=[Q(0)]*(2*len(fs)-1)
  for k,a in enumerate(fs):scalar[2*k]=a*alpha**k
  order=55;ex=[Q(0)]*(2*order+1)
  for k in range(order+1):ex[2*k]=alpha**k/math.factorial(k)
  T=alpha*BREAKS[region+1]**2
  if not (0<T<order+2):raise ValueError('The exponential-tail majorant requires a positive geometric denominator.')
  remainder=[Q(0)]*(2*(order+1)+1);remainder[-1]=bound*alpha**(order+1)/(math.factorial(order+1)*(1-T/(order+2)))
  poly=add(add(scalar,conv(D,ex)),remainder);dd=math.lcm(*(v.denominator for v in poly));C=np.array([int(v*dd) for v in poly],dtype=object)
  res=certify(C,dd,[(BREAKS[region],BREAKS[region+1])]);res.update(region=region,exponential_order=order,positive_part_bound=str(bound));out.append(res)
  print('EXACT pair',region,{k:v for k,v in res.items() if k!='failures'},flush=True)
 return out

if __name__=='__main__':
 import sys
 kind=sys.argv[1] if len(sys.argv)>1 else 'full'
 if kind=='pair':pair_checks()
 else:
  for regions in itertools.combinations_with_replacement(range(NS),3 if kind=='full' else 2):check(kind,regions,collect='--collect' in sys.argv)
