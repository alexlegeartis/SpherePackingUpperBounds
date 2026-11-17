from pathlib import Path
import json, math, hashlib, random, time
import numpy as np
from flint import fmpz_mat
D=Path(__file__).resolve().parent
S=D.parent/'mw80_normal162_20261117'/'parent'
p=41
# Independent arithmetic: encode a+b*z by a+41*b, z*z=3.
def add(a,b):return ((a%p+b%p)%p)+p*((a//p+b//p)%p)
def neg(a):return (-a%p)+p*((- (a//p))%p)
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 x,y=a%p,a//p;u,v=b%p,b//p
 return (x*u+3*y*v)%p+p*((x*v+y*u)%p)
def inv(a):
 assert a
 x,y=a%p,a//p;n=pow((x*x-3*y*y)%p,-1,p)
 return (x*n)%p+p*((-y*n)%p)
INV=[0]+[inv(a) for a in range(1,p*p)]
assert all(mul(a,INV[a])==1 for a in range(1,p*p))
assert pow(3,20,41)==40

def trim(a):
 while a and a[-1]==0:a.pop()
 return a
def ps(a,b):return trim([sub(a[i] if i<len(a) else 0,b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def pm(a,b):
 r=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):r[i+j]=add(r[i+j],mul(x,y))
 return trim(r)
def pp(a,n):
 r=[1]
 for _ in range(n):r=pm(r,a)
 return r
def rem(a,b):
 a=a[:]
 while len(a)>=len(b) and a:
  shift=len(a)-len(b); c=mul(a[-1],INV[b[-1]])
  for j,v in enumerate(b):a[shift+j]=sub(a[shift+j],mul(c,v))
  trim(a)
 return a
def gcd(a,b):
 while b:a,b=b,rem(a,b)
 return a
z=json.loads((S/'unitary_sections80.json').read_text());assert z['prime']==41 and z['field_nonsquare']==3
sections=[];f=[1]+[0]*41+[1]
for row in z['unitaries']:
 assert len(row)==4 and all(type(x)==int and 0<=x<1681 for x in row)
 a,b,c,d=row
 x=[neg(t)for t in pp([b,a],14)];y=pp([d,c],21)
 assert ps(ps(pm(y,y),pm(pm(x,x),x)),f)==[]
 # No singular point: x and y have no simultaneous zero, at any finite fiber.
 assert len(gcd(x,y))==1
 assert len(x)<=15 and len(y)<=22
 sections.append((x,y))
assert len(sections)==190 and len(set((tuple(x),tuple(y))for x,y in sections))==190
# The 42 finite zeros of f are simple; infinity's scaled equation is Y²=X³+1+s^42.
der=[0]*41+[1]; assert len(gcd(f,der))==1
H=[[0]*190 for _ in range(190)]
finite_hist={};infinity_hist={}
for i,(x,y) in enumerate(sections):
 H[i][i]=14
 for j,(u,v) in enumerate(sections[:i]):
  dx,dy=ps(x,u),ps(y,v)
  if not dx and not ps(y,[neg(t)for t in v]): h=-14
  else:
   assert dx or dy
   finite=len(gcd(dx,dy))-1
   # Valuation at s=0 of s^14 dx(1/s), s^21 dy(1/s).
   orders=([15-len(dx)]if dx else [])+([22-len(dy)]if dy else [])
   infinite=min(orders);assert infinite>=0
   h=7-finite-infinite
   finite_hist[finite]=finite_hist.get(finite,0)+1
   infinity_hist[infinite]=infinity_hist.get(infinite,0)+1
  H[i][j]=H[j][i]=h
 print('row',i,flush=True) if i%40==0 else None
z=np.load(S/'parent80.npz',allow_pickle=False)
assert all(np.issubdtype(z[k].dtype,np.integer)for k in z.files)
G=fmpz_mat(z['G'].tolist());Y=fmpz_mat(z['section_coordinates'].tolist())
assert G.nrows()==80 and G.ncols()==80 and Y.nrows()==190 and Y.ncols()==80
assert G==G.transpose() and Y*G*Y.transpose()==fmpz_mat(H)
assert G.det()==41**12
assert all(fmpz_mat([row[:k] for row in G.tolist()[:k]]).det()>0 for k in range(1,81))
# Saturation without HNF: gcd of selected maximal minors is 1. The lattice
# index divides every maximal minor, so this certifies Y rows generate Z^80.
rng=random.Random(78117);g=0;minors=[]
for attempt in range(100):
 rows=list(range(80))if attempt==0 else sorted(rng.sample(range(190),80))
 det=int(fmpz_mat([Y.tolist()[i]for i in rows]).det());g=math.gcd(g,det)
 minors.append({'rows':rows,'det':str(det),'gcd_so_far':str(g)})
 if g==1:break
assert g==1
out={'status':'PASS','sections':190,'pair_heights':17955,'rank':80,'determinant':str(G.det()),'minimum':14,'finite_intersection_histogram':finite_hist,'infinity_intersection_histogram':infinity_hist,'generation_minor_witnesses':minors,'source_hashes':{f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in [S/'parent80.npz',S/'unitary_sections80.json']},'scope':'Independent polynomial/height/integral-generation replay. Published Shioda height theorem is an explicit mathematical input; full MW-group equality is neither proved nor needed.'}
(D/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='generation_minor_witnesses'},indent=2))
