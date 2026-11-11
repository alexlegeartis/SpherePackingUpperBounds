"""Rigorous bounded interval replay of a six-edge tetrahedron candidate.
Not a density verifier. A timeout or remaining box is an explicit failure.
"""
import argparse,itertools,json,time
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as Poly,arb,ctx
from verify_tetra_envelopes import gegenbauer,negative_on
ctx.prec=96
class Jet:
 def __init__(self,v,g=None):self.v=arb(v);self.g=g if g is not None else [arb(0)for _ in range(6)]
 def __add__(self,b):
  if not isinstance(b,Jet):b=Jet(b)
  return Jet(self.v+b.v,[a+c for a,c in zip(self.g,b.g)])
 __radd__=__add__
 def __neg__(self):return Jet(-self.v,[-a for a in self.g])
 def __sub__(self,b):return self+-b
 def __rsub__(self,b):return -self+b
 def __mul__(self,b):
  if not isinstance(b,Jet):b=Jet(b)
  return Jet(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.g,b.g)])
 __rmul__=__mul__
 def inv(self):return Jet(1/self.v,[-a/(self.v*self.v)for a in self.g])
 def __truediv__(self,b):return self*(b.inv()if isinstance(b,Jet)else 1/arb(b))
 def __rtruediv__(self,b):return self.inv()*b
 def __pow__(self,k):
  assert k>=0 and isinstance(k,int)
  if not k:return Jet(1)
  result=Jet(1)
  for _ in range(k):result=result*self
  return result
 def sqrt(self):
  s=self.v.sqrt();return Jet(s,[a/(2*s)for a in self.g])

def evaluate(box,pol,jet):
 X=[]
 for i,(l,r)in enumerate(box):
  v=arb((l+r)/2, (r-l)/2)
  if jet:
   grad=[arb(0)for _ in range(6)];grad[i]=arb(1);v=Jet(v,grad)
  X.append(v)
 edges=list(itertools.combinations(range(4),2));index={e:i for i,e in enumerate(edges)}
 def edge(i,j):return X[index[tuple(sorted((i,j)))]]
 total=0
 for i,j in edges:
  k,l=[z for z in range(4)if z not in(i,j)]
  a=edge(i,j);u=edge(i,k);v=edge(j,k);w=edge(i,l);z=edge(j,l);c=edge(k,l)
  du=u-v;dw=w-z
  R=(u+v)/2-a/4-du*du/(4*a);S=(w+z)/2-a/4-dw*dw/(4*a)
  dd=u-v-w+z;delta=dd*dd/(4*a)
  t=(R+S+delta-c)/(2*(R*S).sqrt())
  value=0
  for co in reversed(pol):value=value*t+arb(co)
  total+=value
 return total/6

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,default=22);ap.add_argument('--q',default='21/20');ap.add_argument('--seconds',type=float,default=90);ap.add_argument('--output',required=True);args=ap.parse_args()
 root=Path(__file__).resolve().parent;raw=json.loads((root/'tetra_full_geometry_screen.json').read_text());q=Q(args.q)
 row=next(r for r in raw if r['dimension']==args.dimension and Q(str(F(str(r['q']))))==q)
 cs=[Q(str(F(str(c))))for c in row['coefficients']];assert cs[0]==1 and all(c>=0 for c in cs)
 basis=gegenbauer(args.dimension-1,len(cs)-1);P=sum((a*p for a,p in zip(cs,basis)),Poly([]));far=1-q/(2*(q-Q(1,4)))
 negative_on(P,Q(-1),far);print('Exact far sign PASS; B',float(sum(cs)),flush=True)
 pol=P.coeffs();stack=[([(Q(1),q)]*6,0)];leaves=0;processed=0;depth=0;start=time.time();failed_positive=False;face_reductions=0
 while stack and time.time()-start<args.seconds:
  box,d=stack.pop();processed+=1;depth=max(depth,d)
  mid=[((l+r)/2,(l+r)/2)for l,r in box]
  J=evaluate(box,pol,True);center=evaluate(mid,pol,False)
  assert center.is_finite() and J.v.is_finite() and all(g.is_finite()for g in J.g), 'Nonfinite interval'
  assert d<100, 'Subdivision limit'
  upper=center+sum((arb(0,(r-l)/2)*g for(l,r),g in zip(box,J.g)),arb(0))
  if upper<0 or J.v<0:leaves+=1;continue
  if center>0:print('Positive center',box,center,flush=True);failed_positive=True;stack.append((box,d));break
  reduced=box.copy();changed=False
  for i,((l,r),g)in enumerate(zip(box,J.g)):
   if l<r and (g>0 or g<0):
    endpoint=r if g>0 else l;reduced[i]=(endpoint,endpoint);changed=True
  if changed:
   face_reductions+=1;stack.append((reduced,d+1));continue
  scores=[float(abs(g).upper())*float(r-l)for(l,r),g in zip(box,J.g)]
  ix=max(range(6),key=lambda i:scores[i]);l,r=box[ix];m=(l+r)/2
  a=box.copy();b=box.copy();a[ix]=(l,m);b[ix]=(m,r);stack.extend([(a,d+1),(b,d+1)])
  if processed%1000==0:print('processed',processed,'leaves',leaves,'pending',len(stack),'depth',depth,'seconds',time.time()-start,flush=True)
 out={'dimension':args.dimension,'q':str(q),'coefficients':list(map(str,cs)),'B_exact':str(sum(cs)),'polynomial':list(map(str,pol)),'exact_far_sign':True,'full_tetra_interval_pass':not stack,'processed':processed,'verified_leaves':leaves,'monotonic_face_reductions':face_reductions,'pending_boxes':len(stack),'maximum_depth':depth,'seconds':time.time()-start,'warning':'Auxiliary aggregate constant only; no density claim.'}
 Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
 assert not stack,'Full tetrahedron verification did not finish'
if __name__=='__main__':main()
