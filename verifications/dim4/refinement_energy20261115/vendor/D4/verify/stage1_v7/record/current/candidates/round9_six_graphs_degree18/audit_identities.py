"""Independent rational audits of the piecewise polynomial construction.
These verify identities, not inequality feasibility.
"""
from fractions import Fraction as Q
from pathlib import Path
import itertools,math,json,sys,importlib.util
import numpy as np
from functools import lru_cache

def load(path):
 existing=sys.modules.get('certificate_core')
 if existing is not None and Path(existing.__file__).resolve()==(path/'certificate_core.py').resolve():return existing
 s=importlib.util.spec_from_file_location('core_audit',path/'certificate_core.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def verify(path,only=None):
 core=load(path);As,den,L=core.gram();breaks=tuple(core.BREAKS);ns=len(breaks)-1;deg=core.DEG;nc=deg+1
 def local(r,s):
  spec=json.loads((path/'kernel.json').read_text()).get('radial_basis_intervals')
  left,right=tuple(map(Q,spec[s])) if spec is not None else breaks[s:s+2]
  t=(r-left)/(right-left);return [math.comb(deg,i)*t**i*(1-t)**(deg-i) for i in range(nc)]
 def gegen(t):
  p=[Q(1),t]
  for l in range(1,L):p.append(((2*l+2)*t*p[-1]-l*p[-2])/(l+2))
  return p[:L+1]
 @lru_cache(maxsize=None)
 def K(r,s,t,ir,is_):
  br=local(r,ir);bs=local(s,is_);ps=gegen(t)
  dr=math.lcm(*(x.denominator for x in br));ds=math.lcm(*(x.denominator for x in bs));dp=math.lcm(*(x.denominator for x in ps))
  vr=[int(x*dr) for x in br];vs=[int(x*ds) for x in bs];vp=[int(x*dp) for x in ps]
  answer=0
  for ell,A in enumerate(As):
   shift=1 if ell==0 else 0
   answer+=vp[ell]*sum(int(A[shift+nc*ir+i,shift+nc*is_+j])*vr[i]*vs[j] for i,j in itertools.product(range(nc),repeat=2))
  return Q(answer,den*dr*ds*dp)
 def integer_contract(C,weights):
  ans=C
  for w in reversed(weights):ans=np.tensordot(ans,np.array(w,dtype=object),axes=(-1,0))
  return int(ans)
 def evaluate(C,xyz):
  # Exact homogenization at rational coordinates: one final denominator,
  # rather than millions of Fraction reductions. All intermediate values are integers.
  weights=[];denominator=1
  for axis,x in enumerate(xyz):
   degree=C.shape[axis]-1;p,q=x.numerator,x.denominator
   weights.append([p**i*q**(degree-i) for i in range(degree+1)])
   denominator*=q**degree
  return Q(integer_contract(C,weights),denominator)

 count=0
 for kind in ['full','outside']:
  for regions in itertools.combinations_with_replacement(range(ns),3 if kind=='full' else 2):
   if only is not None and (kind,tuple(regions))!=only:continue
   print('AUDIT branch',kind,regions,flush=True)
   for eta in [Q(0),Q(1,10**8)]:
    C,dd,bounds=core.build(kind,regions,eta)
    for frac in [(Q(1,3),Q(2,5),Q(3,7)),(Q(0),Q(1),Q(1,2))]:
     pt=[a+(b-a)*u for (a,b),u in zip(bounds,frac)];a,b,c=pt
     if kind=='full':
      ia,ib,ic=regions
      value=K(a,b,(a*a+b*b-c*c)/(2*a*b),ia,ib)+K(a,c,(a*a+c*c-b*b)/(2*a*c),ia,ic)+K(b,c,(b*b+c*c-a*a)/(2*b*c),ib,ic)
      expected=(2*a*b*c)**L*(value-eta)
     else:
      ia,ib=regions;t=-1+c*((a+b)**2-breaks[-1]**2)/(2*a*b)
      expected=(2*a*b)**L*(K(a,b,t,ia,ib)-eta)
     got=evaluate(C,pt)/dd
     if got!=expected:raise ValueError('Power polynomial identity failed: '+repr((kind,regions,eta,frac)))
     count+=1
    # Independent Bernstein evaluation equals direct power evaluation.
    B,bd=core.to_bernstein(C.copy(),dd,bounds)
    u=[Q(1,3),Q(2,5),Q(3,7)];basis=[];rden=bd
    for x,n in zip(u,B.shape):
     p,q=x.numerator,x.denominator
     basis.append([math.comb(n-1,j)*p**j*(q-p)**(n-1-j) for j in range(n)])
     rden*=q**(n-1)
    value=Q(integer_contract(B,basis),rden)
    pt=[a+(b-a)*t for (a,b),t in zip(bounds,u)]
    if value!=evaluate(C,pt)/dd:raise ValueError('Bernstein basis transform failed.')
    count+=1
 print('PASS',count,'rational polynomial and Bernstein identities on requested piecewise branches.',flush=True)
 # Explicit Chebyshev recurrence checks and radial endpoint identities.
 P,pden=core.gegen(L)
 for ell in range(L+1):
  for t in [Q(-1),Q(-1,2),Q(0),Q(2,3),Q(1)]:
   val=sum(Q(v,pden)*t**j for j,v in enumerate(P[ell]))
   if val!=gegen(t)[ell]:raise ValueError('Normalized harmonic polynomial mismatch.')
 for region in range(ns):
  B=core.radial(region)
  for t in [Q(0),Q(1,3),Q(1)]:
   r=breaks[region]+t*(breaks[region+1]-breaks[region]);computed=[sum(int(v)*r**j for j,v in enumerate(row)) for row in B]
   if computed!=local(r,region) or sum(computed)!=1:raise ValueError('Radial basis mismatch.')
 print('PASS harmonic recurrence and radial partition identities.',flush=True)
 return {'polynomial_identity_checks':count,'angular_degree':L,'status':'VERIFIED'}
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='--branch':
  verify(Path(__file__).resolve().parent,(sys.argv[2],tuple(map(int,sys.argv[3:]))))
 else:verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
