"""Exact arithmetic crosschecks using a contracted Bernstein transform.

The same four power/direct and two Bernstein-formula/direct evaluations,
harmonic recurrence checks and radial endpoint checks as the pinned audit
are retained. Finite evaluations are not coefficientwise identity proofs.
No dense exact Bernstein tensor is materialized by this alternate audit,
and no implementation of an interval converter is tested or certified.
This module cannot establish any sign condition or a global packing bound.
"""
from fractions import Fraction as Q
from pathlib import Path
import itertools,math,json,sys,importlib.util
import numpy as np
from functools import lru_cache
import hashlib

LEGACY_AUDIT_SHA256='4680c7ffc7d2fefd236785c3f5a42aeb8731dc4340684d254ea29bbc8dd69259'
CONTRACTED_HELPER_SHA256='2d3450fdef53e9dc2570c5e8b38aa02113c154a6de3eae779b22f20956c0ffb1'
_helper_path=Path(__file__).with_name('contracted_bernstein_identity_v2.py')
assert hashlib.sha256(_helper_path.read_bytes()).hexdigest()==CONTRACTED_HELPER_SHA256
_helper_spec=importlib.util.spec_from_file_location('_exact_contracted_bernstein_identity',_helper_path)
_helper=importlib.util.module_from_spec(_helper_spec);_helper_spec.loader.exec_module(_helper)

def _hash_inputs(path):
 names=['certificate_core.py','kernel.json','audit_identities.py']
 if (path/'exact_bernstein_fast.py').is_file():names.append('exact_bernstein_fast.py')
 result={name:hashlib.sha256((path/name).read_bytes()).hexdigest() for name in names}
 result['contracted_bernstein_identity_v2.py']=hashlib.sha256(_helper_path.read_bytes()).hexdigest()
 result['audit_identities_contracted_v2.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 return result

def load(path):
 existing=sys.modules.get('certificate_core')
 if existing is not None and Path(existing.__file__).resolve()==(path/'certificate_core.py').resolve():return existing
 s=importlib.util.spec_from_file_location('core_audit',path/'certificate_core.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def verify(path,only=None):
 path=Path(path).resolve();before=_hash_inputs(path)
 assert before['audit_identities.py']==LEGACY_AUDIT_SHA256
 assert before['contracted_bernstein_identity_v2.py']==CONTRACTED_HELPER_SHA256
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
    # The exact transformed-basis evaluation is contracted before C.
    # This checks the basis-transform formula, not an interval converter.
    u=[Q(1,3),Q(2,5),Q(3,7)]
    value=_helper.contracted_bernstein_evaluation(C,dd,bounds,u)
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
 assert count>0 and (only is None or count==6)
 assert before==_hash_inputs(path),'Exact arithmetic audit input changed'
 return {'polynomial_identity_checks':count,'angular_degree':L,'status':'VERIFIED',
         'verification_scope':'EXACT_ARITHMETIC_CROSSCHECKS_ONLY_NOT_A_GLOBAL_CERTIFICATE',
         'basis_transform_method':'EXACT_TRANSPOSE_AXIS_CONTRACTION',
         'interval_converter_implementation_checked':False,
         'source_sha256':before}
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='--branch':
  verify(Path(__file__).resolve().parent,(sys.argv[2],tuple(map(int,sys.argv[3:]))))
 else:verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
