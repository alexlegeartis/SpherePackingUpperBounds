"""Independent analytic-gradient and rational-cover check of fixed angular geometry."""
if not __debug__:raise RuntimeError('Assertions required')
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json,time,hashlib
from flint import arb,arb_poly,fmpq,fmpq_poly,ctx
ctx.prec=256;D=Path(__file__).resolve().parent;start=time.monotonic()
fit=json.loads((D/'fit.json').read_text(),parse_float=str);source=json.loads((D/'angular.json').read_text());part=json.loads((D/'partition.json').read_text());claimed=json.loads((D/'RESULT_GEOMETRY.json').read_text())
def qa(x):
 q=F(x);return arb(fmpq(q.numerator,q.denominator))
def qp(p):return fmpq_poly([fmpq(x.numerator,x.denominator)for x in p])
p=list(map(F,fit['parameters']));assert len(p)==14 and min(p)>=0
VK=list(map(F,['1.02','1.05','1.1','1.2','1.4','2']));LK=list(map(F,['1.2','1.4','1.5','1.6','1.65','1.75']))
# Knots are fixed by the construction, identically to primary verifier.
assert source['dimension']==10
c=list(map(F,source['gegenbauer_coefficients']));assert c[0]==1 and min(c)>=0;B=sum(c);assert B==F(source['B_exact'])
def rising(a,m):
 out=F(1)
 for j in range(m):out*=a+j
 return out
coeff=[F(0)]*len(c)
for ell,cl in enumerate(c):
 norm=rising(F(8),ell)/factorial(ell)
 for j in range(ell//2+1):coeff[ell-2*j]+=cl*((-1)**j)*rising(F(4),ell-j)*2**(ell-2*j)/(factorial(j)*factorial(ell-2*j)*norm)
while len(coeff)>1 and coeff[-1]==0:coeff.pop()
assert coeff==list(map(F,source['polynomial_coefficients']))
Q=qp(coeff);AQ=arb_poly([qa(x)for x in coeff]);DQ=AQ.derivative()
def sturm(P):
 out=[P,P.derivative()]
 while out[-1].degree()>0:
  rem=-(out[-2]%out[-1]);assert rem;out.append(rem)
 return out

def count(P,a,b):
 seq=sturm(P)
 def variation(x):
  signs=[1 if q(x)>0 else -1 for q in seq if q(x)!=0];return sum(signs[i]!=signs[i+1]for i in range(len(signs)-1))
 return variation(fmpq(a.numerator,a.denominator))-variation(fmpq(b.numerator,b.denominator))
assert Q(-1)+fmpq(1,20000)<0 and Q(fmpq(1,2))+fmpq(1,20000)<0
assert count(Q+fmpq(1,20000),F(-1),F(1,2))==0
assert Q(fmpq(1,2))<0 and Q(1)>0 and count(Q,F(1,2),F(1))==1
gamma=sum(p[:7]);assert gamma>=7
assert qa(B)*((-3*qa(gamma)).exp()+(-14*qa(gamma)).exp())<qa(F(1,20000))
# End-point range evaluation uses monotonicity, not AD positive-part code.
def maxzero(x):
 if x<=0:return arb(0)
 if x>=0:return x
 return arb(0).union(x.upper())

def point_phi(v):
 ans=qa(p[0])*(v-1)+qa(p[7])*v.log()
 for b,k in zip(p[1:7],VK):ans+=qa(b)*maxzero(v-qa(k))
 for b,k in zip(p[8:],LK):ans+=qa(b)*maxzero((v/qa(k)).log())
 return ans

def phi(v):return point_phi(v.lower()).union(point_phi(v.upper()))
def slope(v):
 ans=qa(p[0])+qa(p[7])/v
 for b,k in zip(p[1:7],VK):
  if v.upper()<=qa(k):continue
  term=qa(b)
  if v.lower()<qa(k):term=arb(0).union(term)
  ans+=term
 for b,k in zip(p[8:],LK):
  if v.upper()<=qa(k):continue
  term=qa(b)/v
  if v.lower()<qa(k):term=arb(0).union(term)
  ans+=term
 return ans

def poly_range(poly,x):
 lo=max(arb(-1),x.lower());hi=min(arb(1),x.upper());assert lo<=hi;mid=(lo+hi)/2;rad=(hi-lo)/2;expanded=poly(arb_poly([mid,rad]));error=sum((abs(z)for z in expanded.coeffs()[1:]),arb(0));return expanded[0]+arb(0,error.upper())

def evaluate(t,u,derivative=True):
 y=u*u;z=1+u*(u-1)+u*(1-2*t)
 if not y>=1:y=arb(1).union(y.upper())
 if not z>=1:z=arb(1).union(z.upper())
 w=z.sqrt();inv=1/w;inv3=inv*inv*inv;one=1-u*t;gap=u-t;v1=one*inv;v2=gap*inv
 r1=(phi(y)-phi(z)).exp();r2=(-phi(z)).exp();q0=poly_range(AQ,t);q1=poly_range(AQ,v1);q2=poly_range(AQ,v2);value=q0+r1*q1+r2*q2
 if not derivative:return value
 sy=slope(y);sz=slope(z);q1d=poly_range(DQ,v1);q2d=poly_range(DQ,v2)
 c1t=-u*inv+one*u*inv3;c1u=-t*inv-one*gap*inv3;c2t=-inv+gap*u*inv3;c2u=inv-gap*gap*inv3
 d1t=2*u*sz;d1u=2*u*sy-2*gap*sz;d2t=2*u*sz;d2u=-2*gap*sz
 dt=poly_range(DQ,t)+r1*(d1t*q1+q1d*c1t)+r2*(d2t*q2+q2d*c2t)
 du=r1*(d1u*q1+q1d*c1u)+r2*(d2u*q2+q2d*c2u)
 assert value.is_finite()and dt.is_finite()and du.is_finite();return value,dt,du

accepted={tuple(map(F,x['box']))for x in part['accepted']};discarded={tuple(map(F,x))for x in part['discarded']};assert len(accepted)==len(part['accepted'])and len(discarded)==len(part['discarded']);assert not accepted&discarded
# Reconstruct the full exact binary tree; area alone would be insufficient.
stack=[(F(-1),F(1,2),F(1),F(4),0)];seen=set();tested=0
while stack:
 l,r,v,w,depth=stack.pop();tested+=1;box=(l,r,v,w)
 if box in accepted or box in discarded:
  seen.add(box)
  if box in discarded:assert 2*l*v>1
  continue
 assert depth<48
 if 3*(r-l)>=2*(w-v):
  m=(l+r)/2;stack.extend([(l,m,v,w,depth+1),(m,r,v,w,depth+1)])
 else:
  m=(v+w)/2;stack.extend([(l,r,v,m,depth+1),(l,r,m,w,depth+1)])
assert seen==accepted|discarded and tested==claimed['tested']
assert sum(((r-l)*(w-v)for l,r,v,w in seen),F(0))==F(9,2)
print('EXACT STRUCTURE AND COVER PASS',len(accepted),len(discarded),flush=True)
# Check every accepted rectangle with a separately derived analytic gradient.
pending=[(*box,0)for box in accepted];checks=0;refinements=0;maximum_refinement=0
while pending:
 l,r,v,w,depth=pending.pop();checks+=1;assert time.monotonic()-start<600
 value,dt,du=evaluate(qa(l).union(qa(r)),qa(v).union(qa(w)))
 upper=value
 if not upper<0:
  center=evaluate(qa((l+r)/2),qa((v+w)/2),False);upper=center+abs(dt)*qa((r-l)/2)+abs(du)*qa((w-v)/2)
 if upper<0:maximum_refinement=max(maximum_refinement,depth);continue
 assert depth<8,('INDEPENDENT_ENCLOSURE_UNRESOLVED',str(l),str(r),str(v),str(w))
 refinements+=1
 if 3*(r-l)>=2*(w-v):
  m=(l+r)/2;pending.extend([(l,m,v,w,depth+1),(m,r,v,w,depth+1)])
 else:
  m=(v+w)/2;pending.extend([(l,r,v,m,depth+1),(l,r,m,w,depth+1)])
out={'status':'PASS','scope':'Independent frozen log-hinge angular geometry and theorem checks only; no radial or density bound asserted','exact_Gegenbauer_closed_form':True,'exact_cap_and_one_transition':True,'all_scale_hypotheses':True,'rigorous_tail':True,'exact_binary_cover':{'accepted':len(accepted),'discarded':len(discarded),'tested_tree_nodes':tested,'area':'9/2'},'independent_analytic_gradient_checks':checks,'additional_refinements':refinements,'maximum_extra_depth':maximum_refinement,'seconds':time.monotonic()-start,'input_hashes':{fn:hashlib.sha256((D/fn).read_bytes()).hexdigest()for fn in ['fit.json','angular.json','partition.json']}}
(D/'VERIFIED.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
