"""Independent exact angular signs and global triangle-domain enclosure."""
import json,math
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as Poly,arb,arb_poly,ctx
assert __debug__
ctx.prec=256;D=Path(__file__).resolve().parent;d=json.loads((D/'angular.json').read_text());assert d['dimension']==38 and d['degree']==17
p=Poly(list(map(Q,d['polynomial_coefficients'])));B=Q(d['B_exact']);c=Q(d['cutoff']);assert c==Q(5001,10000) and p(1)==B

def chain(p):
 s=[p,p.derivative()]
 while s[-1].degree()>0:
  r=-(s[-2]%s[-1]);assert r;z=r.numer();s.append(Poly(z/abs(z.content())))
 return s
def variation(s,t):
 v=[p(t)for p in s];sg=[int(x>0)-int(x<0)for x in v if x];return sum(a!=b for a,b in zip(sg,sg[1:]))
s=chain(p);assert p(-1)<0 and p(c)<0 and p(1)>0
assert variation(s,Q(-1))==variation(s,c) and variation(s,c)-variation(s,Q(1))==1
cap=p+Q(1,100000);s=chain(cap);assert cap(-1)<0 and cap(Q(1,2))<0 and variation(s,Q(-1))==variation(s,Q(1,2))
# PSD coefficient identity is separately verified by the radial/algebra verifier.
# Its positive harmonic expansion implies |Q(t)|<=B on [-1,1].
derivatives=[];z=p
for j in range(20):derivatives.append(arb_poly([arb(x) for x in (z/math.factorial(j)).coeffs()]));z=z.derivative()
cache={}
def exact_bernstein(lo,hi):
 key=(lo,hi)
 if key not in cache:
  powers=[sum((p[k]*Q(math.comb(k,j))*lo**(k-j)*(hi-lo)**j for k in range(j,20)),Q(0))for j in range(20)]
  cache[key]=max(sum((powers[j]*Q(math.comb(i,j),math.comb(19,j))for j in range(i+1)),Q(0))for i in range(20))
 return cache[key]
stack=[(Q(0),Q(1,2),Q(1),Q(5),0,0)];passed=excluded=maxdepth=0
while stack:
 tl,th,ul,uh,dt,du=stack.pop()
 if 2*tl*ul>1:excluded+=1;continue
 T=arb(tl).union(arb(th));U=arb(ul).union(arb(uh));w=(1+U*U-2*T*U).sqrt();bound=exact_bernstein(tl,th)
 okay=False
 if w.is_finite():
  g1=1-2*T*U;g2=U*U-2*T*U
  # Exponent intervals multiplied by upper Q are safe only when Q's sign
  # is known; use a full Q enclosure rather than a scalar upper here.
  def polynomial_interval(interval):
   lo=max(arb(-1),interval.lower());hi=min(arb(1),interval.upper());m=(lo+hi)/2;r=(hi-lo)/2
   center=derivatives[0](m);e=sum((abs(derivatives[j](m))*r**j for j in range(1,20)),arb(0))
   return (center-e).union(center+e)
  total=arb(bound)+polynomial_interval((1-T*U)/w)*(-40*g1).exp()+polynomial_interval((U-T)/w)*(-40*g2).exp()
  okay=total<0
 if okay:passed+=1;maxdepth=max(maxdepth,dt,du);continue
 assert max(dt,du)<26
 if bound>=0 or dt<=du:
  tm=(tl+th)/2;stack.extend([(tl,tm,ul,uh,dt+1,du),(tm,th,ul,uh,dt+1,du)])
 else:
  um=(ul+uh)/2;stack.extend([(tl,th,ul,um,dt,du+1),(tl,th,um,uh,dt,du+1)])
assert Q(40**30,math.factorial(30))>2*B*100000
out={'pass':True,'dimension':38,'angular_Sturm_signs':True,'compact_triangles_passed':passed,'exact_sorting_exclusions':excluded,'max_depth':maxdepth,'analytic_tail_exponential_bound':True,'scale_lemma_requires_single_root_verified':True}
(D/'angular_cyclic_audit.json').write_text(json.dumps(out,indent=2));print(out)
