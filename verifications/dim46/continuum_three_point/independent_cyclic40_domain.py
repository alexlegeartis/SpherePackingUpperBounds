import sys,json,math
sys.path.insert(0,'solver_packages')
from flint import fmpq_poly,fmpq,arb,arb_poly,ctx
ctx.prec=192
from pathlib import Path
p=Path('continuum_three_point');d=json.loads((p/'cyclic_angular40_exact.json').read_text());Q=fmpq_poly(list(map(fmpq,d['polynomial_coefficients'])));AQ=arb_poly([arb(x) for x in Q.coeffs()]);B=arb(fmpq(d['B_exact']));c=fmpq(5001,10000);cache={}
def tupper(a,b):
 if (a,b) not in cache:
  z=Q(fmpq_poly([a,b-a]));bs=[sum((z[j]*fmpq(math.comb(i,j),math.comb(19,j)) for j in range(i+1)),fmpq(0)) for i in range(20)];cache[a,b]=max(bs)
 return cache[a,b]
def plus(x):
 lo=max(arb(-1),x.lower());hi=min(arb(1),x.upper());mid=(lo+hi)/2;rad=(hi-lo)/2
 coeff=AQ(arb_poly([mid,rad]));err=sum((abs(coeff[j]) for j in range(1,len(coeff))),arb(0))
 return (coeff[0]-err).union(coeff[0]+err)
stack=[(fmpq(0),fmpq(1,2),fmpq(1),fmpq(5),0,0)];passed=skipped=0;depth=0
while stack:
 a,b,s,t,l,m=stack.pop()
 if 2*a*s>1:skipped+=1;continue
 U=tupper(a,b);T=arb(a).union(arb(b));V=arb(s).union(arb(t));W=(1+V*V-2*T*V).sqrt()
 okay=False
 if W.is_finite():
  first=(1-V*T)/W;second=(V-T)/W
  g=1-2*T*V;h=V*V-2*T*V
  val=arb(U)+plus(first)*(-40*g).exp()+plus(second)*(-40*h).exp();okay=val<0
 if okay:passed+=1;depth=max(depth,l,m);continue
 assert max(l,m)<25
 if U>=0 or l<=m:
  mid=(a+b)/2;stack.extend([(a,mid,s,t,l+1,m),(mid,b,s,t,l+1,m)])
 else:
  mid=(s+t)/2;stack.extend([(a,b,s,mid,l,m+1),(a,b,mid,t,l,m+1)])
assert fmpq(40**30,math.factorial(30))>2*fmpq(d['B_exact'])*10**5
out={'pass':True,'precision':192,'certified_rectangles':passed,'discarded_rectangles':skipped,'max_depth':depth,'tails_exponential_check':True}
(p/'independent_cyclic40_domain.json').write_text(json.dumps(out,indent=2));print(out)
