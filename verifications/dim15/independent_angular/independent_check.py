"""Independent angular reconstruction, Bernstein cap, root isolation, AD cyclic proof.
Reads only frozen exact angular data. No radial candidate or density bound.
"""
if not __debug__:
    raise RuntimeError('Assertions are required; do not run this verifier with -O or -OO.')

from pathlib import Path
from fractions import Fraction as F
import json,math,time,hashlib
import sympy as sp
from flint import fmpq,fmpq_poly,arb,arb_poly,ctx
ctx.prec=320
D=Path(__file__).resolve().parent;path=D/'inputs/angular_exact.json';assert hashlib.sha256(path.read_bytes()).hexdigest()=='ebeb0fd1b0532445cf274f01f6f01f1996f608fc9794eb09e3cb43f3932f6354'
j=json.loads(path.read_text());n=j['dimension'];assert n==15 and F(j['a_exact'])==12 and F(j['b_exact'])==16
q=lambda v:fmpq(F(v).numerator,F(v).denominator)
rise=lambda x,k:math.prod((x+i for i in range(k)),start=F(1))
c=list(map(F,j['gegenbauer_coefficients']));p=list(map(F,j['polynomial_coefficients']));assert c[0]==1 and min(c)>=0
# Closed coefficient formula, independent of the source's three-term recurrence.
got=[F(0)]*len(c)
for ell,v in enumerate(c):
 norm=rise(F(n-2),ell)/math.factorial(ell)
 for k in range(ell//2+1):got[ell-2*k]+=v*(-1)**k*rise(F(n-2,2),ell-k)*2**(ell-2*k)/math.factorial(k)/math.factorial(ell-2*k)/norm
while got and got[-1]==0:got.pop()
assert got==p and sum(c)==F(j['B_exact']);Q=fmpq_poly(list(map(q,p)));B=q(sum(c));margin=fmpq(1,20000);cap=Q+margin
# Exact Bernstein coefficients certify the cap without the source's Sturm code.
cap_stack=[(fmpq(-1),fmpq(1,2),0)];cap_leaves=[];deg=cap.degree()
while cap_stack:
 l,r,depth=cap_stack.pop();poly=cap(fmpq_poly([l,r-l]));bb=[sum((poly[k]*fmpq(math.comb(i,k),math.comb(deg,k))for k in range(i+1)),fmpq(0))for i in range(deg+1)]
 if max(bb)<0:cap_leaves.append([str(l),str(r),str(max(bb))]);continue
 assert depth<32;mid=(l+r)/2;cap_stack.extend([(l,mid,depth+1),(mid,r,depth+1)])
roots=Q.numer().complex_roots();assert sum(m for z,m in roots)==Q.degree();possible=[]
for z,m in roots:
 if z.imag.contains(0) and not(z.real<arb(fmpq(1,2))or z.real>1):
  assert z.real>arb(fmpq(1,2))and z.real<1;possible.append((z,m))
assert sum(m for z,m in possible)==1 and Q(fmpq(1,2))<0 and Q(1)>0
# Independent symbolic audit of all manual angle/exponent partials.
t,u=sp.symbols('t u',real=True);dd=1-2*t*u;ee=u*u-2*t*u;ww=sp.sqrt(1+ee);t1=(1-u*t)/ww;t2=(u-t)/ww
w_t=-u/ww;w_u=(u-t)/ww;E1=-12*dd-16*dd*(2*u*u-2+dd);E2=-12*ee-16*ee*ee
manual=[(-u/ww-t1*w_t/ww,sp.diff(t1,t)),(-t/ww-t1*w_u/ww,sp.diff(t1,u)),(-1/ww-t2*w_t/ww,sp.diff(t2,t)),(1/ww-t2*w_u/ww,sp.diff(t2,u)),(-(12+32*ee)*(-2*u),sp.diff(E1,t)),(-(12+32*ee)*(-2*t)-64*dd*u,sp.diff(E1,u)),(-(12+32*ee)*(-2*u),sp.diff(E2,t)),(-(12+32*ee)*(2*(u-t)),sp.diff(E2,u))]
assert all(sp.simplify(a-b)==0 for a,b in manual)
assert sp.expand(1+u*(u-1)+u*(1-2*t)-(1+ee))==0
print('EXACT CHECKS PASS',len(cap_leaves),'Bernstein cap leaves; one transition root; eight symbolic derivatives',flush=True)
AQ=arb_poly([arb(q(v))for v in p]);DQ=AQ.derivative()
# Centered Taylor range enclosure. Clipping is valid because these quantities
# are cosines of a real triangle for every t in[-1,.5],u>=1, not only sorted ones.
def prange(poly,x):
 lo=max(arb(-1),x.lower());hi=min(arb(1),x.upper());assert lo<=hi
 mid=(lo+hi)/2;radius=(hi-lo)/2;cc=poly(arb_poly([mid,radius]));err=sum((abs(v)for v in cc.coeffs()[1:]),arb(0));return cc[0]+arb(0,err.upper())
class AD:
 def __init__(self,v,dt=0,du=0):self.v=arb(v);self.dt=arb(dt);self.du=arb(du)
 @staticmethod
 def make(x):return x if isinstance(x,AD)else AD(x)
 def __add__(self,o):o=AD.make(o);return AD(self.v+o.v,self.dt+o.dt,self.du+o.du)
 __radd__=__add__
 def __neg__(self):return AD(-self.v,-self.dt,-self.du)
 def __sub__(self,o):return self+-AD.make(o)
 def __rsub__(self,o):return AD.make(o)+-self
 def __mul__(self,o):o=AD.make(o);return AD(self.v*o.v,self.dt*o.v+self.v*o.dt,self.du*o.v+self.v*o.du)
 __rmul__=__mul__
 def inv(self):
  iv=1/self.v;return AD(iv,-self.dt*iv*iv,-self.du*iv*iv)
 def __truediv__(self,o):return self*AD.make(o).inv()
 def sqrt(self):v=self.v.sqrt();return AD(v,self.dt/(2*v),self.du/(2*v))
 def exp(self):v=self.v.exp();return AD(v,v*self.dt,v*self.du)
 def angular(self):d=prange(DQ,self.v);return AD(prange(AQ,self.v),d*self.dt,d*self.du)
def evaluate(T,U):
 T=AD(T,1,0);U=AD(U,0,1)
 # This algebraic form makes W²>=1 transparent throughout the entire box.
 Wsq=1+U*(U-1)+U*(1-2*T)
 # Arb ball multiplication can lose one-sided positivity. Intersect with the
 # proven analytic lower bound1 before taking the square root.
 if not Wsq.v>=1:Wsq.v=arb(1).union(Wsq.v.upper())
 assert Wsq.v>0
 W=Wsq.sqrt();e=Wsq-1;d=1-2*T*U
 c1=(1-U*T)/W;c2=(U-T)/W
 return T.angular()+(-d*(12+16*(2*U*U-2+d))).exp()*c1.angular()+(-e*(12+16*e)).exp()*c2.angular()
start=time.time();stack=[(fmpq(-1),fmpq(1,2),fmpq(1),fmpq(4),0)];accept=[];skip=[];tested=0;maxdepth=0;area=fmpq(0)
while stack:
 l,r,v,w,depth=stack.pop();tested+=1
 assert time.time()-start<120 and tested<250000,('BUDGET',tested,len(stack))
 if 2*l*v>1:skip.append([str(x)for x in [l,r,v,w]]);area+=(r-l)*(w-v);continue
 f=evaluate(arb(l).union(arb(r)),arb(v).union(arb(w)));assert all(z.is_finite()for z in [f.v,f.dt,f.du]),('NONFINITE',l,r,v,w,str(f.v),str(f.dt),str(f.du));method='direct';upper=f.v
 if not upper<0:
  center=evaluate(arb((l+r)/2),arb((v+w)/2)).v;upper=center+abs(f.dt)*arb((r-l)/2)+abs(f.du)*arb((w-v)/2);method='AD_mean_value'
 if upper<0:accept.append({'box':[str(x)for x in [l,r,v,w]],'method':method,'upper_enclosure':str(upper)});area+=(r-l)*(w-v);maxdepth=max(maxdepth,depth);continue
 assert depth<48,('DEPTH',l,r,v,w)
 # Rational width-based splitting uses a different aspect ratio from the primary.
 if 3*(r-l)>=2*(w-v):
  m=(l+r)/2;stack.extend([(l,m,v,w,depth+1),(m,r,v,w,depth+1)])
 else:
  m=(v+w)/2;stack.extend([(l,r,v,m,depth+1),(l,r,m,w,depth+1)])
assert area==fmpq(9,2)
small=(-arb(3780)).exp();tail1=arb(B)*small;tail2=arb(B)*((-arb(36)).exp()+small);assert tail1<arb(margin)and tail2<arb(margin)
# For x<=y<=z, all >=s, ratio h(lambda z)/h(lambda y) has
# log derivative -(z-y)/s*(12-32+32*lambda*(z+y)/s)<=0.
# Minimum bracket is12+32=44, so both required ratios decrease.
out={'pass':True,'scope':'Angular PSD, exact cap/sign transition, full cyclic geometry and all-scale theorem only. No radial or density claim.','angular_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'precision_bits':ctx.prec,'closed_Gegenbauer_reconstruction_pass':True,'exact_Bernstein_cap_leaves':len(cap_leaves),'one_sign_transition_root':str(possible[0][0]),'all_eight_symbolic_partials_match':True,'denominator_lower_bound_squared_exact':'1','independent_AD_interval_cyclic_pass':True,'accepted_rectangles':len(accept),'discarded_sorting_rectangles':len(skip),'tested':tested,'max_binary_depth':maxdepth,'exact_partition_area':str(area),'tail_case_t_ge_minus_one_quarter_upper_positive_part':str(tail1),'tail_case_t_le_minus_one_quarter_upper_positive_part':str(tail2),'all_scale_ratio_log_derivative_bracket_lower_bound':'44','all_scale_extension_pass':True,'seconds':time.time()-start}
(D/'independent_result.json').write_text(json.dumps(out,indent=2)+'\n');(D/'interval_partition.json').write_text(json.dumps({'accepted':accept,'discarded':skip,'exact_covered_area':str(area)},indent=2)+'\n');(D/'cap_partition.json').write_text(json.dumps(cap_leaves,indent=2)+'\n');print(out,flush=True)
