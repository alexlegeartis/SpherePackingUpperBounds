"""Standalone all-degree angular lower-bound verifier. python-flint 0.9.0 only.
All certificate coordinates, weights, and tail data are rational.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math
from flint import arb,fmpq,ctx
if not __debug__:
 raise RuntimeError("This verifier must not run with Python -O")
ctx.prec=256
q=lambda x:fmpq(F(x).numerator,F(x).denominator)
def verify(path,mutate=None):
 c=json.loads(Path(path).read_text());n=c['dimension'];assert n>=3
 assert c['alpha_pi_multiplier']=='51/4'
 K=c['tail_even_degree']//2;assert 2*K==c['tail_even_degree'] and K>=1
 bound=F(c['absolute_angle_bound']);assert 0<bound<1
 def gegen(t):
  p=[arb(1),t]
  for l in range(1,2*K):p.append(((2*l+n-2)*t*p[-1]-l*p[-2])/(l+n-2))
  return p
 z=[arb(1)]*(2*K+1);W=F(0)
 for j,row in enumerate(c['dual_rows']):
  w=F(row['weight']);assert w>=0
  point=row['point'];t=F(point[1]);T=arb(q(t))
  if point[0]=='cap':
   assert -1<=t<=F(1,2);angles=[T];exps=[arb(1)];W+=w
  else:
   assert point[0]=='cyclic';u=F(point[2]);U=arb(q(u))
   assert 0<=t<=F(1,2) and u>=1 and 2*t*u<=1
   # 1,u,w are the actual lengths divided by the minimum allowed length.
   # Positive triangle area since t in[0,1/2]; w>=u>=1.
   V=(1+U**2-2*U*T).sqrt();assert V>0
   angles=[T,(1-U*T)/V,(U-T)/V]
   alpha=arb(q(c['alpha_pi_multiplier']))*arb.pi()
   exps=[arb(1),(-alpha*arb(q(1-2*t*u))).exp(),(-alpha*arb(q(u*u-2*t*u))).exp()]
   assert 1-2*t*u>=0 and u*u-2*t*u>=0;W+=3*w
  assert all(abs(t)<=arb(q(bound)) for t in angles)
  for t,e in zip(angles,exps):
   p=gegen(t)
   for l in range(2*K+1):z[l]+=arb(q(w))*e*p[l]
 assert all(x>0 for x in z[1:2*K])
 # Integral representation: modulus bound for every harmonic l>=2K.
 moment=F(0);m=F(1)
 for j in range(K+1):
  if j:m*=F(2*j-1,n+2*j-3)
  moment+=math.comb(K,j)*bound**(2*(K-j))*(1-bound**2)**j*m
 assert W*moment<1
 low=F(c['lower_bound_B']);assert z[0]>arb(q(low))
 out={'pass':True,'dimension':n,'lower_bound_B_exact':str(low),'dual_value_interval':str(z[0]),'finite_degrees_checked':[1,2*K-1],'all_degrees_from':2*K,'tail_error_bound_exact':str(W*moment),'tail_error_bound_float':float(W*moment),'all_harmonic_degrees_certified':True,'rows':len(c['dual_rows'])}
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('certificate');ap.add_argument('--out');a=ap.parse_args();r=verify(a.certificate)
 print(json.dumps(r,indent=2))
 if a.out:Path(a.out).write_text(json.dumps(r,indent=2)+'\n')
