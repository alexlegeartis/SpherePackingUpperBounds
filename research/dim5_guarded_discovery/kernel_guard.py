"""Exact spherical-polynomial kernel guard; search and finite export helpers.

F(t)<=-1/2 on [-1,7/9], F has nonnegative Gegenbauer coefficients.
The annulus-constant PSD kernel has pair correction identically zero and
strictly negative surviving geometric terms, at a known isolated-origin cost.
"""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
import argparse,json,math
from pathlib import Path
from fractions import Fraction as Fr
import numpy as np
from scipy.special import eval_gegenbauer,poch,comb,factorial
from scipy.optimize import linprog
from flint import fmpq as Q,fmpq_poly as P

def make_guard(n,degree=12):
 s=Q(7,9);lam=(n-2)/2;mats=[]
 grid=np.linspace(-1,float(s),81)
 for left,right in zip(grid[:-1],grid[1:]):
  jets=np.zeros((degree+1,degree+1));bc=np.zeros_like(jets)
  for k in range(degree+1):
   for ell in range(k,degree+1):jets[k,ell]=2**k*poch(lam,k)*eval_gegenbauer(ell-k,lam+k,left)/eval_gegenbauer(ell,lam,1)/factorial(k)
  for i in range(degree+1):
   for k in range(i+1):bc[i,k]=comb(i,k)/comb(degree,k)*(right-left)**k
  mats.append(bc@jets)
 A=np.vstack(mats);sol=linprog(np.ones(degree),A_ub=A[:,1:],b_ub=-A[:,0]-1e-5,bounds=[(0,None)]*degree,method='highs');assert sol.success,sol.message
 cs=[Q(1)]+[Q(str(Fr(str(v))))for v in sol.x];assert min(cs)>=0
 X=P([0,1]);G=[P([1]),X]
 for k in range(1,degree):G.append(((2*k+n-2)*X*G[-1]-k*G[-2])/Q(k+n-2))
 pol=sum((q*g for q,g in zip(cs,G)),P([]));st=[pol,pol.derivative()]
 while st[-1].degree()>0:
  rr=-(st[-2]%st[-1])
  if not rr:break
  z=rr.numer();st.append(P(z/abs(z.content())))
 def var(t):
  sg=[1 if q(t)>0 else -1 for q in st if q(t)];return sum(a!=b for a,b in zip(sg,sg[1:]))
 assert pol(-1)<0 and pol(s)<0 and var(Q(-1))==var(s)
 cs[0]-=Q(1,2);total=sum(cs,Q(0));origin=total*total/(4*cs[0]);cross=-total/2
 assert cs[0]>0 and min(cs)>=0 and origin*cs[0]==cross*cross
 return {'dimension':n,'angular_degree':degree,'max_angle_inner_product':'7/9','strict_sign_upper':'-1/2','gegenbauer_coefficients':list(map(str,cs)),'kernel_origin':str(origin),'kernel_cross':str(cross),'kernel_diagonal':str(total),'pair_correction':'0','exact_sturm_pass':True}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--degree',type=int,default=12);args=ap.parse_args();r=make_guard(args.dimension,args.degree);args.output.write_text(json.dumps(r,indent=2)+'\n');print('EXACT_KERNEL_GUARD',args.dimension,float(Q(r['kernel_diagonal'])),float(Q(r['kernel_origin'])))
