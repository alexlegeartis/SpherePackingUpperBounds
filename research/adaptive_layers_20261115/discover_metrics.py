"""Optional finite numerical metric discovery; not a proof auditor."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
import numpy as np,json,time
from scipy.optimize import linprog
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
P=args.output;P.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(202611155285)
R5=np.array([[int(j==i)-int(j==i+1)for j in range(5)]for i in range(4)]+[[int(j in(3,4))for j in range(5)]],dtype=np.int64)
Q5=R5@R5.T;ips=R5@R5[0];m85=1.25*Q5+.25*np.outer(ips,ips)
cases={
 'p48_r3':(np.array([[2,-1,0],[-1,2,-1],[0,-1,2.5]]),lambda v:6 if not any(v%2)else 2.5 if v[2]%2 else 2),
 'p48_r4':(np.array([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2.5,-.25],[0,-1,-.25,2.5]]),lambda v:6 if not any(v%2)else 2.5 if any(v[2:]%2)else 2),
 'adaptive85':(m85,lambda v:8 if not any(v%2)else 3.5 if list(v%2)==[1,0,0,0,0]else 2.5 if int(v@Q5@v)//2%2 else 3)}
for name,(initial,target)in cases.items():
 r=len(initial);ij=[(i,j)for i in range(r)for j in range(i,r)]
 V=np.array([v for v in product(range(-2,3),repeat=r)if any(v)and next(x for x in v if x)!=-abs(next(x for x in v if x))],dtype=np.int64)
 A=np.array([[v[i]*v[j]*(1 if i==j else 2)for i,j in ij]for v in V],dtype=float)
 b=np.array([target(v)for v in V],dtype=float)
 def matrix(x):
  m=np.zeros((r,r))
  for t,(i,j)in zip(x,ij):m[i,j]=m[j,i]=t
  return m
 best=initial.copy();bestdet=np.linalg.det(best);print('START',name,bestdet,flush=True)
 for trial in range(32):
  if trial==0:m=best.copy()
  else:
   X=rng.normal(size=(r,r));m=X@X.T+np.eye(r)
   x=np.array([m[i,j]for i,j in ij]);m*=max(1.,float(np.max(b/(A@x))))
  for it in range(18):
   inv=np.linalg.inv(m);objective=np.array([inv[i,j]*(1 if i==j else 2)for i,j in ij])
   res=linprog(objective,A_ub=-A,b_ub=-b,bounds=[(.01,12)if i==j else(-12,12)for i,j in ij],method='highs')
   if not res.success:break
   new=matrix(res.x)
   if np.linalg.eigvalsh(new)[0]<.01:break
   d=np.linalg.det(new)
   if d<bestdet-1e-7:
    best=new.copy();bestdet=d
    out={'name':name,'determinant_float':d,'gram':[[str(F(float(x)).limit_denominator(10000))for x in v]for v in best]}
    (P/(name+'_metric.json')).write_text(json.dumps(out,indent=2)+'\n');print('BEST',name,trial,it,d,out['gram'],flush=True)
   if abs(d-np.linalg.det(m))<1e-8:break
   m=new
 print('FINAL',name,bestdet,flush=True)
