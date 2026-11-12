"""Discovery-only exact-numerator CVP scores with fplll floating enumeration.
A high score requires a separate rigorous CVP certificate before any density claim.
"""
import argparse,json,time,math
from pathlib import Path
import numpy as np
from fpylll import Enumeration,GSO,IntegerMatrix
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=13);ap.add_argument('--denominator',type=int,default=3);ap.add_argument('--restarts',type=int,default=20);ap.add_argument('--seconds',type=float,default=300);ap.add_argument('--random',type=int,default=0);args=ap.parse_args()
Q=np.loadtxt(ROOT/'constructions/dim36/Q_reduced.txt',dtype=np.int64);n=len(Q);M=GSO.Mat(IntegerMatrix.from_matrix(Q.tolist()),gram=True,update=True)
mu=np.array([[M.get_mu(j,i) if j>i else int(j==i) for i in range(n)]for j in range(n)])
p=args.denominator; rng=np.random.default_rng(args.seed);start=time.time();cache={};calls=0;best=-1;bestrow=None;trace=[]
def score(a):
 global calls,best,bestrow
 key=bytes(a)
 if key in cache:return cache[key][0]
 s=a/p;target=s@mu;E=Enumeration(M);r=E.enumerate(0,n,200.,0,target=target.tolist());z=np.rint(r[0][1]).astype(np.int64);u=a-p*z;num=int(u@Q@u);assert abs(num/(p*p)-r[0][0])<1e-6
 calls+=1;cache[key]=(num,z)
 if num>best:
  best=num;bestrow={'shift_numerators':a.tolist(),'denominator':p,'Q_divisor':20,'nearest_witness':z.tolist(),'distance_squared_num':num,'distance_squared_den':20*p*p,'squared_distance':num/(20*p*p),'calls':calls,'seconds':time.time()-start,'scope':'CVP result uses floating fplll enumeration; not a rigorous minimum certificate'}
  trace.append(bestrow);(P/f'torsion_p{p}_seed{args.seed}.json').write_text(json.dumps({'best':bestrow,'trace':trace,'calls':calls},indent=2));print(json.dumps(bestrow),flush=True)
 return num
for restart in range(args.restarts):
 if time.time()-start>args.seconds:break
 a=rng.integers(0,p,n,dtype=np.uint8);sc=score(a)
 if args.random:
  for i in range(args.random):
   if time.time()-start>args.seconds:break
   score(rng.integers(0,p,n,dtype=np.uint8))
  break
 tabu=set();flat=0
 for step in range(300):
  if time.time()-start>args.seconds:break
  candidates=[]
  for j in rng.permutation(n):
   for d in range(1,p):
    b=a.copy();b[j]=(int(b[j])+d)%p
    if bytes(b) in tabu:continue
    q=score(b);candidates.append((q,rng.random(),b))
    if time.time()-start>args.seconds:break
   if time.time()-start>args.seconds:break
  if not candidates:break
  q,_,b=max(candidates,key=lambda r:r[:2]);tabu.add(bytes(a))
  if q<sc:
   if flat>10:break
   # Allow deterministic-frequency downhill moves to cross short plateaus.
   if rng.random()<.1:a=b;sc=q;flat+=1
   else:break
  else:
   flat=flat+1 if q==sc else 0;a=b;sc=q
  if flat>25:break
 print('restart',restart,'calls',calls,'best',best/(20*p*p),'elapsed',time.time()-start,flush=True)
print('DONE',calls,'best',best/(20*p*p),flush=True)
