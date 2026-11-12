import numpy as np,sympy as sp,json,math,time
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parent;J=np.load(P/'spectator_reduced.npz');G=J['G'];L,D=sp.Matrix(G).LDLdecomposition(hermitian=False);L=[[F(L[i,j])for j in range(12)]for i in range(12)];d=[F(D[i,i])for i in range(12)];assert all(x>0 for x in d);z=[0]*12;mins=[999]*4096;wit=[None]*4096;nodes=0
start=time.time()
def go(k,q):
 global nodes
 nodes+=1
 if k<0:
  assert q.denominator==1;mask=sum((v&1)<<i for i,v in enumerate(z))
  if q.numerator<mins[mask]:mins[mask]=q.numerator;wit[mask]=z.copy()
  return
 c=-sum((L[j][k]*z[j]for j in range(k+1,12)),F(0));r=(16-q)/d[k]
 if r<0:return
 b=c.denominator;a=c.numerator;rad=math.isqrt((r.numerator*b*b)//r.denominator);lo=-((-(a-rad))//b);hi=(a+rad)//b
 for x in range(lo,hi+1):z[k]=x;go(k-1,q+d[k]*(x-c)**2)
go(11,F(0));assert max(mins)<=16
old=np.loadtxt(P/'spectator_witnesses.txt',dtype=np.int64);assert list(old[:,1])==mins
labels=J['labels'];idx=np.sum(labels*(1<<np.arange(12)),axis=1);gluemins=np.array(mins)[idx];np.savez(P/'spectator_certified.npz',norm_numerators=gluemins,denominator=np.array(4),witnesses=np.array(wit)[idx],labels=labels,G=G,U=J['U'])
json.dump({'exact':True,'bound':16,'nodes':nodes,'seconds':time.time()-start,'histogram':{str(m):mins.count(m)for m in sorted(set(mins))},'scope':'all4096 half-integral spectatorcosets;normnumeratorsdividedby4'},open(P/'spectator_exact.json','w'),indent=2);print('EXACT PASS',nodes,time.time()-start,sorted(set(mins)),flush=True)
