"""Exact algebra of the candidate six-dimensional 48D deformation pencil."""
from pathlib import Path
import numpy as np,sympy as s,json,time
P=Path(__file__).resolve().parent
a=np.load(P.parent/'seed_deformations/dim48mod3/flexes.npz')
H=s.Matrix(a['H']);Ds=[s.Matrix(x)for x in a['Ds']]
start=time.time();Hi=H.inv();As=[Hi*x for x in Ds]
print('inverse',time.time()-start,flush=True)
vectors=[s.eye(48)[:,0],s.Matrix(a['V'][0]),s.Matrix(a['V'][101])]
records=[]
for t,v in enumerate(vectors):
 C=s.Matrix.hstack(v,*[A*v for A in As]);piv=C.rref()[1];B=C[:,list(piv)]
 for iteration in range(3):
  big=s.Matrix.hstack(B,*[A*B for A in As]);pp=big.rref()[1]
  if len(pp)==B.cols:break
  B=big[:,list(pp)]
 print('cyclic_dimension',t,B.cols,flush=True)
 records.append({'vector':t,'cyclic_dimension':B.cols})
 if B.cols!=4:continue
 rows=list(B.T.rref()[1]);inv=B[rows,:].inv();Rs=[inv*(A*B)[rows,:]for A in As]
 assert all(A*B==B*R for A,R in zip(As,Rs))
 E=B.T*H*B
 out={'B':[[str(x)for x in row]for row in B.tolist()],'metric':[[str(x)for x in row]for row in E.tolist()],'representations':[[[str(x)for x in row]for row in R.tolist()]for R in Rs]}
 (P/f'reduced48_{t}.json').write_text(json.dumps(out,indent=2)+'\n')
 print('reduced representations saved',t,flush=True)
 break
(P/'algebra48.json').write_text(json.dumps({'probes':records,'elapsed':time.time()-start},indent=2)+'\n')
