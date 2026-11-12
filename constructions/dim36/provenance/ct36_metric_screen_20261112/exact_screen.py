import numpy as np,itertools,json
from pathlib import Path
ROOT=str(Path(__file__).resolve().parent)
from fractions import Fraction as F
P=np.load(ROOT+'/parts.npy');D=np.load(ROOT+'/types.npz');types=D['types'];reps=D['representatives'];# first19 reps, appendedtype24 lacksrep
reps=np.vstack([reps,np.load(ROOT+'/kernel_sat_basis.npy')[0],np.load(ROOT+'/newshort.npy')[0]])
Q=np.array([[sum(int(v[i])*int(M[i,j])*int(v[j])for i in range(36)for j in range(36)) for M in P]for v in reps],dtype=np.int64)
assert np.all(Q%1536==0);Q=Q//1536
rows=[[F(int(a),12)for a in r]for r in Q]
def det(A):return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
verts=set()
for ids in itertools.combinations(range(len(rows)),3):
 A=[rows[i]for i in ids];d=det(A)
 if not d:continue
 x=[]
 for j in range(3):
  M=[r.copy()for r in A]
  for r in M:r[j]=F(8)
  x.append(det(M)/d)
 if min(x)>0 and all(sum(a*b for a,b in zip(r,x))>=8 for r in rows):verts.add(tuple(x))
best=min(verts,key=lambda x:x[0]*x[1]*x[2]);prod=best[0]*best[1]*best[2];print(len(verts),best,prod,'densityupper',3/prod**6)
assert prod**12>9
json.dump({'vertices':[[str(a)for a in x]for x in sorted(verts)],'best':[str(a)for a in best],'product':str(prod),'density_ratio_upper':str(3/prod**6),'constraints':[[str(a)for a in r]for r in rows],'witness_vectors':reps.tolist(),'scope':'three scalar frame eigenspace deformation only'},open(ROOT+'/exact_screen.json','w'),indent=2)
