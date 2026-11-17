"""Exact nearest-distance to 2E8, using integer and half-integer D8 cosets."""
import numpy as np,itertools,json
from pathlib import Path
P=Path(__file__).resolve().parent
E=np.zeros((8,8),dtype=np.int64)
for i in range(6):E[i,i]=2;E[i,i+1]=-2
E[6,5:7]=2;E[7,:]=1
G=E@E.T//4;assert np.array_equal(E@E.T,4*G)
def mod2rank(A):
 A=np.array(A,dtype=np.int64)%2;row=0
 for c in range(A.shape[1]):
  ids=np.flatnonzero(A[row:,c])
  if not len(ids):continue
  j=row+ids[0];A[[row,j]]=A[[j,row]]
  for j in range(A.shape[0]):
   if j!=row and A[j,c]:A[j]^=A[row]
  row+=1
 return row
# Hyperbolic basis in E8/2E8, with all vectors singular.
V=np.array(list(itertools.product(range(2),repeat=8)),dtype=np.int64);q=np.einsum('bi,ij,bj->b',V,G,V)//2%2;es=[];fs=[]
for k in range(4):
 cand=[v for v,qq in zip(V,q)if qq==0 and any(v)and all(int(v@G@w)%2==0 for w in es+fs)]
 e=next(v for v in cand if mod2rank(es+fs+[v])==2*k+1);f=next(v for v in cand if int(v@G@e)%2==1);es.append(e);fs.append(f)
H=np.array(es+fs);assert mod2rank(H)==8
states=np.array([np.array([(s>>i)&1 for i in range(8)])@H for s in range(256)],dtype=np.int64);S2=states@E
# Input numerators have denominator16 after dividing the 2E8 problem by2.
def nearest(num):
 best=None
 for off in [0,8]:
  a=num-off;z=np.floor_divide(a+8,16);res=a-16*z;cost=int(res@res)
  if int(z.sum())%2:cost+=min(256-32*abs(int(x))for x in res)
  best=cost if best is None else min(best,cost)
 return best
if __name__=='__main__':
 from flint import fmpz_mat as Z
 assert Z(G.tolist()).det()==1
 rng=np.random.default_rng(2026111780);rows=[]
 for j in range(10):
  # x=quarter of an E8 vector, Euclidean denominator8.
  coeff=rng.integers(0,8,8);target=coeff@E
  costs=[nearest(target-4*s)for s in S2];rows.append({'target_num_over8':target.tolist(),'costs_numerator_over128':costs})
  print(j,{v:costs.count(v)for v in sorted(set(costs))},flush=True)
 (P/'local_tables.json').write_text(json.dumps({'E8_basis_times2':E.tolist(),'E8_gram':G.tolist(),'polarization_basis':H.tolist(),'residue_representatives':states.tolist(),'rows':rows}));
