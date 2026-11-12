import numpy as np,json
from pathlib import Path
P=Path(__file__).resolve().parent
# Grassl q4 n18 k10 construction intermediate [20,13,6], w=2 w²=3
B=np.array([[2,0,2,3,1,1,1],[2,0,1,2,2,2,0],[1,1,1,3,0,1,0],[3,0,0,2,2,3,1],[2,2,1,0,3,2,2],[2,0,1,1,3,0,1],[1,3,2,2,1,2,0],[2,3,1,2,1,3,1],[0,1,1,1,1,2,0],[3,2,0,3,2,2,0],[2,2,1,3,1,0,0],[0,0,2,2,3,2,1],[2,0,3,0,3,1,3]],dtype=np.uint8)
M=np.concatenate([np.eye(13,dtype=np.uint8),B],axis=1);mul=np.array([[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]],dtype=np.uint8)
r=0
for c in [18,19]:
 pivot=next(i for i in range(r,13)if M[i,c]);M[[pivot,r]]=M[[r,pivot]];iv=next(i for i in range(1,4)if mul[i,M[r,c]]==1);M[r]=mul[iv,M[r]]
 for i in range(13):
  if i!=r and M[i,c]:M[i]^=mul[M[i,c],M[r]]
 r+=1
assert not M[2:,18:].any();G=M[2:,:18];np.save(P/'outer18_G11.npy',G)
C=np.zeros((1,18),dtype=np.uint8)
for g in G:C=np.concatenate([C^mul[z,g]for z in range(4)])
w=np.count_nonzero(C,axis=1);assert w[1:].min()==6
bad=C[(w==6)|(w==7)];support=np.unique(np.sum((bad!=0).astype(np.uint32)*(1<<np.arange(18,dtype=np.uint32)),axis=1));np.save(P/'outer18_bad_supports.npy',support)
json.dump({'source':'https://www.codetables.de/BKLC/BKLC.php?q=4&n=18&k=11','G20':np.concatenate([np.eye(13,dtype=np.uint8),B],axis=1).tolist(),'G18':G.tolist(),'weight_enumerator':np.bincount(w,minlength=19).tolist(),'bad_support_count':len(support)},open(P/'outer18.json','w'),indent=2)
print('PASS size',len(C),'bad supports',len(support),'hist',np.bincount(w).tolist(),flush=True)
