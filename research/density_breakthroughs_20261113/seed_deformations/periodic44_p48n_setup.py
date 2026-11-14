"""Embed the det605 normal lattice directly in the published P48n Gram."""
import json
import pathlib
import numpy as np
from prepare_sources import parse

root=pathlib.Path(__file__).resolve().parent;p=root/'periodic44_p48n';p.mkdir(exist_ok=True)
d=np.load(root.parent/'gluing'/'P48n_glue_pool.npz');P,V=d['G'],d['V']
catalogue,_=parse('CQ48a');assert np.array_equal(P,catalogue)
rng=np.random.default_rng(195)
trial=V[rng.choice(len(V),min(16000,len(V)),replace=False)]
trial=np.vstack([trial,-trial]);TP=trial@P
target=np.array([[6,-1,-1,-2],[-1,6,-2,-2],[-1,-2,6,-1],[-2,-2,-1,6]],np.int64)
found=None
for ia in rng.choice(len(trial),100,replace=False):
    dots=TP@trial[ia];Bs=np.flatnonzero(dots==-1);rng.shuffle(Bs)
    for ib in Bs[:50]:
        dotb=TP@trial[ib];Cs=np.flatnonzero((dots==-1)&(dotb==-2));Ds=np.flatnonzero((dots==-2)&(dotb==-2))
        if not len(Cs) or not len(Ds):continue
        inner=TP[Cs]@trial[Ds].T;ii,jj=np.where(inner==-1)
        if len(ii):found=trial[[ia,ib,Cs[ii[0]],Ds[jj[0]]]];break
    if found is not None:break
assert found is not None and np.array_equal(found@P@found.T,target)
w=found@P;F=V@w.T
S=np.vstack([np.zeros((1,4),int),np.eye(4,dtype=int)])
allowed={tuple(S[i]-S[j]) for i,j in [(2,0),(2,1),(3,1),(4,0),(4,3),(0,2),(1,2),(1,3),(0,4),(3,4)]}
keep=np.array([tuple(f) in allowed for f in F]);active=V[keep]
assert np.all(np.einsum('vi,ij,vj->v',active,P,active)==6)
np.savez(p/'antipode44_embedding.npz',P=P,N=found,w=w,Gtarget=target,source_name=np.array('CQ48a'))
np.savez_compressed(p/'periodic44_initial.npz',V=active)
(p/'setup.json').write_text(json.dumps({'source':'https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ48a.html','mother_name':'P48n','published_minimum':6,'catalogue_gram_identity_exact':True,'source_vectors':len(V),'active_contact_vectors':len(active),'normal_gram':target.tolist()},indent=2))
print('Published P48n embedding found; active vectors',len(active),flush=True)
