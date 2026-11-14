"""Cutting-plane periodic contact screen with a preconditioned FISTA dual.

Numerical discovery only.  All input rays are checked exactly for P-norm6;
the final contact set can be incomplete, and density is never certified here.
"""
import json
import pathlib
import time
import sys
import numpy as np
from scipy.linalg import null_space

p=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parent/'dim48mod3'
e=np.load(p/'antipode44_embedding.npz');P,N,w=e['P'],e['N'],e['w']
V=np.load(p/'periodic44_initial.npz')['V']
for source in sorted(p.glob('periodic44_shell*.npy')):
    V=np.vstack([V,np.load(source)])
V=np.unique(V,axis=0)
assert np.all(np.einsum('vi,ij,vj->v',V,P,V)==6)
F=V@w.T;S=np.vstack([np.zeros((1,4),int),np.eye(4,dtype=int)])
pairs={tuple(S[i]-S[j]):(i,j) for i in range(5) for j in range(5) if i!=j}
pairidx=np.array([pairs[tuple(f)] for f in F])
C=np.linalg.cholesky(P);Z=null_space(N@C);Y=V@C@Z
m=314/55;assert np.max(np.abs((Y*Y).sum(axis=1)-m))<1e-9
ii,jj=np.triu_indices(44);nmetric=len(ii);nvar=nmetric+176
goal=np.zeros(nvar);goal[:nmetric][ii==jj]=1
rng=np.random.default_rng(42)
selected=np.sort(rng.choice(len(V),min(10000,len(V)),replace=False))

def constraints(indices):
    y=Y[indices];A=np.zeros((len(y),nvar));A[:,:nmetric]=y[:,ii]*y[:,jj]*(1+(ii!=jj))
    for i in range(1,5):
        signs=(pairidx[indices,0]==i).astype(float)-(pairidx[indices,1]==i)
        A[:,nmetric+44*(i-1):nmetric+44*i]=2*y*signs[:,None]
    return A

history=[]
for stage in range(12):
    began=time.time();A=constraints(selected)
    eig,U=np.linalg.eigh(A.T@A);keep=eig>eig[-1]*1e-11
    assert keep.sum()==nvar
    X=(A@U[:,keep])/np.sqrt(eig[keep]);target=(U[:,keep].T@goal)/np.sqrt(eig[keep])
    weights=np.maximum(0,X@target);y=weights.copy();t=1.;previous=None;stable=0
    for iteration in range(5001):
        residual=X.T@y-target
        nxt=np.maximum(0,y-X@residual)
        nt=(1+np.sqrt(1+4*t*t))/2;y=nxt+(t-1)/nt*(nxt-weights);weights=nxt;t=nt
        if iteration%100==0:
            err=float(np.linalg.norm(X.T@weights-target))
            stable=stable+1 if previous is not None and abs(err-previous)<1e-11 else 0
            previous=err
            if err<1e-10 or stable>=5:break
    residual=X.T@weights-target
    h=U[:,keep]@(residual/np.sqrt(eig[keep]))
    # A near-zero residual is a eutaxy result, not a usable descent ray.
    # Do not magnify numerical noise by normalizing that residual.
    if err<1e-9:h[:]=0
    E=np.zeros((44,44));E[ii,jj]=h[:nmetric];E[jj,ii]=h[:nmetric]
    norm=np.linalg.norm(E)
    if norm>1e-15:h/=norm;E/=norm
    shifts=np.vstack([np.zeros((1,44)),h[nmetric:].reshape(4,44)])
    deriv=np.einsum('vi,ij,vj->v',Y,E,Y)+2*np.einsum('vi,vi->v',Y,shifts[pairidx[:,0]]-shifts[pairidx[:,1]])
    inflation=max(0,-float(deriv.min())/m)
    adjusted_trace=float(np.trace(E)+44*inflation)
    record={'stage':stage,'contacts_used':len(selected),'contacts_in_pool':len(V),'fista_iterations':iteration,
            'transformed_residual':err,'trace':float(np.trace(E)),
            'minimum_derivative_sample':float((A@h).min()),'minimum_derivative_pool':float(deriv.min()),
            'negative_pool_derivatives':int((deriv<-1e-8).sum()),
            'trace_after_uniform_repair':adjusted_trace,'seconds':time.time()-began}
    history.append(record);print(json.dumps(record),flush=True)
    (p/'periodic44_adaptive_report.json').write_text(json.dumps({'history':history,'scope':'sampled numerical cone only'},indent=2))
    np.savez(p/'periodic44_adaptive_state.npz',selected=selected,V=V,weights=weights,direction=h,E=E,shifts=shifts,
             Z=Z,C=C,derivatives=deriv,uniform_repair=inflation)
    if err<1e-9:print('Numerical positive eutaxy found on sample.',flush=True);break
    if adjusted_trace<-1e-7:
        print('Full collected pool admits numerical descent after uniform repair.',flush=True);break
    order=np.argsort(deriv)
    selected_set=set(selected.tolist())
    new=[int(i) for i in order if i not in selected_set and deriv[i]<-1e-9][:2000]
    if not new:print('No further sampled cuts.',flush=True);break
    selected=np.sort(np.concatenate([selected,np.array(new,dtype=np.int64)]))
