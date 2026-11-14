import numpy as np,pathlib,sys,json
from scipy.optimize import linprog,minimize,LinearConstraint
p=pathlib.Path(__file__).resolve().parent/'dim48mod3';f=np.load(p/'flexes.npz');H=f['H'];D=f['Ds']/206.;P=np.load(p/'P.npy');V=np.vstack([f['V'],-f['V']]);xstar=np.array([-2,4,6,4,2,28])/3;raw=np.load(p/'reduced_cuts.npz')['raw'];rng=np.random.default_rng(29)
Gtarget=np.array([[6,-1,-1,-2],[-1,6,-2,-2],[-1,-2,6,-1],[-2,-2,-1,6]])
found=None
for ia in rng.choice(len(V),100,replace=False):
 dots=V@P@V[ia];Bs=np.flatnonzero(dots==-1);rng.shuffle(Bs)
 for ib in Bs[:50]:
  dotb=V@P@V[ib];Cs=np.flatnonzero((dots==-1)&(dotb==-2));Ds=np.flatnonzero((dots==-2)&(dotb==-2))
  if not len(Cs)or not len(Ds):continue
  dc=V[Cs]@P@V[Ds].T;ii,jj=np.where(dc==-1)
  if len(ii):found=V[[ia,ib,Cs[ii[0]],Ds[jj[0]]]];break
 if found is not None:break
assert found is not None;assert np.array_equal(found@P@found.T,Gtarget);w=found@P;np.savez(p/'antipode44_embedding.npz',N=found,w=w,P=P,Gtarget=Gtarget)
print('embedding',found.tolist(),flush=True)

def calc(x):
 Q=H+np.einsum('i,ijk->jk',x,D);ev=np.linalg.eigvalsh(Q)
 if ev[0]<1e-10:return -1e8,None
 Qi=np.linalg.inv(Q);M=w@Qi@w.T;B=np.linalg.inv(M);ds=B.diagonal();ij=np.triu_indices(4,1);beta=max(ds.max(),(ds[ij[0]]+ds[ij[1]]-2*B[ij]).max());delta=5*((8-beta)/4)**22/np.exp((np.linalg.slogdet(Q)[1]+np.linalg.slogdet(M)[1])/2)
 return delta,beta
for t in [0.,.25,.5,.75,.9,.95,.99,1.]:print('segment',t,calc(t*xstar),flush=True)
base=calc(xstar)[0]
vertices=np.load(p/'polytope_vertices.npz')['positive_vertices'];scores=[calc(x)[0]for x in vertices];bestidx=np.argmax(scores);print('vertexbest',scores[bestidx],vertices[bestidx],flush=True)
con=LinearConstraint(raw[:,1:],raw[:,0],np.inf)
res=minimize(lambda x:-np.log(calc(x)[0]) if calc(x)[0]>0 else 1e8,xstar,method='SLSQP',constraints=[con],options={'maxiter':200,'ftol':1e-12})
print('opt',res.success,res.message,calc(res.x),res.x,flush=True)
(p/'antipode44_metric_report.json').write_text(json.dumps({'target_G':Gtarget.tolist(),'baseline':base,'segment':[(t,calc(t*xstar))for t in [0,.5,.9,.99,1]],'best_vertex':vertices[bestidx].tolist(),'best_vertex_value':scores[bestidx],'optimizer':{'success':bool(res.success),'message':res.message,'x':res.x.tolist(),'density_conditional':calc(res.x)[0],'beta':calc(res.x)[1]},'scope':'numeric necessary-constraint relaxation; no verified new minimum'},indent=2))
