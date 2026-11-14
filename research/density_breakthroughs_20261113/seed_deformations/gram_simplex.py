import numpy as np,itertools,json,pathlib,sys
p=pathlib.Path(__file__).resolve().parent
k=int(sys.argv[1]) if len(sys.argv)>1 else 4
n=48-k;batch=20000;best=[];ii,jj=np.triu_indices(k,1);diags=[6] if k>4 else [6,8,10]
# Every dual simplex is represented by 0,e1,...,ek for this finite search.
# Float screen only. Exact certificate / embedding remains mandatory.
for diag in itertools.combinations_with_replacement(diags,k):
 rows=itertools.product(range(-3,4),repeat=len(ii));buf=[];tested=0
 while True:
  for r in itertools.islice(rows,batch-len(buf)):buf.append(r)
  if not buf:break
  G=np.zeros((len(buf),k,k));G[:,np.arange(k),np.arange(k)]=diag
  G[:,ii,jj]=buf;G[:,jj,ii]=buf;good=np.linalg.eigvalsh(G)[:,0]>.05;G=G[good];tested+=len(buf);buf=[]
  if not len(G):continue
  inv=np.linalg.inv(G);ds=np.diagonal(inv,axis1=1,axis2=2);beta=np.maximum(ds.max(axis=1),(ds[:,ii]+ds[:,jj]-2*inv[:,ii,jj]).max(axis=1));det=np.linalg.det(G)
  val=np.where(beta<6,(k+1)*np.maximum(0,(6-beta)/4)**(n/2)/np.sqrt(det),0)
  idx=np.argsort(val)[-5:]
  for a in idx:
   best.append({'delta':float(val[a]),'beta':float(beta[a]),'det':float(det[a]),'G':G[a].astype(int).tolist()})
  best.sort(key=lambda r:-r['delta']);best=best[:30]
 print('diag',diag,'tested',tested,'best',best[0],flush=True)
 (p/f'simplex{k}_float.json').write_text(json.dumps(best,indent=2))
