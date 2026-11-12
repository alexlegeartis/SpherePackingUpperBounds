import numpy as np,json
x=np.load('research/ct36_contact_walk_20261112/flex.npz');V=x['V'];N=x['G4'];D=x['D'];VN=V@N;VD=V@D
best=[]
for a in [-1,1]:
 Q=32+2*a*VN+np.diag(N);R=2*a*VD+np.diag(D)
 for sg in [-1,1]:
  Z=Q/4+sg*.2*R;Z[Q==0]=np.inf;idx=np.unravel_index(np.argmin(Z),Z.shape);v=V[idx[0]].copy();v[idx[1]]+=a
  best.append({'t':sg*.2,'norm':float(Z[idx]),'q':int(v@N@v),'r':int(v@D@v),'v':v.tolist()})
print(json.dumps(best));print('threshold',((15/16)*(1-(12288/15625)**(1/6)))**.5)
json.dump(best,open('research/audit_ct36_flex/neighbors.json','w'))
