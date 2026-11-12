exec(open('research/audit_ct36_flex/neighbors.py').read().split('best=[]')[0])
best=(1e99,None)
for i in range(36):
 for j in range(i+1,36):
  for a in [-1,1]:
   for b in [-1,1]:
    Q=32+2*a*VN[:,i]+2*b*VN[:,j]+N[i,i]+N[j,j]+2*a*b*N[i,j]
    R=2*a*VD[:,i]+2*b*VD[:,j]+D[i,i]+D[j,j]+2*a*b*D[i,j]
    mask=(Q>0)&(R!=0);inds=np.where(mask)[0]
    if not len(inds):continue
    ts=(Q[mask]/4-8)/abs(R[mask]);k=inds[np.argmin(ts)];t=float(ts.min())
    if t<best[0]:
     v=V[k].copy();v[i]+=a;v[j]+=b;best=(t,{'q':int(v@N@v),'r':int(v@D@v),'v':v.tolist()})
print(best);json.dump(best,open('research/audit_ct36_flex/pairs.json','w'))
