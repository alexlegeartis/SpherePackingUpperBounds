if not __debug__:
 raise RuntimeError("Do not run with python -O")
"""Exact certificate replay: no LP solver and no floating-point arithmetic."""
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
out=[]
for p in [13,17,19,23,29,31,37,41,43]:
 h=(p-1)//2;R=[[min(i*k%p,p-i*k%p) for i in range(1,h+1)] for k in range(1,h+1)]
 # Regenerate all parity supports and certify the complements by exact average.
 pop=np.zeros(1,dtype=np.int16)
 for j in range(h):pop=np.concatenate([pop,pop+1])
 corr=np.zeros(1<<h,dtype=np.int64)
 for k in range(h):
  par=np.zeros(1,dtype=np.int8);mx=np.zeros(1,dtype=np.int16)
  for j in range(h):
   r=R[k][j];par=np.concatenate([par,par^(r%2)]);mx=np.concatenate([mx,np.maximum(mx,r)])
  corr+=par.astype(np.int64)*(p-2*mx.astype(np.int64))
 m=16-pop%2;score=m.astype(np.int64)*(p+1)*h+12*corr
 bad=set(map(int,np.where((pop<=m)&(score>18*p*h))[0]))
 data=json.load(open(Path(__file__).with_name(f'prime_support_branch_p{p}.json')))
 assert data['p']==p and set(r['mask'] for r in data['records'])==bad
 total=0;leafmax=F(0)
 for rec in data['records']:
  assert rec['bad'] is None
  mask=rec['mask'];bits=[mask>>j&1 for j in range(h)];budget=(16-sum(bits))//2;nodes=rec['nodes'];visited=set()
  def walk(index,expected_lo,expected_up):
   nonlocal_dummy=None
   assert index not in visited;visited.add(index)
   nd=nodes[index];lo=nd['lo'];up=nd['up'];assert lo==expected_lo and up==expected_up
   remaining=budget-sum(lo)
   if nd.get('infeasible'):
    assert remaining<0 or any(a>b for a,b in zip(lo,up));return []
   assert remaining>=0 and all(0<=a<=b<=budget for a,b in zip(lo,up))
   if 'branch_var' in nd:
    j=nd['branch_var'];cut=nd['branch_cut'];u=up.copy();u[j]=min(u[j],cut);l=lo.copy();l[j]=max(l[j],cut+1)
    return walk(nd['left'],lo,u)+walk(nd['right'],l,up)
   active=[bool(e or z) for e,z in zip(bits,lo)];cnt=[e+2*z for e,z in zip(bits,lo)]
   if 'alpha' in nd:
    aa=list(map(F,nd['alpha']));assert len(aa)==h and min(aa)>=0 and sum(aa)==1
    bs=[]
    for k in range(h):
     parity=sum(R[k][j]*bits[j] for j in range(h))%2;mx=max(R[k][j] for j in range(h) if active[j])
     bs.append(F(sum(R[k][j]**2*cnt[j] for j in range(h)),p*p)+parity*(1-F(2*mx,p)))
    weight=[sum((aa[k]*F(2*R[k][j]**2,p*p) for k in range(h)),F(0)) for j in range(h)]
    bound=sum((a*b for a,b in zip(aa,bs)),F(0));left=remaining
    for j in sorted(range(h),key=lambda j:weight[j],reverse=True):
     take=min(left,up[j]-lo[j]);bound+=take*weight[j];left-=take
    assert bound==F(nd['bound']) and bound<=F(3,2);return [bound]
   assert 'point_counts' in nd and lo==up and cnt==nd['point_counts']
   vals=[]
   for row in R:
    par=sum(a*b for a,b in zip(row,cnt))%2;mx=max(r for r,z in zip(row,cnt) if z)
    vals.append(F(sum(r*r*z for r,z in zip(row,cnt)),p*p)+par*(1-F(2*mx,p)))
   bound=min(vals);assert bound==F(nd['point_bound']) and bound<=F(3,2);return [bound]
  bounds=walk(0,[0]*h,[budget]*h);assert len(visited)==len(nodes)
  total+=len(nodes);leafmax=max(leafmax,max(bounds,default=F(0)))
 out.append(dict(p=p,bad_supports=len(bad),tree_nodes=total,max_leaf_bound=str(leafmax),all_parity_supports=1<<h,passed=True))
 print(p,'PASSED',len(bad),'supports',total,'nodes','maxleaf',float(leafmax),flush=True)
Path(__file__.replace('.py','.json')).write_text(json.dumps(out,indent=2))
print('ALL p13..43 WITNESSES RAW <=3/2 CERTIFIED',flush=True)
