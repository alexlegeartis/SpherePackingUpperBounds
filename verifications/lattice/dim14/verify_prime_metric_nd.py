if not __debug__:
 raise RuntimeError('Do not use -O for verification')
from fractions import Fraction as F
from math import isqrt
import json,sys
from pathlib import Path

n=int(sys.argv[1]);p=int(sys.argv[2]);mode=f'n{n}_p{p}';out=[]
assert n>=4 and p>=3 and p%2==1 and all(p%d for d in range(2,isqrt(p)+1)), 'Requires n>=4 and an odd prime p'
target=F(sys.argv[3]) if len(sys.argv)>3 else None
stored=json.load(open(Path(__file__).with_name(mode+'_metric.json')))
assert stored['mode']==mode
records={tuple(z['counts']):z for z in stored['cases']}
assert len(records)==len(stored['cases'])
types=list(range((p+1)//2));chars=[[F((k*i)%p if (k*i)%p<=p//2 else (k*i)%p-p,p) for i in types] for k in range(1,(p+1)//2)]
def compositions(total,k):
 if k==1:yield[total];return
 for x in range(total+1):
  for tail in compositions(total-x,k-1):yield[x]+tail
for case,counts in enumerate(compositions(n,len(types))):
 if sum(counts[1:])==0:continue
 groups=[i for i in range(len(types)) if counts[i]];keys=[('a',g) for g in groups]+[('c',g,h) for g in groups for h in groups if g<=h and(g!=h or counts[g]>=2)];idx={key:i for i,key in enumerate(keys)};dim=len(keys)
 rows=[];rhs=[]
 def blank():return[F(0)]*dim
 for g in groups:
  r=blank();r[idx['a',g]]=-1;rows.append(r);rhs.append(F(1))
 for key in keys:
  if key[0]!='c':continue
  _,g,h=key
  for sg in [1,-1]:
   r=blank();r[idx['a',g]]-=F(1,2);r[idx['a',h]]-=F(1,2);r[idx[key]]=sg;rows.append(r);rhs.append(F(0))
 groupids=[g for g in groups for _ in range(counts[g])]
 for rr in chars:
  base=[rr[g] for g in groupids];par=sum(int(p*x) for x in base)%2;variants=[]
  if not par:variants.append(base)
  else:
   for g in groups:
    ix=groupids.index(g)
    for sg in([1,-1] if rr[g]==0 else[1 if rr[g]>0 else -1]):
     v=base[:];v[ix]-=sg;variants.append(v)
  for v in variants:
   r=blank();sums={};sq={}
   for g in groups:
    vv=[x for x,gg in zip(v,groupids) if gg==g];sums[g]=sum(vv);sq[g]=sum(x*x for x in vv);r[idx['a',g]]=-sq[g]
   for key in keys:
    if key[0]!='c':continue
    _,g,h=key;r[idx[key]]=-(sums[g]*sums[h]*(2 if g!=h else 1)-(sq[g] if g==h else 0))
   rows.append(r);rhs.append(2*sum(sq.values())-4)
 obj=[F(counts[key[1]]) if key[0]=='a' else F(0) for key in keys]
 rec=records.pop(tuple(counts));assert rec['keys']==[list(k) for k in keys]
 y=list(map(F,rec['dual']));assert len(y)==len(rhs)
 assert all(v<=0 for v in y)
 assert all(sum(row[j]*v for row,v in zip(rows,y))==obj[j] for j in range(dim)),counts
 bound=sum(b*v for b,v in zip(rhs,y))
 if target is not None:assert bound>target,(counts,bound)
 out.append({'counts':counts,'keys':keys,'dual':list(map(str,y)),'bound':str(bound)})
 if len(out)%500==0:print(mode,len(out),'min',min(F(z['bound']) for z in out),flush=True)
assert not records
print('ALLPASS',mode,len(out),'min',min((F(z['bound']),z['counts']) for z in out),flush=True)
json.dump({'mode':mode,'passed':True,'case_count':len(out),'minimum_bound':str(min(F(z['bound']) for z in out)),'target':str(target)},open(Path(__file__).with_name(mode+'_metric_replay.json'),'w'),indent=2)
