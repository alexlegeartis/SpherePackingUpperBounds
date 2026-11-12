import itertools,json
from pathlib import Path
P=Path(__file__).resolve().parent
cols=[(1,0),(0,1),(1,1),(1,2)]
def comps(n,k):
 if k==1:yield(n,);return
 for i in range(n+1):
  for t in comps(n-i,k-1):yield(i,)+t
out={}
for h in [6,7]:
 good=[]
 for counts in comps(h,4):
  sums={(0,0)}
  for c,v in zip(counts,cols):
   for _ in range(c):sums={((x+s*v[0])%3,(y+s*v[1])%3)for x,y in sums for s in(-1,1)}
  if(0,0)not in sums:good.append(counts)
 out[h]=good;print(h,good)
json.dump(out,open(P/'rank2_allowed_patterns.json','w'),indent=2)
