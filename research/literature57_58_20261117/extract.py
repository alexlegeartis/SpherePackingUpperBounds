"""Extract numerical LP comparators; not a certification of the prior bounds."""
import sys,re,json,hashlib
from pathlib import Path
import mpmath as m
m.mp.dps=80
if len(sys.argv)!=2: raise SystemExit('Usage: python extract.py /path/to/dataminusone.txt')
p=Path(sys.argv[1])
assert hashlib.md5(p.read_bytes()).hexdigest()=='2c062199a2305ae17f512dc00f68d735'
records={}; f=p.open()
for line in f:
 a=re.match(r'Scaling dimensions for c_total=(\S+), N=(\d+):',line)
 if a and a[1] in ['57','58']:
  d,N=int(a[1]),int(a[2]); gap=next(f).strip(); c=m.mpf(d)/2
  den=(m.pi*m.mpf(gap)/2)**c/m.gamma(c+1)
  records.setdefault(d,[]).append(dict(dimension=d,N=N,gap=gap,numerical_density=str(den)))
result=[]
for d,rs in records.items():
 r=max(rs,key=lambda r:r['N']); new=m.mpf({57:'1.5718326e-8',58:'1.073633e-8'}[d])
 r.update(certified_repo_bound=str(new),percent_below_numerical_benchmark=str(100*(1-new/m.mpf(r['numerical_density']))))
 result.append(r)
print(json.dumps(result,indent=2))
