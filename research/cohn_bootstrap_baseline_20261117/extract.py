import re,json,csv,hashlib
from pathlib import Path
import mpmath as mp
mp.mp.dps=90;D=Path(__file__).parent;path=D/'dataminusone.txt'
assert hashlib.md5(path.read_bytes()).hexdigest()=='2c062199a2305ae17f512dc00f68d735'
allrows={}
with path.open() as f:
 for line in f:
  m=re.fullmatch(r'Scaling dimensions for c_total=(\d+), N=(\d+):\s*',line)
  if m:
   n,N=map(int,m.groups());gap=next(f).strip()
   if 1<=n<=96:allrows.setdefault(n,[]).append({'N':N,'gap':gap})
rows=[]
for n in range(49,97):
 vals=sorted(allrows[n],key=lambda v:v['N']);v=vals[-1];gap=mp.mpf(v['gap']);rho=(mp.pi*gap/2)**(mp.mpf(n)/2)/mp.gamma(mp.mpf(n)/2+1)
 prior=vals[-2]if len(vals)>1 else None
 row=dict(dimension=n,central_charge=mp.nstr(mp.mpf(n)/2,5),N=v['N'],delta1=v['gap'],density=mp.nstr(rho,45),previous_N=None if prior is None else prior['N'],relative_gap_change=None if prior is None else mp.nstr(abs(mp.mpf(prior['gap'])/gap-1),15),status='published numerical spinless-bootstrap/LP estimate; not independently certified',source='https://hdl.handle.net/1721.1/125646')
 rows.append(row)
(D/'baseline49_96.json').write_text(json.dumps(rows,indent=2));(D/'gap_history1_96.json').write_text(json.dumps(allrows,indent=2))
with (D/'baseline49_96.csv').open('w')as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys(),lineterminator='\n');w.writeheader();w.writerows(rows)
for row in rows:print(row['dimension'],row['N'],mp.nstr(mp.mpf(row['density']),14),row['relative_gap_change'])
for n in [8,24,48]:
 v=max(allrows[n],key=lambda v:v['N']);rho=(mp.pi*mp.mpf(v['gap'])/2)**(mp.mpf(n)/2)/mp.gamma(mp.mpf(n)/2+1);print('CONTROL',n,v['N'],mp.nstr(rho,30))
