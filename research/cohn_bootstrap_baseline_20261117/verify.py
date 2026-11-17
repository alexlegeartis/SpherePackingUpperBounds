from pathlib import Path
import json,csv
import mpmath as mp
D=Path(__file__).resolve().parent
rows=json.loads((D/'baseline49_96.json').read_text())
history=json.loads((D/'gap_history1_96.json').read_text())
assert [r['dimension'] for r in rows]==list(range(49,97))
csvrows=list(csv.DictReader((D/'baseline49_96.csv').open()))
assert len(csvrows)==48
for precision in (90,120):
 mp.mp.dps=precision
 for r,c in zip(rows,csvrows):
  n=r['dimension']; vals=sorted(history[str(n)],key=lambda v:v['N']);v=vals[-1]
  assert r['N']==v['N'] and r['delta1']==v['gap']
  assert mp.mpf(r['central_charge'])==mp.mpf(n)/2
  rho=(mp.pi*mp.mpf(v['gap'])/2)**(mp.mpf(n)/2)/mp.gamma(mp.mpf(n)/2+1)
  assert mp.nstr(rho,45)==r['density']
  for key,value in r.items(): assert c[key]==('' if value is None else str(value)),key
  assert r['previous_N']==vals[-2]['N']
  change=abs(mp.mpf(vals[-2]['gap'])/mp.mpf(v['gap'])-1)
  assert abs(mp.mpf(r['relative_gap_change'])-change)<mp.mpf('1e-14')*max(change,mp.mpf('1e-80'))
 for n,target in [(8,mp.pi**4/384),(24,mp.pi**12/mp.factorial(12)),(48,mp.mpf('4.529067790473726e-7'))]:
  v=max(history[str(n)],key=lambda v:v['N']);gap=mp.mpf(v['gap'])
  rho=(mp.pi*gap/2)**(mp.mpf(n)/2)/mp.gamma(mp.mpf(n)/2+1)
  assert abs(rho/target-1)<mp.mpf('1e-14'),(n,rho,target)
text=(D.parent.parent/'NUMERICAL_BENCHMARKS49_96.md').read_text()
expected=['| %d | %s | %d |'%(r['dimension'],mp.nstr(mp.mpf(r['density']),14),r['N']) for r in rows]
actual=[line for line in text.splitlines() if line.startswith('| ') and line.split('|')[1].strip().isdigit()]
assert actual==expected
print('PASS: 48 highest-N rows, CSV/JSON, 90/120-digit conversions, controls and displayed table. Numerical benchmarks only; no bound certification.')
