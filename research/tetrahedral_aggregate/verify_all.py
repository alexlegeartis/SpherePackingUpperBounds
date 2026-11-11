"""Portable complete replay. No density bound is asserted."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,shutil,subprocess,sys,tempfile,time
if not __debug__:raise RuntimeError('Assertions are required; do not use Python -O')
D=Path(__file__).resolve().parent
for line in (D/'SHA256SUMS').read_text().splitlines():
 expected,name=line.split('  ',1)
 assert hashlib.sha256((D/name).read_bytes()).hexdigest()==expected, name
sources=json.loads((D/'primary/tetra_full_geometry_screen.json').read_text())
for q,numeric,primary,independent,upper in [
 ('21/20',1.05,'tetra_full_geometry22_check.json','independent_q21_20','3341.14916082'),
 ('11/10',1.1,'tetra_full_geometry22_q11_10_monotone.json','independent_q11_10','5086.98083741')]:
 p=json.loads((D/'primary'/primary).read_text());i=json.loads((D/independent/'certificate.json').read_text())
 for key in ['dimension','q','coefficients','polynomial','B_exact']:assert p[key]==i[key],key
 raw=next(r for r in sources if r['dimension']==22 and r['q']==numeric)
 assert [F(str(c))for c in raw['coefficients']]==list(map(F,p['coefficients']))
 assert F(p['B_exact'])<F(upper)
 q=F(q);lo=1-q/(2*(1-q/4));hi=1-1/(2*(q-F(1,4)))
 def evaluate(t):
  result=F(0)
  for c in reversed(p['polynomial']):result=result*t+F(c)
  return result
 assert lo<F(1,3)<hi
 assert evaluate(hi)/(hi-F(1,3))>evaluate(lo)/(lo-F(1,3))
print('Hashes, exact witness agreement, outward rounding and both affine contradictions PASS',flush=True)
start=time.time()
with tempfile.TemporaryDirectory(prefix='tetrahedral_replay_')as td:
 work=Path(td)/'package';shutil.copytree(D,work)
 commands=[
  ['primary/verify_tetra_envelopes.py'],
  ['primary/verify_tetra_separation.py'],
  ['primary/verify_tetra_full_geometry.py','--dimension','22','--q','21/20','--seconds','120','--output',str(Path(td)/'primary_q21.json')],
  ['primary/verify_tetra_full_geometry_monotone.py','--dimension','22','--q','11/10','--seconds','120','--output',str(Path(td)/'primary_q11.json')],
  ['independent_q21_20/verify.py'],
  ['independent_q21_20/verify_vertex.py'],
  ['independent_q11_10/verify.py']]
 for args in commands:
  print('RUN',args[0],flush=True)
  subprocess.run([sys.executable,'-B',*args],cwd=work,check=True)
print('ALL COMPLETE REPLAYS PASS in',round(time.time()-start,2),'seconds. Auxiliary inequalities only; no density claim.',flush=True)
